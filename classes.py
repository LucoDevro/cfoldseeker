#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import time
import requests
import kegg_pull.pull as kgp
import polars as pl
from pathlib import Path
from itertools import chain
from copy import deepcopy
from concurrent.futures import ThreadPoolExecutor


class Hit:
    
    def __init__(self, uniprot, query, kegg_id = '', name = '', scaff = '', start = 0, end = 0, strand = ''):
        
        self.uniprot: str = uniprot #Uniprot ID
        self.query: str = query #ID of the homologous query protein
        self.kegg_id: list = kegg_id #KEGG ID
        self.name: list = name #name of the hit in the KEGG entry
        self.scaff: str = scaff #scaffold this hit is encoded in
        self.start: int = int(start) #starting position of the coding gene on the scaffold
        self.end: int = int(end) #end position of the coding gene on the scaffold
        self.strand: str = strand #DNA strand the encoding gene is part from
    
    def __repr__(self):
        return f"{self.query} Hit {self.uniprot}\t{self.scaff} {self.start}-{self.end} ({self.strand})"
        

class Cluster:
    
    def __init__(self, hits, number, score = 0):
        
        self.hits: list(Hit) = hits
        self.number: int = number
        self.score: int = score
        
        self.start: int = min([h.start for h in self.hits])
        self.end: int = max([h.end for h in self.hits])
        
        self.strand: str = list(set([h.strand for h in self.hits]))[0]
        if self.strand == "-":
            self.start, self.end = self.end, self.start
        
        self.length: int = abs(self.end - self.start)
        
        self.scaff: str = list(set([h.scaff for h in self.hits]))[0]
        
        return None
    
    def __repr__(self):
        
        return f"Cluster {self.number}: {len(self.hits)} proteins from {self.scaff} ({self.start} - {self.end}), ({self.strand})\tScore: {self.score}"
    
    
class Search:
    
    def __init__(self, query, params, mapping_table_path, hits = {}, clusters = []):
        
        self.query: dict = query # dictionary of query names as keys and file paths as values
        self.params: dict = params
        self.hits: list = hits
        self.clusters: list = clusters
        # The full Uniprot crossref mapping table casted in a polars DataFrame
        self.mapping_table: dict = pl.read_csv("uniprot_kegg_nucl.gz", has_header = False, separator = "\t",
                                               new_columns = ['Uniprot_ID', 'Target_DB', 'ID'])
                
        return None
    
    def __repr__(self):
        pass
    
    def run_foldseek(self):

        FOLDSEEK_SUBMISSION_URL = "https://search.foldseek.com/api/ticket"
        FOLDSEEK_RESULTS_URL = "https://search.foldseek.com/api/result"
        
        """
        Submits one structure file to the FoldSeek API and returns the submission ticket in dictionary form.
        """
        def submit_foldseek_query(query_path):
            with open(query_path, "rb") as f:
                files = {"q": f}
                data = [
                    ("mode", "3diaa"),
                    ("database[]", "afdb-proteome"),
                    ("database[]", "afdb-swissprot")
                ]
                response = requests.post(FOLDSEEK_SUBMISSION_URL, files=files, data=data)
                if response.status_code == 200:
                    return response.json()
                else:
                    raise Exception(f"Error submitting query {query_path}!")
                    
            return None
        
        """
        Checks the status of the job associated with the provided job ID. Returns the status flag.
        """
        def check_query_status(job_id):
            url = f"{FOLDSEEK_SUBMISSION_URL}/{job_id}"
            results = requests.get(url).json()
            status = results['status']
            
            return status
        
        """
        Waits for the FoldSeek job associated with the provided job ID to complete, and retrieves its result.
        Returns the parsed results in a dictionary.
        """
        def retrieve_foldseek_results(job_id):
            while True:
                status = check_query_status(job_id)
                if status == "COMPLETE":
                    entry = 0
                    url = f"{FOLDSEEK_RESULTS_URL}/{job_id}/{entry}"
                    results = requests.get(url).json()
                    break
                else:
                    time.sleep(2)
                    
            return results
            
        
        # First submit all query proteins to FoldSeek.
        with ThreadPoolExecutor(max_workers = 5) as executor:
            tickets = dict(zip(self.query.keys(), executor.map(submit_foldseek_query, self.query.values())))
                
        # Then wait for the results and retrieve them automatically when completed
        all_job_ids = [ticket['id'] for ticket in tickets.values()]
        with ThreadPoolExecutor(max_workers = 5) as executor:
            all_results = dict(zip(self.query.keys(), executor.map(retrieve_foldseek_results, all_job_ids)))
        
        return all_results
            
                
    def crossref(self):
        
        """
        Extracts the entry name from a full pulled KEGG record.
        """
        def extract_entry_name(entry: str) -> str:
            return [line for line in entry.split('\n') if 'NAME' in line][0][12:]
        
        """
        Extracts the genomic positional information from a full pulled KEGG record.
        """
        def extract_positional_information(entry: str) -> str:
            return [line for line in entry.split('\n') if 'POSITION' in line][0][12:]
        
        """
        Extracts a mapping dictionary from the full Uniprot crossref mapping table, linking the Uniprot ID (keys) with a list of
        one or more multiple crossrefs in the provided target database (values).
        """
        def prepare_mapping_dict(crossref_df: pl.DataFrame, target_db: str) -> dict:
            res = crossref_df.filter(crossref_df['Target_DB'] == target_db).drop('Target_DB')
            res = res.group_by('Uniprot_ID').all()
            res = dict(zip(res['Uniprot_ID'], res['ID'].to_list()))
            
            return res
        
        """
        Auxiliary function that splits a hit with multiple crossrefs for a given attribute into multiple hits that are identical
        except for that exact crossref attribute. The split hits are first introduced at the index of the old hit before the full
        hit list is being flattened.
        
        Mutates:
            self.hits: list: the hit list containing Hit objects that represent proteins that are structurally homologous to
                             one of the query proteins.
        """
        def split_multiple_crossrefs(self, attr):
            # Loop over all hits, but keep track of the index
            for idx, h in enumerate(self.hits):
                # Determine whether we should split a record
                nb_copies = len(getattr(h, attr))
                # If we need to split a record, make as much copies as there are crossrefs and adjust each crossref to one of these
                if nb_copies > 1:
                    split_records = []
                    copy_counter = 0
                    while copy_counter < nb_copies:
                        new_record = deepcopy(h)
                        setattr(new_record, attr, getattr(h, attr)[copy_counter])
                        split_records.append(new_record)
                        copy_counter += 1
                    # Insert at the index of the old record
                    self.hits[idx] = split_records
                # If we don't need to split the record, just unlist it
                else:
                    setattr(h, attr, getattr(h, attr)[0])
                    # Insert at the index of the old record
                    self.hits[idx] = [h]
            
            # Flatten the full hit list again
            self.hits = list(chain(*self.hits))
            
            return None
        
        
        ## First, fill all KEGG IDs
        # Extract a mapping table for KEGG IDs from the Uniprot crossref mapping table
        all_uniprot_ids = [h.uniprot for h in self.hits]
        all_cross_refs = self.mapping_table.filter(self.mapping_table['Uniprot_ID'].is_in(all_uniprot_ids))
        all_uniprot_kegg = prepare_mapping_dict(all_cross_refs, 'KEGG')
        # Fill the KEGG IDs
        for h in self.hits:
            h.kegg_id = all_uniprot_kegg[h.uniprot]
        # Split records with multiple crossrefs
        split_multiple_crossrefs(self, 'kegg_id')
            
        ## Then, pull all KEGG records
        all_kegg_ids = [h.kegg_id for h in self.hits]
        pull = kgp.MultiProcessMultiplePull(n_workers = 2)
        _, records = pull.pull_dict(all_kegg_ids)
        
        ## Extract the record names from the pulled KEGG entries
        names = {k: extract_entry_name(v) for k,v in records.items()}
        for h in self.hits:
            h.name = names[h.kegg_id]
        
        ## Extract genomic positions from the pulled KEGG entries
        positions = {k: extract_positional_information(v) for k,v in records.items()}
        # Extract a mapping table for RefSeq Nucleotide IDs from the Uniprot crossref mapping table
        all_uniprot_refseq = prepare_mapping_dict(all_cross_refs, 'RefSeq_NT')
        # Start sorting out the genomic information
        processed_hits = []
        for h in self.hits:
            position_info = positions[h.kegg_id]
            # If there are no numeric characters, there is no positional information, which makes this hit useless.
            # Therefore, it is not included in the processed hit set
            if not(any(char.isdigit() for char in position_info)):
                continue
            # If there is no colon, then there probably is no scaffold information,
            # but we have a second data source in our Uniprot crossref mapping table
            elif ':' not in position_info:
                refseq_id = all_uniprot_refseq[h.uniprot] # This might result in a double crossref, which will be fixed later on
                first_coor, second_coor = position_info.split('..')
            # Check whether it has the format of a complete genomic position field as we would like to have it from KEGG
            elif len(re.findall(r'.+:[complement]?\(?[0-9]+\.\.[0-9]+\)?', position_info)) > 1:
                refseq_id, coords = [position_info.split(':')] # wrapped in a list for consistency with the output of the case above
                first_coor, second_coor = coords.split('..')
            # Something else not encountered so far
            else:
                print(f'Invalid positional record found for {h.kegg_id}! Skipping...')
                continue
                        
            # Set the scaffold information      
            h.scaff = refseq_id
            
            # Need to do some further processing in case of a reverse strand gene
            if 'complement' in first_coor:
                h.strand = "-"
                h.start = int(first_coor[11:])
                h.end = int(second_coor[:-1])
            else:
                h.strand = '+'
                h.start = int(first_coor)
                h.end = int(second_coor)
                
            # Stack the processed hits
            processed_hits.append(h)
                
        # Update the hit set and split double crossrefs in the scaff attribute
        self.hits = processed_hits
        split_multiple_crossrefs(self, 'scaff')   
    
    def identify_clusters(self):
        pass
    
    def generate_output(self):
        pass
    
a = Hit('K0EJ92', 'query1')
b = Hit('K0EVU8', 'query2')
c = Cluster([a,b], 1)
s = Search({'query1': '/home/lucas/bin/cfoldseeker_old/AF3_models/fold_sco5088_model_0.cif',
            'query2': '/home/lucas/bin/cfoldseeker_old/AF3_models/fold_sco5089_model_0.cif'},
           {}, "uniprot_kegg_nucl.gz", hits = [a,b])
# s.crossref()
res = s.run_foldseek()
