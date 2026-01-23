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
    
    def __init__(self, uniprot, query, kegg_id = [], name = '', taxon = '', 
                 db = "", evalue = 1, prob = 1, score = 0, seqid = 0, qcov = 0,
                 scaff = '', start = 0, end = 0, strand = ''):
        
        self.uniprot: str = uniprot #Uniprot ID
        self.query: str = query #ID of the homologous query protein
        self.kegg_id: list = kegg_id #KEGG ID
        self.name: list = name #name of the hit in the KEGG entry
        self.taxon: str = taxon #name of the taxon in which this hit was found
        self.db: str = db #Structure database the hit was found in
        self.evalue: float = evalue #evalue of the FoldSeek hit
        self.prob: float = prob #FoldSeek hit probability score
        self.score: float = score #FoldSeek score
        self.seqid: float = seqid #Sequence identity with the query protein
        self.qcov: float = qcov #Query coverage with the query protein
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
        self.taxon: str = list(set([h.taxon for h in self.hits]))[0]
        
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
        
        self.hits = all_results
        
        return None
    
    
    def parse_foldseek_results(self, max_eval, min_prob, min_score, min_seqid, min_qcov):
        all_hits = []
        for query, task in self.hits.items():
            for results_by_db in task['results']:
                db = results_by_db['db']
                for hit_entry in results_by_db['alignments'][0]:
                    target = hit_entry['target']
                    uniprot = target.split('-')[1]
                    name = ' '.join(target.split(' ')[1:])
                    taxon = hit_entry['taxName']
                    evalue = float(hit_entry['eval'])
                    prob = float(hit_entry['prob'])
                    score = float(hit_entry['score'])
                    seqid = float(hit_entry['seqId'])
                    qcov = (int(hit_entry['qEndPos']) - int(hit_entry['qStartPos']))/int(hit_entry['qLen'])*100
                    if evalue <= max_eval and prob >= min_prob and score >= min_score and seqid >= min_seqid and qcov >= min_qcov:
                        hit = Hit(uniprot, query, name = name, taxon = taxon, db = db, evalue = evalue,
                                  prob = prob, score = score, seqid = seqid, qcov = qcov)
                        all_hits.append(hit)
            
        self.hits = all_hits
        
        return None
            
                
    def crossref(self):
        """
        Extracts the genomic positional information from a full pulled KEGG record.
        """
        def extract_positional_information(entry: str) -> str:
            try:
                return [line for line in entry.split('\n') if 'POSITION' in line][0][12:]
            except IndexError:
                return None
        
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
        Auxiliary function that standardises the hits. In case of multiple values for a given attribute, it splits it into 
        multiple hits that are identical except for that exact crossref attribute. The split hits are first introduced instead of the old hit before the full
        hit list is being flattened.
        In case of an absent or default attribute value, the hit is discarded.
        
        Mutates:
            self.hits: list: the hit list containing Hit objects that represent proteins that are structurally homologous to
                             one of the query proteins.
        """
        def sanitise_attr(self, attr):
            standardised_hits = []
            # Loop over all hits, but keep track of the index
            for h in self.hits:
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
                    standardised_hits.append(split_records)
                # If we don't need to split the record, just unlist it
                elif nb_copies == 1:
                    setattr(h, attr, getattr(h, attr)[0])
                    # Insert at the index of the old record
                    standardised_hits.append([h])
                # If it's empty, don't include it in the new hit set
                else:
                    continue
            
            # Flatten the full hit list again
            self.hits = list(chain(*standardised_hits))
            
            return None
        
        
        ## First, fill all KEGG IDs
        # Extract a mapping table for KEGG IDs from the Uniprot crossref mapping table
        all_uniprot_ids = [h.uniprot for h in self.hits]
        all_cross_refs = self.mapping_table.filter(self.mapping_table['Uniprot_ID'].is_in(all_uniprot_ids))
        all_uniprot_kegg = prepare_mapping_dict(all_cross_refs, 'KEGG')
        # Fill the KEGG IDs, if possible
        for h in self.hits:
            try:
                h.kegg_id = all_uniprot_kegg[h.uniprot]
            except KeyError:
                continue
        # Split records with multiple crossrefs and discard the ones without crossref
        sanitise_attr(self, 'kegg_id')
            
        ## Then, pull all KEGG records
        all_kegg_ids = [h.kegg_id for h in self.hits]
        pull = kgp.MultiProcessMultiplePull(n_workers = 2)
        _, records = pull.pull_dict(all_kegg_ids)
        
        ## Extract genomic positions from the pulled KEGG entries
        positions = {k: extract_positional_information(v) for k,v in records.items()}
        # Extract a mapping table for RefSeq Nucleotide IDs from the Uniprot crossref mapping table
        all_uniprot_refseq = prepare_mapping_dict(all_cross_refs, 'RefSeq_NT')
        # Start sorting out the genomic information
        processed_hits = []
        for h in self.hits:
            position_info = positions[h.kegg_id]
            # If there is not even position info, there is no point of keeping the hit
            if position_info == None:
                continue
            # If there are no numeric characters, there is no positional information
            if not(any(char.isdigit() for char in position_info)):
                continue
            # If there is no colon, then there probably is no scaffold information,
            # but we have a second data source in our Uniprot crossref mapping table
            elif ':' not in position_info:
                try:
                    refseq_id = all_uniprot_refseq[h.uniprot] # This might result in a double crossref, which will be fixed later on
                    first_coor, second_coor = position_info.split('..')
                # Bad luck, there's no RefSeq crossref in Uniprot either
                except KeyError:
                    continue
            # TODO: decide on whether we will support non-bacterial analyses. The cross-referencing is more cumbersome in this case
            # Check whether it has the format of a complete genomic position field as we would like to have it from KEGG
            elif len(re.findall(r'.+:.*[0-9]+\.\.[0-9]+.*', position_info)) > 0:
                refseq_id, coords = position_info.split(':') # wrapped in a list for consistency with the output of the case above
                # Just take the minimal and maximum coordinate. For introns, this might give issues with the gene/exon length.
                coord_groups = re.findall(r'\d+\.\.\d+', coords)
                coord_groups = list(chain(*[i.split('..') for i in coord_groups]))
                first_coor, second_coor = min(coord_groups), max(coord_groups)
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
                
        # Update the hit set and sanitise the scaff attribute
        self.hits = processed_hits
        sanitise_attr(self, 'scaff')   
    
    def identify_clusters(self):
        pass
    
    def generate_output(self):
        pass
    
a = Hit('K0EJ92', 'query1')
b = Hit('K0EVU8', 'query2')
c = Cluster([a,b], 1)
s = Search({'query1': '/home/lucas/bin/cfoldseeker_old/AF3_models/fold_sco5088_model_0.cif',
            'query2': '/home/lucas/bin/cfoldseeker_old/AF3_models/fold_sco5089_model_0.cif'},
           {}, mapping_table_path = "uniprot_kegg_nucl.gz")
s.run_foldseek()
s.parse_foldseek_results(1, 0, 0, 25, 70)
s.crossref()