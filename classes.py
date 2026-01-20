#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import kegg_pull.pull as kgp
import polars as pl
from pathlib import Path
from itertools import chain
from copy import deepcopy

class Hit:
    
    def __init__(self, uniprot, query, kegg_id = '', name = '', scaff = '', start = 0, end = 0, strand = ''):
        
        self.uniprot: str = uniprot
        self.query: str = query
        self.kegg_id: list = kegg_id
        self.name: list = name
        self.scaff: str = scaff
        self.start: int = int(start)
        self.end: int = int(end)
        self.strand: str = strand
    
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
        
        self.query: dict = query
        self.params: dict = params
        self.hits: list = hits
        self.clusters: list = clusters
        self.mapping_table: dict = pl.read_csv("uniprot_kegg_nucl.gz", has_header = False, separator = "\t",
                                               new_columns = ['Uniprot_ID', 'Target_DB', 'ID'])
        
        
        return None
    
    def __repr__(self):
        pass
    
    def foldseek(self):
        pass
    
    def crossref(self):
        
        def extract_entry_name(entry):
            return [line for line in entry.split('\n') if 'NAME' in line][0][12:]
        
        def extract_positional_information(entry):
            return [line for line in entry.split('\n') if 'POSITION' in line][0][12:]
        
        def prepare_mapping_dict(crossref_df, target_db):
            res = crossref_df.filter(crossref_df['Target_DB'] == target_db).drop('Target_DB')
            res = res.group_by('Uniprot_ID').agg('ID')
            res = dict(zip(res['Uniprot_ID'], res['ID']))
            return res
        
        # First, fill all KEGG IDs
        all_uniprot_ids = [h.uniprot for h in self.hits]
        all_cross_refs = self.mapping_table.filter(self.mapping_table['Uniprot_ID'].is_in(all_uniprot_ids))
        all_uniprot_kegg = prepare_mapping_dict(all_cross_refs, 'KEGG')
        for h in self.hits:
            h.kegg_id = all_uniprot_kegg[h.uniprot]
        # Split records with multiple crossrefs
        for idx, h in enumerate(self.hits):
            nb_copies = len(h.kegg_id)
            if nb_copies > 1:
                split_records = []
                copy_counter = 0
                while copy_counter < nb_copies:
                    new_record = deepcopy(h)
                    new_record.kegg_id = h.kegg_id[copy_counter]
                    split_records.append(new_record)
                    copy_counter += 1
                self.hits[idx] = split_records
            else:
                h.kegg_id = h.kegg_id[0]
                self.hits[idx] = [h]
        self.hits = list(chain(*self.hits))
            
        # Then, pull all KEGG records
        all_kegg_ids = [h.kegg_id for h in self.hits]
        pull = kgp.MultiProcessMultiplePull(n_workers = 2)
        _, records = pull.pull_dict(all_kegg_ids)
        
        # Extract the record names
        names = {k: extract_entry_name(v) for k,v in records.items()}
        for h in self.hits:
            h.name = names[h.kegg_id]
        
        # Extract genomic positions from the KEGG records
        positions = {k: extract_positional_information(v) for k,v in records.items()}
        all_uniprot_refseq = prepare_mapping_dict(all_cross_refs, 'RefSeq_NT')
        for h in self.hits:
            position_info = positions[h.kegg_id]
            # If there are no numeric characters, there is no positional information
            if not(any(char.isdigit() for char in position_info)):
                h.scaff = h.start = h.stop = None
            # If there is no colon, then there probably is no scaffold information
            elif ':' not in position_info:
                refseq_id = all_uniprot_refseq[h.uniprot][0]
                first_coor, second_coor = position_info.split('..')
            # Check whether it has the format of a complete genomic position field
            elif len(re.findall(r'.+:[complement]?\(?[0-9]+\.\.[0-9]+\)?', position_info)) > 1:
                refseq_id, coords = position_info.split(':')
                first_coor, second_coor = coords.split('..')
            else:
                print('Invalid positional record found!')
                continue
                        
            h.scaff = refseq_id
            if 'complement' in first_coor:
                h.strand = "-"
                h.start = int(first_coor[11:])
                h.end = int(second_coor[:-1])
            else:
                h.strand = '+'
                h.start = int(first_coor)
                h.end = int(second_coor)
            
    
    def identify_clusters(self):
        pass
    
    def generate_output(self):
        pass
    
a = Hit('K0EJ92', 'query1')
b = Hit('K0EVU8', 'query2')
c = Cluster([a,b], 1)
s = Search({}, {}, "uniprot_kegg_nucl.gz", hits = [a,b])
s.crossref()
