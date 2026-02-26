#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import subprocess
import polars as pl

from classes import Search, Hit

class LocalSearch(Search):
    
    def __init__(self, query, db_path, coord_db_path, params = {}, hits = [], clusters = []):
        
        super().__init__(query, params, hits, clusters)
        
        self.db_path = db_path
        
        self.coord_db = pl.scan_csv(coord_db_path, has_header = False, separator = "\t", 
                                    new_columns = ['gene_tag', 'name', 'contig', 'coords',
                                                   'strand', 'taxon_id', 'taxon_name'])
        
        return None
    
    def __repr__(self):
        return f"Local Search of {','.join(list(self.query.keys()))} with {len(self.clusters)} clusters identified"
    
    
    def run_foldseek(self):
        """
        Runs FoldSeek on the local DB.
        """
        
        cmd = ['foldseek', 'easy-search', 
               ' '.join(self.query.values()), 
               self.db_path,
               self.TEMP_DIR / 'foldseek_result.txt',
               self.TEMP_DIR / 'foldseek_tmp',
               "--format-mode", '4',
               "--input-format", '2',
               '--format-output', 'query,target,qstart,qend,tstart,tend,pident,qcov,tcov,prob,evalue,bits',
               '-v', '1',
               "--min-seq-id", str(self.params['min_seqid']),
               '-e', str(self.params['max_eval']),
               '--threads', self.params['cores']
               ]
        subprocess.run(cmd, check = True)
        
        return None
    
    
    def parse_foldseek_results(self):
        """
        Parses the FoldSeek result table.
        """
        
        ## Load the thresholds from params
        min_prob = self.params['min_prob']
        min_score = self.params['min_score']
        min_qcov = self.params['min_qcov']
        min_tcov = self.params['min_tcov']
        
        ## Parse results table
        results = pl.scan_csv(self.TEMP_DIR / 'foldseek_result.txt', has_header = True, separator = "\t")
        
        ## Convert tcov and qcov to percentages
        results = results.with_columns([pl.col("qcov") * 100, pl.col('tcov') * 100])
        
        ## Filter hits
        results = results.filter(pl.col('prob') >= min_prob)
        results = results.filter(pl.col('bits') >= min_score)
        results = results.filter(pl.col('qcov') >= min_qcov)
        results = results.filter(pl.col('tcov') >= min_tcov)
        
        ## Join with coordinates DB
        results = results.join(self.coord_db, left_on = 'target', right_on = 'gene_tag', maintain_order = "left_right")
        
        ## Discard the alignment coordinates columns
        results = results.drop("tstart", "tend", "qstart", "qend")
        
        ## Materialise the LazyFrame and turn into a row iterator
        results = results.collect().iter_rows(named = True)

        ## Generate the Hit objects
        all_hits = []
        for result in results:
            # Parse the genomic coordinates on the fly
            coord_groups = re.findall(r'\d+\.\.\d+', result['coords'])
            coord_groups = [i.split('..') for i in coord_groups]
            coord_groups = [[int(j) for j in i] for i in coord_groups]
            
            hit = Hit(db_id = result['target'],
                      crossref_id = result['target'],
                      query = result['query'],
                      name = result['name'],
                      taxon_name = result['taxon_name'],
                      taxon_id = result['taxon_id'],
                      db = ["local"],
                      crossref_method = 'local',
                      evalue = result["evalue"],
                      prob = result['prob'],
                      score = result["bits"],
                      seqid = result["pident"],
                      qcov = result["qcov"],
                      tcov = result["tcov"],
                      scaff = result['contig'],
                      strand = result['strand'],
                      coords = coord_groups)
            all_hits.append(hit)
        
        self.hits = all_hits
        
        return None
    
    
    def run(self):
        """
        Complete local search workflow run
        """
        
        self.run_foldseek()
        self.parse_foldseek_results()
        self.identify_clusters()
        
        return None
        
    