#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import re
import polars as pl
from pathlib import Path

from cfoldseeker.local import LocalSearch
from cfoldseeker.classes import Hit


LOG = logging.getLogger(__name__)


class LocalClusteredSearch(LocalSearch):
    
    def __init__(self, query, db_path, coord_db_path, params = {}, hits = [], clusters = [],
                 output_folder = Path('.'), temp_folder = Path('.'), seq_clust_tsv = Path('.')):
           
        super().__init__(query, db_path, coord_db_path, params, hits, clusters,
                         output_folder, temp_folder)
           
        self.seq_clust = pl.scan_csv(seq_clust_tsv, has_header = False, separator = "\t",
                                     new_columns = ['representative', 'protein']).unique()
           
        return None
   
    
    def __repr__(self) -> str:
        return f"Local Clustered Search of {','.join(list(self.query.keys()))} with {len(self.clusters)} clusters identified"
    
    
    def expand_sequence_clusters(self, foldseek_result: pl.LazyFrame) -> pl.LazyFrame:
        """
        Expand a given FoldSeek result table with all members of the original hits' sequence clusters.
        
        Includes all members of the sequence clusters of which the original FoldSeek hits are the 
        representatives of by joining the FoldSeek result table with the MMseqs2 clustering table.
        
        Drops duplicate protein/query pairs.
        
        Hit metadata are taken over from the representative protein.
        
        Args:
            foldseek_result (polars.LazyFrame): Original FoldSeek result table with only the representative proteins
            
        Returns:
            expanded_results (polars.LazyFrame): Result table with all added non-representative proteins 
            for each representative
        """
        # Include all sequence cluster members by joining with the MMseqs2 clustering table
        expanded_results = self.seq_clust.join(foldseek_result, left_on = "representative", right_on = 'target')
        
        # Keep only one hit in case a protein/query pair reoccurs
        # (e.g. different proteins in different assemblies with identical labels at NCBI)
        expanded_results = expanded_results.drop('representative').unique(subset = ['protein', 'query'])
        
        return expanded_results
    
    
    def parse_foldseek_results(self):
        """
        Parse the FoldSeek result table, expand it to include all members of the sequence clusters,
        and generate Hit objects with filled genomic coordinates.
        
        Reads the FoldSeek result table, expands it with all members of each sequence cluster of
        which the original FoldSeek hits are the representatives (by joining with the clustering table),
        applies filtering thresholds (bit score, query coverage, target coverage), removes duplicate hits,
        and joins results with the CDS coordinates database. Parses genomic coordinates from the
        coordinate string and creates Hit objects for each match.
        
        The following filtering thresholds are applied:
        1. Sequence identity >= min_seqid
        2. E-value <= max_eval
        3. Bit score >= min_score
        4. Query coverage >= min_qcov (converted to percentage)
        5. Target coverage >= min_tcov (converted to percentage)
        
        Returns:
            None
            
        Note:
            Stores generated Hit objects in self.hits as a list. Genomic coordinates
            are parsed from a comma-separated string of joined range pairs (e.g., "10..50",
            "join(150..200,250..300)") into nested lists of integers.
        """
        ## Load the thresholds from params
        min_score = self.params['min_score']
        min_qcov = self.params['min_qcov']
        min_tcov = self.params['min_tcov']
        min_seqid = self.params['min_seqid']
        max_eval = self.params['max_eval']
        
        ## Parse results table
        LOG.debug(f"Scanning FoldSeek result table at {str(self.TEMP_DIR / 'foldseek_result.txt')}")
        results = pl.scan_csv(self.TEMP_DIR / 'foldseek_result.txt', has_header = True, separator = "\t")
        results = results.unique() # Discard duplicate hits
        
        ## Expanding sequence clusters
        LOG.info('Expanding the sequence clusters')
        results = self.expand_sequence_clusters(results)
        
        ## Convert tcov and qcov to percentages
        results = results.with_columns([pl.col("qcov") * 100, pl.col('tcov') * 100])
        
        ## Filter hits
        LOG.debug('Applying the following hit filtering thresholds:')
        LOG.debug(f'sequence identity >= {min_seqid}')
        LOG.debug(f'evalue <= {max_eval}')
        LOG.debug(f'bitscore >= {min_score}')
        LOG.debug(f'query coverage >= {min_qcov}')
        LOG.debug(f'target coverage >= {min_tcov}')
        
        results = results.filter(pl.col('pident') >= min_seqid)
        results = results.filter(pl.col('evalue') <= max_eval)
        results = results.filter(pl.col('bits') >= min_score)
        results = results.filter(pl.col('qcov') >= min_qcov)
        results = results.filter(pl.col('tcov') >= min_tcov)
        
        ## Join with coordinates DB
        LOG.info('Fetching CDS coordinates from local CDS DB')
        results = results.join(self.coord_db, left_on = 'protein', right_on = 'gene_tag')
        
        ## Discard the alignment coordinates columns
        results = results.drop("tstart", "tend", "qstart", "qend")
        
        ## Materialise the LazyFrame
        results = results.collect()
        LOG.info(f"Found {results.height} gene hits.")

        ## Generate the Hit objects
        LOG.debug('Generating the Hit objects')
        results_it = results.iter_rows(named = True)
        all_hits = []
        for result in results_it:
            # Parse the genomic coordinates on the fly
            coord_groups = re.findall(r'\d+\.\.\d+', result['coords'])
            coord_groups = [i.split('..') for i in coord_groups]
            coord_groups = [[int(j) for j in i] for i in coord_groups]
            
            hit = Hit(db_id = result['protein'],
                      crossref_id = result['protein'],
                      query = result['query'],
                      name = result['name'],
                      taxon_name = result['taxon_name'],
                      taxon_id = result['taxon_id'],
                      db = "local_clustered",
                      crossref_method = 'local_clustered',
                      evalue = result["evalue"],
                      score = result["bits"],
                      seqid = result["pident"],
                      qcov = result["qcov"],
                      tcov = result["tcov"],
                      scaff = result['contig'],
                      strand = result['strand'],
                      coords = coord_groups)
            all_hits.append(hit)
        
        self.hits = all_hits
        
        LOG.info(f'{len(all_hits)} hits have been processed.')
        
        return None
    
    
    def run(self) -> None:
        """
        Execute the complete local search workflow.
        
        Orchestrates all processing steps in sequence: running FoldSeek locally against
        the local database, parsing the FoldSeek results and creating Hit objects for
        hits passing the hit criteria, and identifying gene clusters from the hits.
        
        Returns:
            None
        """
        
        LOG.info('STARTING PART 1: Executing FoldSeek search')
        self.run_foldseek()
        LOG.info("FINISHED PART 1")
        
        LOG.info("STARTING PART 2: Parsing FoldSeek results")
        self.parse_foldseek_results()
        LOG.info('FINISHED PART 2')
        
        LOG.info("STARTING PART 3: Identifying gene clusters")
        self.identify_clusters()
        LOG.info('FINISHED PART 3')
        
        return None
    
    