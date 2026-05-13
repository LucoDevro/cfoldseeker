#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
import polars as pl
from pathlib import Path

from cfoldseeker.local import LocalSearch


LOG = logging.getLogger(__name__)


class LocalClusteredSearch(LocalSearch):
    
    def __init__(self, query, db_path, coord_db_path, params = {}, hits = [], clusters = [],
                 output_flags = {}, output_folder = Path('.'), temp_folder = Path('.'),
                 seq_clust_tsv = Path('.')):
        """
        Initialise a LocalClusteredSearch instance with database paths and parameters.
        
        Calls the parent LocalSearch class constructor and loads the CDS coordinates
        database from a TSV file.
        
        Args:
            query: A dictionary of query protein structures and their filepaths (inherited from Search).
            db_path: Path to the local FoldSeek protein structure target database.
            coord_db_path: Path to the CDS coordinates database file (tab-separated,
                no header).
            params: Dictionary of search parameters (e.g., min_seqid, max_eval,
                min_score). Defaults to empty dict.
            hits: List of Hit objects from previous searches. Defaults to empty list.
            clusters: List of gene clusters from previous analysis. Defaults to
                empty list.
            output_flags (dict, optional): Output parameters dictionary. Defaults to {}.
            output_folder: Path to folder for output files. Defaults to current
                directory.
            temp_folder: Path to folder for temporary files. Defaults to current
                directory.
            seq_clust_tsv: Path to the MMseqs2 clustering table. Defaults to current
                directory.
                
        Note:
            The coord_db and the seq_clust are loaded as LazyFrames for memory efficiency when dealing
            with large databases.
        """
           
        super().__init__(query, db_path, coord_db_path, params, hits, clusters, output_flags,
                         output_folder, temp_folder)
           
        self.seq_clust = pl.scan_csv(seq_clust_tsv, has_header = False, separator = "\t",
                                     new_columns = ['representative', 'protein']).unique()
           
        return None
   
    
    def __repr__(self) -> str:
        return f"Local Clustered Search of {','.join(list(self.query.keys()))} with {len(self.clusters)} clusters identified"
    
    
    def expand_sequence_clusters(self, results: pl.DataFrame) -> pl.DataFrame:
        """
        Expand a given FoldSeek result table with all members of the original hits' sequence clusters.
        
        Includes all members of the sequence clusters of which the original FoldSeek hits are the 
        representatives of by joining the FoldSeek result table with the MMseqs2 clustering table.
        
        Drops duplicate protein/query pairs.
        
        Hit metadata are taken over from the representative protein.
        
        Args:
            results (polars.DataFrame): Original FoldSeek result table with only the representative proteins
            
        Returns:
            expanded_results (polars.DataFrame): Result table with all added non-representative proteins 
            for each representative
        """
        LOG.info('Expanding the hit set with sequence cluster members')
        
        # Identify which genes remained after filtering
        filtered_genes = results.select('target').unique()
        
        # Prefilter clustering table for present genes
        LOG.info('Prefiltering sequence clustering table')
        seq_clust_filt = self.seq_clust.filter(pl.col('representative').is_in(filtered_genes['target']))
        seq_clust_filt = seq_clust_filt.collect() # Materialise for join efficiency
        
        # Include all sequence cluster members by joining with the MMseqs2 clustering table
        LOG.info('Adding sequence cluster members')
        expanded_results = seq_clust_filt.join(results, left_on = "representative", right_on = 'target')
        
        # Keep only one hit in case a protein/query pair reoccurs
        # (e.g. different proteins in different assemblies with identical labels at NCBI)
        expanded_results = expanded_results.drop('representative').unique(subset = ['protein', 'query'])
        
        # Rename the protein column back to target as in the original results dataframe
        expanded_results = expanded_results.rename({'protein': 'target'})
        
        return expanded_results
    
    
    def identify_hits(self) -> None:
        """
        Identify hits passing the hit thresholds from the FoldSeek results.
        
        Parses the FoldSeek results table and applies hit-level filtering, then
        fetches genomic context information for each hit from the context DB,
        and collects freshly instantiated Hit objects to host all metadata.
        
        Returns:
            None
            
        Mutates:
            self.hits: Instantiates the list of identified Hit objects.
        """
        # Parses FoldSeek results and applies hit-level filtering
        parsed_results = self.parse_foldseek_results()
        
        # Adds all sequence cluster members for each hit
        expanded_results = self.expand_sequence_clusters(parsed_results)
        
        # Adds genomic context and instantiates Hit objects
        self.collect_hits(expanded_results)
        
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
        
        LOG.info("STARTING PART 2: Identifying hits in FoldSeek results")
        self.identify_hits()
        LOG.info('FINISHED PART 2')
        
        LOG.info("STARTING PART 3: Identifying gene clusters")
        self.identify_clusters()
        LOG.info('FINISHED PART 3')
        
        return None
    
    