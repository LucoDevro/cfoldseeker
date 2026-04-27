#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import logging
import threading
import subprocess
import shutil
import polars as pl
from pathlib import Path

from cfoldseeker.classes import Search, Hit


LOG = logging.getLogger(__name__)


def _stream_reader(pipe, write_func):
    """
    Reads from a pipe stream and writes decoded text to a callback function.
    
    Continuously reads binary chunks from a pipe, attempts to decode them as UTF-8),
    and passes the decoded text to a write function. Handles encoding errors
    gracefully and logs exceptions.
    
    Args:
        pipe: A binary stream pipe (typically from subprocess.PIPE) to read from.
        write_func: A callable that accepts a string argument. Called with each
            decoded line of text from the pipe.
            
    Note:
        This function is designed to be run in a separate thread to capture
        subprocess output without blocking. Exceptions during reading are logged
        but do not raise; the function exits silently on pipe closure or errors.
    """
    try:
        with pipe:
            for chunk in iter(lambda: pipe.readline(), b''):
                if not chunk:
                    break
                text = chunk.decode('utf-8', 'replace')
                write_func(text)
    except Exception:
        LOG.exception("stream reader error")
        

class LocalSearch(Search):
    """
    Subclass executing the workflow for gene cluster identification from local protein
    searches using FoldSeek.
    
    Extends the Search base class to perform searches against local FoldSeek databases.
    Handles FoldSeek execution, result parsing, Hit object generation, and gene cluster
    identification.
    Uses a TSV of CDS coordinates made beforehand with cfoldseeker-cds.
    
    Attributes:
        db_path (Path): Path to the FoldSeek protein structure target database.
        coord_db (polars.LazyFrame): DataFrame containing CDS coordinates
            with columns: gene_tag, name, contig, strand, coords, taxon_id,
            taxon_name.
    """
    def __init__(self, query, db_path, coord_db_path, params = {}, hits = [], clusters = [],
                 output_folder = Path('.'), temp_folder = Path('.')):
        """
        Initialise a LocalSearch instance with database paths and parameters.
        
        Calls the parent Search class constructor and loads the CDS coordinates
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
            output_folder: Path to folder for output files. Defaults to current
                directory.
            temp_folder: Path to folder for temporary files. Defaults to current
                directory.
                
        Note:
            The coord_db is loaded as a LazyFrame for memory efficiency when dealing
            with large databases.
        """
        super().__init__(query, params, hits, clusters, output_folder, temp_folder)
        
        self.db_path: Path = db_path
        
        LOG.debug(f"Scanning CDS DB from {coord_db_path}")
        self.coord_db: pl.LazyFrame = pl.scan_csv(coord_db_path, has_header = False, separator = "\t",
                                                  new_columns = ['gene_tag', 'name', 'contig', 'strand', 
                                                                 'coords', 'taxon_id', 'taxon_name'])
        
        return None
    
    
    def __repr__(self) -> str:
        return f"Local Search of {','.join(list(self.query.keys()))} with {len(self.clusters)} clusters identified"
    
    
    def run_foldseek(self) -> None:
        """
        Execute a FoldSeek search against the local protein structure database.
        
        Constructs and runs a FoldSeek 'easy-search' command with all query structures
        (in CIF format) against the local database. Applies filters for sequence identity
        and E-value thresholds. Captures stdout and stderr in real-time via
        separate threads and logs them appropriately.
        
        Exhaustive search (no database prefiltering) has been enabled to retrieve
        all hits in the target database.
        
        Returns:
            None
            
        Raises:
            RuntimeError: If FoldSeek returns a non-zero exit code.
            
        Note:
            FoldSeek output is written to a temporary TSV file in TEMP_DIR.
        """

        foldseek_executable = Path(shutil.which('foldseek'))
        foldseek_verbosity = str(min(self.params['verbosity'], 3))
        
        LOG.debug(f"FoldSeek executable: {foldseek_executable}")
        LOG.debug(f"FoldSeek target DB: {str(self.db_path)}")
        LOG.debug(f"FoldSeek workdir: {str(self.TEMP_DIR / 'foldseek_tmp')}")
        LOG.debug(f'Applying minimum sequence identity threshold >= {str(self.params["min_seqid"])}')
        LOG.debug(f'Applying maximum evalue threshold <= {str(self.params["max_eval"])}')
        
        cmd = [str(foldseek_executable), 'easy-search', 
               *[str(q) for q in self.query.values()], 
               str(self.db_path),
               str(self.TEMP_DIR / 'foldseek_result.txt'),
               str(self.TEMP_DIR / 'foldseek_tmp'),
               "--format-mode", '4',
               "--input-format", '2',
               '--format-output', 'query,target,qstart,qend,tstart,tend,pident,qcov,tcov,evalue,bits',
               '-v', foldseek_verbosity,
               "--min-seq-id", str(self.params['min_seqid']),
               '-e', str(self.params['max_eval']),
               '--threads', str(self.params['cores']),
               '--exhaustive-search', '1'
               ]
        
        # Launching search process
        LOG.debug(f'Running command: {" ".join(cmd)}')
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        # Capture stdout and stderr in realtime and wrap it in the logs
        def foldseek_stdout_log(s): return LOG.debug(s.rstrip())
        def foldseek_stderr_log(s): return LOG.warning(s.rstrip())
        
        t_out = threading.Thread(target=_stream_reader, args=(proc.stdout, foldseek_stdout_log))
        t_err = threading.Thread(target=_stream_reader, args=(proc.stderr, foldseek_stderr_log))
        t_out.daemon = True
        t_err.daemon = True
        t_out.start()
        t_err.start()
    
        returncode = proc.wait()
        t_out.join()
        t_err.join()
    
        # Wrap up
        if returncode != 0:
            msg = f"{foldseek_executable.name} failed with return code {returncode}."
            LOG.critical(msg)
            raise RuntimeError(msg)
        else:
            LOG.info(f'{foldseek_executable.name} finished successfully.')
    
        return None        
    
    
    def parse_foldseek_results(self) -> None:
        """
        Parse the FoldSeek result table and generate Hit objects with filled genomic coordinates.
        
        Reads the FoldSeek result table, applies filtering thresholds (bit score,
        query coverage, target coverage), removes duplicate hits, and joins results
        with the CDS coordinates database. Parses genomic coordinates from the
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
            are parsed from a comma-separated string of joined range pairs (e.g., "10..50,
            join(150..200,250..300)") into nested lists of integers.
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
        results = results.join(self.coord_db, left_on = 'target', right_on = 'gene_tag', maintain_order = "left_right")
        
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
            
            hit = Hit(db_id = result['target'],
                      crossref_id = result['target'],
                      query = result['query'],
                      name = result['name'],
                      taxon_name = result['taxon_name'],
                      taxon_id = result['taxon_id'],
                      db = "local",
                      crossref_method = 'local',
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
        
    