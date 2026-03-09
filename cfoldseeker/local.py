#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import logging
import threading
import subprocess
import shutil
import sys
import polars as pl
from pathlib import Path

from cfoldseeker.classes import Search, Hit

LOG = logging.getLogger(__name__)


def _stream_reader(pipe, write_func):
    try:
        with pipe:
            for chunk in iter(lambda: pipe.readline(), b''):
                if not chunk:
                    break
                try:
                    text = chunk.decode('utf-8', 'replace')
                except Exception:
                    text = chunk.decode('latin-1', 'replace')
                write_func(text)
    except Exception:
        LOG.exception("stream reader error")
        

class LocalSearch(Search):
    
    def __init__(self, query, db_path, coord_db_path, params = {}, hits = [], clusters = [],
                 output_folder = Path('.'), temp_folder = Path('.')):
        
        super().__init__(query, params, hits, clusters, output_folder, temp_folder)
        
        self.db_path: Path = db_path
        
        LOG.debug(f"Scanning CDS DB from {coord_db_path}")
        self.coord_db: pl.LazyFrame = pl.scan_csv(coord_db_path, has_header = True, separator = "\t")
        
        return None
    
    def __repr__(self) -> str:
        return f"Local Search of {','.join(list(self.query.keys()))} with {len(self.clusters)} clusters identified"
    
    
    def run_foldseek(self) -> None:
        """
        Runs FoldSeek on the local DB.
        """

        foldseek_executable = shutil.which('foldseek')
        foldseek_verbosity = str(min(self.params['verbosity'], 3))
        
        LOG.debug(f"FoldSeek executable: {foldseek_executable}")
        LOG.debug(f"FoldSeek target DB: {str(self.db_path)}")
        LOG.debug(f"FoldSeek workdir: {str(self.TEMP_DIR / 'foldseek_tmp')}")
        LOG.debug(f'Applying minimum sequence identity threshold >= {str(self.params["min_seqid"])}')
        LOG.debug(f'Applying maximum evalue threshold <= {str(self.params["max_eval"])}')
        
        cmd = [foldseek_executable, 'easy-search', 
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
            LOG.critical(f"FoldSeek exited with code {returncode}")
            sys.exit()
        else:
            LOG.info('FoldSeek finished successfully.')
    
        return None        
    
    
    def parse_foldseek_results(self) -> None:
        """
        Parses the FoldSeek result table.
        """
        
        ## Load the thresholds from params
        min_score = self.params['min_score']
        min_qcov = self.params['min_qcov']
        min_tcov = self.params['min_tcov']
        
        ## Parse results table
        LOG.debug(f"Scanning FoldSeek result table at {str(self.TEMP_DIR / 'foldseek_result.txt')}")
        results = pl.scan_csv(self.TEMP_DIR / 'foldseek_result.txt', has_header = True, separator = "\t")
        results = results.unique() # Discard duplicate hits
        
        ## Convert tcov and qcov to percentages
        results = results.with_columns([pl.col("qcov") * 100, pl.col('tcov') * 100])
        
        ## Filter hits
        LOG.debug('Applying the following hit filtering thresholds:')
        LOG.debug(f'bitscore >= {min_score}')
        LOG.debug(f'query coverage >= {min_qcov}')
        LOG.debug(f'target coverage >= {min_tcov}')
        
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
        Complete local search workflow run
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
        
    