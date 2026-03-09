#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import shutil
import logging
import tempfile
from pathlib import Path
from cblaster.plot import plot_session
from cblaster.plot_clusters import plot_clusters
from importlib.metadata import version

from cfoldseeker.remote import RemoteSearch
from cfoldseeker.local import LocalSearch

__version__ = version("cfoldseeker")

LOG = logging.getLogger()


def getArguments() -> argparse.Namespace:
    """
    This function gets the CLI arguments, without any parsing.
    """
    
    parser = argparse.ArgumentParser(
        prog = 'cfoldseeker',
                epilog = 
                """
                Lucas De Vrieze, Miguel Biltjes
                (c) 2026 Masschelein lab, VIB
                """,
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = 
                """
                cfoldseeker: find clusters of colocalised genes using structural similarity with FoldSeek
                """,
                add_help = False
                )
    
    args_general = parser.add_argument_group('General')
    args_general.add_argument('-m', '--mode', dest = 'mode', default = 'remote', type = str, choices = ['local', 'remote'], help = "Search mode (default: remote)")
    args_general.add_argument('-c', '--cores', dest = 'cores', default = 1, type = int, help = "Number of cores available to use (default: 1)")          
    args_general.add_argument('-f', '--force', dest = 'force', default = False, action = 'store_true', help = "Force overwriting output (default: false).")
    args_general.add_argument('-vv', '--verbosity', dest = 'verbosity', default = 3, type = int, choices = [0,1,2,3,4], help = "Console verbosity level (default: 3 (info))")
    args_general.add_argument('-v', '--version', action = "version", version = "%(prog)s " + __version__)
    args_general.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit")
    
    args_io = parser.add_argument_group('Inputs and outputs')
    args_io.add_argument('-q', '--query', dest = 'query_folder', required = True, type = Path, help = "Path of the folder containing the query proteins.")
    args_io.add_argument('-o', '--output', dest = 'output', type = Path, default = Path('.'), help = "Output directory (default: current location)")
    args_io.add_argument('-t', '--temp', dest = "temp", type = Path, default = tempfile.gettempdir(), help = "Path to store temporary files (default: your OS's default temporary directory).")
    args_io.add_argument('--no-tables', dest = 'output_tables', default = True, action = 'store_false', help = "Don't write overview tables (default: false).")
    args_io.add_argument('--session', dest = 'output_session', default = False, action = 'store_true', help = "Write cblaster session file (default: false).")
    args_io.add_argument('--summary', dest = 'output_summary', default = False, action = 'store_true', help = "Write cblaster summary file (default: false).")
    args_io.add_argument('--binary', dest = 'output_binary', default = False, action = 'store_true', help = "Write cblaster binary file (tab-separated) (default: false).")
    args_io.add_argument('--plot', dest = 'output_plot', default = False, action = 'store_true', help = "Write cblaster clusterplot file (default: false).")
    args_io.add_argument('--clinker', dest = 'output_clinker', default = False, action = 'store_true', help = "Write clinker plot file (default: false).")
    args_io.add_argument('--foldseek', dest = 'output_foldseek', default = False, action = 'store_true', 
                         help = "Save FoldSeek output. In remote mode, this is the raw json from the FoldSeek webserver. In local mode, this is a BLAST-like tabular text file. (default: false).")

    args_search = parser.add_argument_group('General search options')
    args_search.add_argument('--max-eval', dest = "max_eval", type = float, default = 1e-9, help = "Maximum e-value to include a FoldSeek hit (default: 1e-9).")
    args_search.add_argument('--min-score', dest = "min_score", type = float, default = 250, help = "Minimum FoldSeek bitscore to include a hit (default: 250).")
    args_search.add_argument('--min-seqid', dest = "min_seqid", type = float, default = 0, help = "Minimum sequence identity to include a hit (in percentages) (default: 0).")
    args_search.add_argument('--min-qcov', dest = "min_qcov", type = float, default = 50, help = "Minimum query coverage to include a hit (in percentages) (default: 50).")
    args_search.add_argument('--min-tcov', dest = "min_tcov", type = float, default = 50, help = "Minimum target coverage to include a hit (in percentages) (default: 50).")
    args_search.add_argument('--max-gap', dest = "max_gap", type = int, default = 5000, help = "Maximum intergenic gap within a cluster (in bp) (default: 5000).")
    args_search.add_argument('--max-length', dest = "max_length", type = int, default = 1e5, help = "Maximum genomic length of a cluster (in bp) (default: 1e5).")
    args_search.add_argument('--min-hits', dest = "min_hits", type = int, default = 2, help = "Minimum number of members in a cluster (default: 2).")
    args_search.add_argument('--min-cov-qrs', dest = "min_cov_qrs", type = int, default = 2, help = "Minimum different queries covered by a cluster (default: 2).")
    args_search.add_argument('--require', dest = "require", type = str, default = '', nargs = '*', help = "Queries that have to present in a cluster (use filenames without extensions).")
    
    args_remote = parser.add_argument_group("Remote-specific search options")
    args_remote.add_argument('-db', '--database', dest = 'db', type = str, default = ['afdb50'], nargs = '+', choices = ['afdb-proteome', 'afdb-swissprot', 'afdb50'], 
                             help = "Remote target database (default: afdb50) (choices: afdb-proteome, afdb-swissprot, afdb50)")
    args_remote.add_argument('-tf', '--taxon-filter', dest = 'taxfilters', type = str, default = '', nargs = '*',
                             help = "Taxon ID(s) to filter the FoldSeek results table.")
    args_remote.add_argument('-uma', '--uniprot-mapping', dest = 'mapping_table_path', type = Path, default = Path('uniprot_kegg_genpept.gz'),
                             help = "Path to the UniProt AFDB ID mapping table (default: uniprot_kegg_genpept.gz)")
    args_remote.add_argument('-w', '--max-workers', dest = "max_workers", type = int, default = 1, help = "Maximum number of workers to query the remote servers (FoldSeek, KEGG, ENA) (default: 1)")
    
    args_local = parser.add_argument_group('Local-specific search options')
    args_local.add_argument('-ldb', '--local-database', dest = 'local_db_path', type = Path, default = Path('local_db/local_db'), help = "Path to your local FoldSeek DB (format: <path-to-containing-folder>/<DB-prefix>) (default: local_db/local_db).")
    args_local.add_argument('-cm', '--cds-mapping', dest = 'cds_mapping_table_path', type = Path, default = Path('local_cds_db.gz'), help = "Path of the CDS coordinates DB (default: local_cds_db.gz).")
    
    args = parser.parse_args()
    
    return args


def parseArguments(args) -> dict:
    """
    This function parses the arguments and returns a parsed arguments dictionary that is used to call the workflows.
    """
    
    assert args.mode in ['local', 'remote'], 'Invalid search mode. Possible choices: "local" and "remote".'
    assert set(args.db) <= {'local', 'afdb-proteome', 'afdb-swissprot', 'afdb50'}, "Invalid target database choice. Possible choices: 'afdb-proteome', 'afdb-swissprot' and 'afdb50'."
    assert args.query_folder.exists() and args.query_folder.is_dir() and any(args.query_folder.iterdir()), 'Query folder path does not exist or is not a non-empty directory.'
    assert args.cores > 0, 'Number of cores must be strictly positive.'
    assert args.max_workers > 0, 'Number of workers must be positive.'
    assert args.max_eval <= 1 and args.max_eval > 0, 'Maximum e-value should be a number between 0 and 1.'
    assert args.min_seqid >= 0 and args.min_seqid <= 100, "Minimum sequence identity should be a percentage between 0 and 100."
    assert args.min_score >= 0, "Minimum FoldSeek bitscore should be a positive number."
    assert args.min_qcov >= 0 and args.min_qcov <= 100, "Minimum query coverage should be a percentage between 0 and 100."
    assert args.min_tcov >= 0 and args.min_tcov <= 100, "Minimum target coverage should be a percentage between 0 and 100."
    assert args.max_gap >= 0, "Maximum intergenic gap should be a positive number."
    assert args.max_length >= 1, "Maximum cluster length should be strictly positive."
    assert args.min_hits >= 1, "Minimum number of hits in a cluster should be strictly positive."
    assert args.min_cov_qrs >= 1, "Minimum number of covered queries in a cluster should be strictly positive."
    if args.mode == 'remote':
        db = args.db
        assert args.mapping_table_path.exists() and args.mapping_table_path.is_file(), "UniProt mapping table path does not exist or is not a file."
    elif args.mode == 'local':
        db = ["local"]
        assert args.local_db_path.exists(), "Local FoldSeek DB does not seem to exist."
        assert args.cds_mapping_table_path.exists() and args.cds_mapping_table_path.is_file(), "CDS mapping table path does not exist or is not a file."
    
    # Configure the logger
    log_levels = {0: logging.CRITICAL,
                  1: logging.ERROR,
                  2: logging.WARNING,
                  3: logging.INFO,
                  4: logging.DEBUG
                  }
    logging.basicConfig(
        level = log_levels[args.verbosity],
        format = "[%(asctime)s] %(levelname)s [%(filename)s: %(funcName)s] - %(message)s",
        datefmt="%H:%M:%S",
        handlers = [logging.StreamHandler(sys.stdout)]
        )
    
    # Parse the arguments
    params = {'mode': args.mode,
              'cores': args.cores,
              'verbosity': args.verbosity,
              'max_workers': args.max_workers,
              'max_eval': args.max_eval,
              'min_score': args.min_score,
              'min_seqid': args.min_seqid,
              'min_qcov': args.min_qcov,
              'min_tcov': args.min_tcov,
              'max_gap': args.max_gap,
              'max_length': args.max_length,
              'min_hits': args.min_hits,
              'min_cov_qrs': args.min_cov_qrs,
              'require': args.require,
              'db': db,
              'taxfilters': args.taxfilters
              }
    
    paths = {'query' : {q.stem: q.resolve() for q in args.query_folder.iterdir()},
             'uniprot_mapping' : args.mapping_table_path.resolve(),
             'cds_mapping' : args.cds_mapping_table_path.resolve(),
             'local_db_path' : args.local_db_path.resolve(),
             'output_folder' : args.output.resolve(),
             'temp_folder': args.temp.resolve()
             }
    try:
        paths['output_folder'].mkdir(parents = True)
    except FileExistsError:
        if args.force:
            LOG.warning('Output folder already exists, but it will be overwritten.')
        else:
            LOG.error('Output folder already exists! Rerun with -f to overwrite it.')
            sys.exit()
    if str(paths['temp_folder']) != tempfile.gettempdir():
        try:
            paths['temp_folder'].mkdir(parents = True)
        except FileExistsError:
            if args.force:
                LOG.warning('Temporary folder already exists, but it will be overwritten.')
            else:
                LOG.error('Temporary folder already exists! Rerun with -f to overwrite it.')
                sys.exit()
    paths['temp_folder'] = Path(tempfile.mkdtemp(dir = paths['temp_folder']))
    
    
    output_flags = {'tables': args.output_tables,
                    'session': args.output_session,
                    'summary': args.output_summary,
                    'binary': args.output_binary,
                    'plot': args.output_plot,
                    'clinker': args.output_clinker,
                    'foldseek': args.output_foldseek
                    }
    
    parsed_args = {'params': params,
                   'paths': paths,
                   'output_flags': output_flags
                   }
    
    return parsed_args
    

def main():
    # First we parse the arguments:
    args = getArguments()
    parsed_args = parseArguments(args)
    
    # Then we initiate the right workflow
    if parsed_args['params']['mode'] == 'remote':
        LOG.info("Launching cfoldseeker in remote mode")
        the_run = RemoteSearch(query = parsed_args['paths']['query'],
                               mapping_table_path = parsed_args['paths']['uniprot_mapping'],
                               params = parsed_args['params'],
                               output_folder = parsed_args['paths']['output_folder'],
                               temp_folder = parsed_args['paths']['temp_folder']
                               )
    elif parsed_args['params']['mode'] == 'local':
        LOG.info("Launching cfoldseeker in local mode")
        the_run = LocalSearch(query = parsed_args['paths']['query'],
                              db_path = parsed_args['paths']['local_db_path'],
                              coord_db_path = parsed_args['paths']['cds_mapping'],
                              params = parsed_args['params'],
                              output_folder = parsed_args['paths']['output_folder'],
                              temp_folder = parsed_args['paths']['temp_folder']
                              )
    else:
        LOG.critical("Invalid search mode!")
        LOG.critical(f"Search mode requested: {parsed_args['params']['mode']}")
        sys.exit()
    
    # Run the workflow    
    LOG.info("STARTING SEARCH")
    the_run.run()
    
    # Generate requested output
    if any([args.output_binary, args.output_clinker, args.output_plot, args.output_summary, args.output_session]):
        LOG.info("Generating cblaster session")
        cblaster_session = the_run.generate_cblaster_session()
        
        if args.output_session:
            LOG.info("Writing cblaster session file")
            path = the_run.OUTPUT_DIR / "session.json"
            with open(path, "w") as handle:
                cblaster_session.to_json(fp = handle)
            LOG.debug(f'cblaster session file written at {str(path)}')
        
        if args.output_summary:
            LOG.info("Writing cblaster summary file")
            path = the_run.OUTPUT_DIR / 'summary.txt'
            with open(path, 'w') as handle:
                cblaster_session.format(form = "summary", fp = handle)
            LOG.debug(f'cblaster summary file written at {str(path)}')
            
        if args.output_binary:
            LOG.info("Writing cblaster binary table")
            path = the_run.OUTPUT_DIR / 'binary.txt'
            with open(path, 'w') as handle:
                cblaster_session.format(form = "binary", fp = handle, delimiter = "\t")
            LOG.debug(f'cblaster binary table written at {str(path)}')
        
        if args.output_plot:
            LOG.info("Writing cblaster plot")
            path = the_run.OUTPUT_DIR / 'plot.html'
            plot_session(cblaster_session, output = path)
            LOG.debug(f'cblaster plot written at {str(path)}')
        
        if args.output_clinker:
            LOG.info("Writing clinker plot")
            path = the_run.OUTPUT_DIR / "clinker.html"
            with open(the_run.TEMP_DIR / "session.json", "w") as handle:
                cblaster_session.to_json(fp = handle)
            plot_clusters(the_run.TEMP_DIR / "session.json", plot_outfile = path, max_clusters = 10**6)
            LOG.debug(f'clinker plot written at {str(path)}')
        
    if args.output_foldseek:
        LOG.info("Copying FoldSeek output")
        for file in the_run.TEMP_DIR.glob('foldseek_result*'):
            shutil.copy(file, the_run.OUTPUT_DIR / file.name)
        LOG.debug(f'FoldSeek output copied to {the_run.OUTPUT_DIR}')
    
    if args.output_tables:
        LOG.info("Writing output tables")
        the_run.generate_tables(the_run.OUTPUT_DIR)
        LOG.debug(f'Output tables written to {the_run.OUTPUT_DIR}')
    
    # Clean up temporary directory
    LOG.info("Cleaning up temporary files")
    shutil.rmtree(the_run.TEMP_DIR)
    LOG.debug("Temporary files have been removed")
    
    
    LOG.info('DONE!')


if __name__ == "__main__":
    main()
