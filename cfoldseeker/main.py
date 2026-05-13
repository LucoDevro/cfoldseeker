#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import shutil
import logging
import tempfile
from pathlib import Path
from importlib.metadata import version

from cfoldseeker.remote import RemoteSearch
from cfoldseeker.local import LocalSearch
from cfoldseeker.local_clustered import LocalClusteredSearch

__version__ = version("cfoldseeker")

# Setup default logger configuration
logging.basicConfig(
    level = logging.INFO,
    format = "[%(asctime)s] %(levelname)s [%(filename)s: %(funcName)s] - %(message)s",
    datefmt="%H:%M:%S",
    handlers = [logging.StreamHandler(sys.stdout)],
    )
LOG = logging.getLogger(__name__)


def create_parser() -> argparse.ArgumentParser:
    """
    This function creates a parser object that will collect the arguments given through the command line.
    
    Args:
        None
    
    Returns:
        parser (argparse.ArgumentParser): An ArgumentParser object holding the CLI ready to collect the arguments when called
    """
    
    parser = argparse.ArgumentParser(
        prog = 'cfoldseeker',
                epilog = 
                """
                Lucas De Vrieze
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
    args_general.add_argument('-m', '--mode', dest = 'mode', default = 'remote', type = str, choices = ['local', 'local_clustered', 'remote'], help = "Search mode (default: remote)")
    args_general.add_argument('-c', '--cores', dest = 'cores', default = 1, type = int, help = "Number of cores available to use (default: 1)")          
    args_general.add_argument('-f', '--force', dest = 'force', default = False, action = 'store_true', help = "Force overwriting output (default: false).")
    args_general.add_argument('-vv', '--verbosity', dest = 'verbosity', default = 3, type = int, choices = [0,1,2,3,4], help = "Console verbosity level (default: 3 (info))")
    args_general.add_argument('-np', '--no-progress', dest = 'no_progress', default = False, action = 'store_true', help = "Don't show progress bar (default: False).")
    args_general.add_argument('-v', '--version', action = "version", version = "%(prog)s " + __version__)
    args_general.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit")
    
    args_io = parser.add_argument_group('Inputs and outputs')
    args_io.add_argument('-q', '--query', dest = 'query_folder', required = True, type = Path, help = "Path of the folder containing the query proteins.")
    args_io.add_argument('-o', '--output', dest = 'output', type = Path, default = Path('.'), help = "Output directory (default: current location)")
    args_io.add_argument('-t', '--temp', dest = "temp", type = Path, default = Path(tempfile.gettempdir()), help = "Path to store temporary files (default: your OS's default temporary directory).")
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
    args_search.add_argument('--all-layouts', dest = 'all_layouts', default = False, action = 'store_true', help = "Include cluster layouts matching all criteria but with a less-than-best score (default: False).")
    
    args_remote = parser.add_argument_group("Remote-specific search options")
    args_remote.add_argument('-rdb', '--remote-db', dest = 'db', type = str, default = ['afdb50'], nargs = '+', choices = ['afdb-proteome', 'afdb-swissprot', 'afdb50'], 
                             help = "Remote target database (default: afdb50) (choices: afdb-proteome, afdb-swissprot, afdb50)")
    args_remote.add_argument('-tf', '--taxon-filter', dest = 'taxfilters', type = str, default = '', nargs = '*',
                             help = "Taxon ID(s) to filter the FoldSeek results table.")
    args_remote.add_argument('-uma', '--uniprot-mapping', dest = 'mapping_table_path', type = Path, default = Path('uniprot_kegg_genpept.gz'),
                             help = "Path to the UniProt AFDB ID mapping table (default: uniprot_kegg_genpept.gz)")
    args_remote.add_argument('-w', '--max-workers', dest = "max_workers", type = int, default = 1, help = "Maximum number of workers to query the remote servers (FoldSeek, KEGG, ENA) (default: 1)")
    
    args_local = parser.add_argument_group('Local-specific search options')
    args_local.add_argument('-ldb', '--local-db', dest = 'local_db_path', type = Path, default = Path('local_db/local_db'), help = "Path to your local FoldSeek DB (format: <path-to-containing-folder>/<DB-prefix>) (default: local_db/local_db).")
    args_local.add_argument('-cdb', '--cds-coords-db', dest = 'cds_db_path', type = Path, default = Path('local_cds_db.gz'), help = "Path of the CDS coordinates DB (default: local_cds_db.gz).")
    
    args_local_clustered = parser.add_argument_group('Local-clustered-specific search options')
    args_local_clustered.add_argument('-scl', '--seq-clusters', dest = "seq_clusters", type = Path, default = Path('cluster_clustered.tsv'),
                                      help = "Path to MMseqs2 clustering table TSV file (default: cluster_clustered.tsv).")
    
    return parser


def setup_logging(verbosity: int) -> None:
    """
    Set up the root logger if it has not been set up yet.
    
    Args:
        verbosity (int): Verbosity level (choices: 0,1,2,3,4).
        
    Returns:
        None
    """
    log_levels = {0: logging.CRITICAL,
                  1: logging.ERROR,
                  2: logging.WARNING,
                  3: logging.INFO,
                  4: logging.DEBUG
                  }
    logging.basicConfig(
        level = log_levels[verbosity],
        format = "[%(asctime)s] %(levelname)s [%(filename)s: %(funcName)s] - %(message)s",
        datefmt="%H:%M:%S",
        handlers = [logging.StreamHandler(sys.stdout)],
        force = True
        )
    
    return None


def parse_and_validate_arguments(args: argparse.Namespace, skip_context_table_check: bool = False) -> dict:
    """
    This function validates the parsed arguments given through the command line.
    
    Args:
        parser (argparse.NameSpace): A NameSpace object with parsed CLI arguments
        skip_csuite_IO_checks (bool): Skip argument validation for intermediary inputs and outputs
            in the csuite workflows. For compatibility with the csuite validation checker.
    
    Returns:
        parsed_args (dict): A dictionary holding the parsed and validated argument values.
        
    Raises:
        ValueError: if an invalid argument value was given.
    """
    ## Validate arguments
    if args.mode not in ['local', 'remote', 'local_clustered']:
        raise ValueError('Invalid search mode. Possible choices: "local" and "remote".')
    if not(set(args.db) <= {'local', 'afdb-proteome', 'afdb-swissprot', 'afdb50'}):
        raise ValueError("Invalid target database choice. Possible choices: 'afdb-proteome', 'afdb-swissprot' and 'afdb50'.")
    if not(args.query_folder.is_dir() and any(args.query_folder.glob('*cif'))):
        raise ValueError('Query folder path does not exist or does not contain cif files.')
    if not(args.cores > 0):
        raise ValueError('Number of cores must be strictly positive.')
    if not(args.max_workers > 0): 
        raise ValueError('Number of workers must be positive.')
    if not(args.max_eval <= 1 and args.max_eval > 0): 
        raise ValueError('Maximum e-value should be a number between 0 and 1.')
    if not(args.min_seqid >= 0 and args.min_seqid <= 100):
        raise ValueError("Minimum sequence identity should be a percentage between 0 and 100.")
    if not(args.min_score >= 0):
        raise ValueError("Minimum FoldSeek bitscore should be a positive number.")
    if not(args.min_qcov >= 0 and args.min_qcov <= 100):
        raise ValueError("Minimum query coverage should be a percentage between 0 and 100.")
    if not(args.min_tcov >= 0 and args.min_tcov <= 100):
        raise ValueError("Minimum target coverage should be a percentage between 0 and 100.")
    if not(args.max_gap >= 0): 
        raise ValueError("Maximum intergenic gap should be a positive number.")
    if not(args.max_length >= 1): 
        raise ValueError("Maximum cluster length should be strictly positive.")
    if not(args.min_hits >= 1): 
        raise ValueError("Minimum number of hits in a cluster should be strictly positive.")
    if not(args.min_cov_qrs >= 1): 
        raise ValueError("Minimum number of covered queries in a cluster should be strictly positive.")
    if not(set(args.require) <= {f.stem for f in args.query_folder.glob('*cif')}):
        raise ValueError("A required query cannot be found in your query folder. Please check the filenames.")
    
    # Mode-specific checks
    match args.mode:
        case 'remote':
            db = args.db
            if not(args.mapping_table_path.is_file()):
                raise ValueError("UniProt mapping table path does not exist or is not a file.")
        case 'local':
            db = ["local"]
            if not(args.local_db_path.is_file()):
                raise ValueError("Local FoldSeek DB does not exist.")
            if not skip_context_table_check:
                if not(args.cds_db_path.is_file()):
                    raise ValueError("CDS mapping table path does not exist or is not a file.")
        case 'local_clustered':
            db = ['local']
            if not(args.local_db_path.is_file()):
                raise ValueError("Local FoldSeek DB does not exist.")
            if not(args.seq_clusters.is_file()):
                raise ValueError("MMseqs2 clustering table does not exist.")
            if not skip_context_table_check:
                if not(args.cds_db_path.is_file()):
                    raise ValueError("CDS mapping table path does not exist or is not a file.")
    
    ## Parse the arguments
    # Search parameters
    params = {'mode': args.mode,
              'cores': args.cores,
              'verbosity': args.verbosity,
              'no_progress': args.no_progress,
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
              'taxfilters': args.taxfilters,
              'all_layouts': args.all_layouts
              }
    
    # File paths
    paths = {'query' : {q.stem: q.resolve() for q in args.query_folder.glob('*.cif')},
             'uniprot_mapping' : args.mapping_table_path.resolve(),
             'cds_db_path' : args.cds_db_path.resolve(),
             'local_db_path' : args.local_db_path.resolve(),
             'output_folder' : args.output.resolve(),
             'temp_folder': args.temp.resolve(),
             'seq_clusters': args.seq_clusters.resolve()
             }
    # Check the most important paths
    try:
        paths['output_folder'].mkdir(parents = True)
    except FileExistsError as err:
        if args.force:
            LOG.warning('Output folder already exists, but it will be overwritten.')
        else:
            msg = 'Output folder already exists! Rerun with -f to overwrite it.'
            LOG.error(msg)
            raise err
            
    if paths['temp_folder'] != Path(tempfile.gettempdir()):
        try:
            paths['temp_folder'].mkdir(parents = True)
        except FileExistsError as err:
            if args.force:
                LOG.warning('Temporary folder already exists, but it will be overwritten.')
            else:
                msg = 'Temporary folder already exists! Rerun with -f to overwrite it.'
                LOG.error(msg)
                raise err
    paths['temp_folder'] = Path(tempfile.mkdtemp(dir = paths['temp_folder']))
    
    # Output flags
    output_flags = {'tables': args.output_tables,
                    'session': args.output_session,
                    'summary': args.output_summary,
                    'binary': args.output_binary,
                    'plot': args.output_plot,
                    'clinker': args.output_clinker,
                    'foldseek': args.output_foldseek
                    }
    
    # Gather all parameters
    parsed_args = {'params': params,
                   'paths': paths,
                   'output_flags': output_flags
                   }
    
    return parsed_args


def init_search(parsed_args):
    """
    Initialise the correct search class and pass it on the necessary arguments.
    
    Args:
        parsed_args (dict): nested dictionary holding the arguments as parsed by parse_and_validate_arguments
        
    Returns:
        the_run (RemoteSearch | LocalSearch | LocalClusteredSearch): A Search workflow object ready to run
    """
    match parsed_args['params']['mode']:
        case 'remote':
            LOG.info("Launching cfoldseeker in remote mode")
            the_run = RemoteSearch(query = parsed_args['paths']['query'],
                                   mapping_table_path = parsed_args['paths']['uniprot_mapping'],
                                   params = parsed_args['params'],
                                   output_folder = parsed_args['paths']['output_folder'],
                                   temp_folder = parsed_args['paths']['temp_folder'],
                                   output_flags = parsed_args['output_flags'],
                                   )
        case 'local':
            LOG.info("Launching cfoldseeker in local mode")
            the_run = LocalSearch(query = parsed_args['paths']['query'],
                                  db_path = parsed_args['paths']['local_db_path'],
                                  coord_db_path = parsed_args['paths']['cds_db_path'],
                                  params = parsed_args['params'],
                                  output_folder = parsed_args['paths']['output_folder'],
                                  temp_folder = parsed_args['paths']['temp_folder'],
                                  output_flags = parsed_args['output_flags'],
                                  )
        case 'local_clustered':
            LOG.info("Launching cfoldseeker in local-clustered mode")
            the_run = LocalClusteredSearch(query = parsed_args['paths']['query'],
                                           db_path = parsed_args['paths']['local_db_path'],
                                           coord_db_path = parsed_args['paths']['cds_db_path'],
                                           params = parsed_args['params'],
                                           output_folder = parsed_args['paths']['output_folder'],
                                           temp_folder = parsed_args['paths']['temp_folder'],
                                           seq_clust_tsv = parsed_args['paths']['seq_clusters'],
                                           output_flags = parsed_args['output_flags'],
                                           )
    
    return the_run


def run_workflow(parsed_args: dict) -> None:
    """
    Execute the complete cfoldseeker workflow.
    
    Initialises the appropriate Run instance, executes it, generates the output,
    and cleans up the temporary files.
    
    Returns:
        None
    """
    # Instantiate the right workflow type
    the_run = init_search(parsed_args)
    
    # Run the workflow    
    LOG.info("STARTING SEARCH")
    the_run.run()
    
    # Generate requested output
    the_run.generate_output()
    
    # Clean up temporary directory
    LOG.info("Cleaning up temporary files")
    shutil.rmtree(the_run.TEMP_DIR)
    LOG.debug("Temporary files have been removed")
    
    return None


def main():
    """
    Main entry point of cfoldseeker.
    
    Oversees the complete workflow: parses command-line arguments, sets up the
    logger and the run, and calls the workflow.
    """
    # First we parse and validate the arguments:
    parser = create_parser()
    args = parser.parse_args()
    parsed_args = parse_and_validate_arguments(args)
    
    # Configure the logger
    setup_logging(args.verbosity)
    
    # Set up and run the appropriate workflow
    run_workflow(parsed_args)
    
    LOG.info('DONE!')


if __name__ == "__main__":
    main()
