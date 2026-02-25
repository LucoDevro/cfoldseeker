#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from pathlib import Path
from importlib.metadata import version
from importlib import resources

from remote import RemoteSearch
from local import LocalSearch

# __version__ = version("cfoldseeker")


def getArguments():
    """
    This function gets the CLI arguments, without any additional parsing.
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
                cfoldseeker: find clusters of colocalised genes via structural similarity using FoldSeek
                """,
                add_help = False
                )
    
    args_general = parser.add_argument_group('General')
    args_general.add_argument('-m', '--mode', dest = 'mode', default = 'remote', type = str, choices = ['local', 'remote'], help = "Search mode (default: remote)")
    args_general.add_argument('-c', '--cores', dest = 'cores', default = 1, type = int, help = "Number of cores available to use (default: 1)")    
    # args_general.add_argument('-v', '--version', action = "version", version = "%(prog)s " + __version__)
    args_general.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit")      
    args_general.add_argument('--verbose', dest = 'verbose', default = False, action = 'store_true', help = "Enable verbose logging")
    
    args_io = parser.add_argument_group('Inputs and outputs')
    args_io.add_argument('-q', '--query', dest = 'query_folder', required = True, type = Path, help = "Path of the folder containing the query proteins.")
    args_io.add_argument('-o', '--output', dest = 'output', type = Path, default = Path('.'), help = "Output directory (default: current location)")
    
    args_search = parser.add_argument_group('General search options')
    args_search.add_argument('--max-eval', dest = "max_eval", type = float, default = 1e-3, help = "Maximum e-value to include a FoldSeek hit.")
    args_search.add_argument('--min-prob', dest = "min_prob", type = float, default = 0, help = "Minimum FoldSeek probability to include a hit (range: 0-1).")
    args_search.add_argument('--min-score', dest = "min_score", type = float, default = 0, help = "Minimum FoldSeek bitscore to include a hit.")
    args_search.add_argument('--min-seqid', dest = "min_seqid", type = float, default = 0, help = "Minimum sequence identity to include a hit (in percentages).")
    args_search.add_argument('--min-qcov', dest = "min_qcov", type = float, default = 0, help = "Minimum query coverage to include a hit (in percentages).")
    args_search.add_argument('--min-tcov', dest = "min_tcov", type = float, default = 0, help = "Minimum target coverage to include a hit (in percentages).")
    args_search.add_argument('--max-gap', dest = "max_gap", type = int, default = 5000, help = "Maximum intergenic gap within a cluster (in bp).")
    args_search.add_argument('--max-length', dest = "max_length", type = int, default = 100000, help = "Maximum genomic length of a cluster (in bp).")
    args_search.add_argument('--min-hits', dest = "min_hits", type = int, default = 0, help = "Minimum number of members in a cluster.")
    args_search.add_argument('--min-cov-qrs', dest = "min_cov_qrs", type = int, default = 0, help = "Minimum different queries covered by a cluster.")
    args_search.add_argument('--require', dest = "require", type = str, default = '', nargs = '*', help = "Queries that have to present in a cluster (use filename stems as labels).")
    
    args_remote = parser.add_argument_group("Remote-specific search options")
    args_remote.add_argument('-db', '--database', dest = 'db', type = str, default = 'local', nargs = '*', choices = ['afdb-proteome', 'afdb-swissprot', 'afdb50'], 
                             help = "Target database to include (choices: afdb-proteome, afdb-swissprot, afdb50)")
    # args_remote.add_argument('-uma', '--uniprot-mapping', dest = 'mapping_table_path', type = Path, default = Path(resources.files(__name__)) / 'uniprot_kegg_genpept.gz',
    args_remote.add_argument('-uma', '--uniprot-mapping', dest = 'mapping_table_path', type = Path, default = 'uniprot_kegg_genpept.gz',
                             help = "Path of the UniProt AFDB ID mapping table (default: use builtin one)")
    args_remote.add_argument('--max-workers', dest = "max_workers", type = int, default = 1, help = "Maximum number of workers to query the remote servers (FoldSeek, KEGG, ENA) (default: 1)")
    
    args_local = parser.add_argument_group('Local-specific search options')
    args_local.add_argument('-ldb', '--local-database', dest = 'local_db', type = Path, default = Path('local_db'), help = "Path to your local protein structures folder (default: local_db).")
    args_local.add_argument('-cm', '--cds-mapping', dest = 'cds_mapping_table_path', type = Path, default = Path('local_db_cds.tsv'), help = "Path of the CDS coordinates mapping table of your local database (default: local_db_cds.tsv).")
    args_local.add_argument('--fs-db-prefix', dest = 'fs_db_prefix', type = str, default = 'local_db', help = "Prefix of the FoldSeek DB files that will be constructed from your local database (default: local_db).")
    args_local.add_argument('--foldseek-temp', dest = 'foldseek_temp', type = Path, default = Path('./tmp'), help = "Location that FoldSeek may use for its temporary files (default: tmp subfolder in the current folder)")
    
    args = parser.parse_args()
    
    return args


def parseArguments(args):
    """
    This function parses the arguments and returns a parsed arguments dictionary that can be used to call the workflows.
    """
    
    assert args.mode in ['local', 'remote'], 'Invalid search mode. Possible choices: "local" and "remote".'
    assert set(args.db) <= {'local', 'afdb-proteome', 'afdb-swissprot', 'afdb50'}, "Invalid target database choice. Possible choices: 'afdb-proteome', 'afdb-swissprot' and 'afdb50'."
    assert args.query_folder.exists() and args.query_folder.is_dir() and any(args.query_folder.iterdir()), 'Query folder path does not exist or is not a non-empty directory.'
    assert args.cores > 0, 'Number of cores must be positive.'
    assert args.max_workers > 0, 'Number of workers must be positive.'
    assert args.max_eval < 1 and args.max_eval > 0, 'Maximum e-value should be a number between 0 and 1.'
    assert args.min_seqid >= 0 and args.min_seqid <= 100, "Minimum sequence identity should be a percentage between 0 and 100."
    assert args.min_prob >= 0 and args.min_prob <= 1, "Minimum hit probability score should be a number between 0 and 1."
    assert args.min_score >= 0, "Minimum FoldSeek bitscore should be a positive number."
    assert args.min_qcov >= 0 and args.min_qcov <= 100, "Minimum query coverage should be a percentage between 0 and 100."
    assert args.min_tcov >= 0 and args.min_tcov <= 100, "Minimum target coverage should be a percentage between 0 and 100."
    assert args.max_gap >= 0, "Maximum intergenic gap should be a positive number."
    assert args.max_length >= 0, "Maximum cluster length should be a positive number."
    assert args.min_hits >= 0, "Minimum number of hits in a cluster should be a positive number."
    assert args.min_cov_qrs >= 0, "Minimum number of covered queries in a cluster should be a positive number."
    if args.mode == 'remote':
        assert args.mapping_table_path.exists() and args.mapping_table_path.is_file(), "UniProt mapping table path does not exist or is not a file."
    elif args.mode == 'local':
        assert args.local_db.exists() and args.local_db.is_dir() and any(args.local_db.iterdir()), "Local structure folder does not exist or is not a non-empty directory."
        assert args.cds_mapping_table_path.exists() and args.cds_mapping_table_path.is_file(), "CDS mapping table path does not exist or is not a file."
        assert len(args.fs_db_prefix) > 0, "FoldSeek DB prefix cannot be empty."
    
    params = {'mode': args.mode,
              'cores': args.cores,
              'max_workers': args.max_workers,
              'max_eval': args.max_eval,
              'min_prob': args.min_prob,
              'min_score': args.min_score,
              'min_seqid': args.min_seqid,
              'min_qcov': args.min_qcov,
              'min_tcov': args.min_tcov,
              'max_gap': args.max_gap,
              'max_length': args.max_length,
              'min_hits': args.min_hits,
              'min_cov_qrs': args.min_cov_qrs,
              'require': args.require,
              'db': args.db
              }
    
    fs_db_prefix = args.fs_db_prefix
    
    query = {q.stem: q.resolve() for q in args.query_folder.iterdir()}
    uniprot_mapping = args.mapping_table_path.resolve()
    cds_mapping = args.cds_mapping_table_path.resolve()
    local_db = args.local_db.resolve()
    output = args.output.resolve()
    output.mkdir(parents = True, exist_ok = True)
    foldseek_temp = args.foldseek_temp.resolve()
    foldseek_temp.mkdir(parents = True, exist_ok = True)
    
    parsed_args = {'params': params,
                   'fs_db_prefix': fs_db_prefix,
                   'query': query,
                   'uniprot_mapping': uniprot_mapping,
                   'cds_mapping': cds_mapping,
                   'output': output,
                   'foldseek_temp': foldseek_temp,
                   'local_db': local_db}
    
    return parsed_args
    

def main():
    # First we parse the arguments:
    args = getArguments()
    parsed_args = parseArguments(args)
    
    # Then we initiate the right workflow
    if parsed_args['params']['mode'] == 'remote':
        the_run = RemoteSearch(query = parsed_args['query'],
                               mapping_table_path = parsed_args['uniprot_mapping'],
                               params = parsed_args['params'])
    elif parsed_args['params']['mode'] == 'local':
        the_run = LocalSearch(query = parsed_args['query'],
                              db_path = parsed_args['local_db'],
                              coord_db_path = parsed_args['cds_mapping'],
                              db_prefix = parsed_args['fs_db_prefix'],
                              params = parsed_args['params'])
    else:
        sys.exit("Invalid search mode!")
    
    # Run the workflow
    the_run.run()
    

if __name__ == "__main__":
    main()
