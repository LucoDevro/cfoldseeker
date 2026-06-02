#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import polars as pl
import gzip
import logging
import sys
import argparse
from pathlib import Path
from tqdm.contrib.concurrent import thread_map
from tqdm.contrib.logging import logging_redirect_tqdm
from cblaster.classes import Session, Scaffold
from cblaster.extract_clusters import get_sorted_cluster_hierarchies, cluster_to_record
from Bio import SeqIO


LOG = logging.getLogger(__name__)

logging.basicConfig(
    level = logging.INFO,
    format = "[%(asctime)s] %(levelname)s [%(filename)s: %(funcName)s] - %(message)s",
    datefmt="%H:%M:%S",
    handlers = [logging.StreamHandler(sys.stdout)]
    )


def create_parser() -> argparse.ArgumentParser:
    """
    This function creates a parser object that will collect the arguments given through the command line.
    
    Args:
        None
    
    Returns:
        parser (ArgumentParser): An ArgumentParser object holding the CLI ready to collect the arguments when called
    """
    
    parser = argparse.ArgumentParser(
        prog = 'cfoldseeker-seqs',
                epilog = 
                """
                Lucas De Vrieze
                (c) 2026 Masschelein lab, VIB
                """,
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = 
                """
                Helper tool to extract gene cluster Genbanks from a cfoldseeker search session
                """,
                add_help = False
                )
    args_general = parser.add_argument_group('General')
    args_general.add_argument('-c', '--cores', dest = 'cores', type = int, default = 1,
                        help = 'Number of parallel workers (default: 1).')
    args_general.add_argument('-f', '--force', dest = 'force', default = False, action = 'store_true', 
                        help = "Force overwriting output (default: false).")
    args_general.add_argument('-np', '--no-progress', dest = 'no_progress', default = False, action = "store_true", 
                        help = "Don't show progress bar (default: False).")
    args_general.add_argument('-vv', '--verbosity', dest = 'verbosity', default = 3, type = int, choices = [0,1,2,3,4],
                        help = "Console verbosity level (default: 3 (info))")
    args_general.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit")    
    
    args_io = parser.add_argument_group('File inputs and outputs')
    args_io.add_argument('-s', '--session', dest = "session", type = Path, required = True,
                        help = "Path to cfoldseeker session file.")
    args_io.add_argument('-o', '--output', dest = 'output_dir', type = Path, default = Path('.'),
                        help = 'Path to output folder (default: current workdir).')
    args_io.add_argument('-fna', '--nucleotide-fasta', dest = 'nucl_fastas_path', type = Path, required = True,
                        help = 'Path to folder with genomic nucleotide fasta files.')
    args_io.add_argument('-faa', '--protein-fasta', dest = 'prot_fastas_path', type = Path, required = True,
                        help = "Path to folder with genomic protein fasta files.")
    args_io.add_argument('--prefix', dest = 'prefix', type = str, default = '',
                         help = "String to start the file name of each cluster with (default: '').")
    args_io.add_argument('--flavour', dest = 'flavour', type = str, choices = ['genbank', 'bigscape'], default = 'genbank',
                         help = 'The flavour that the extracted cluster genbank should have (choices: genbank, bigscape) (default: genbank).')
    
    args_filt = parser.add_argument_group('Cluster filters')
    args_filt.add_argument('--cluster-numbers', dest = 'cluster_numbers', type = str, nargs = '*', default = None,
                           help = "cluster numbers to include.")
    args_filt.add_argument("--score-threshold", dest = "score_threshold", type = float, default = None,
                           help = "minimum score for a cluster to be included")
    args_filt.add_argument('--organisms', dest = "organisms", type = str, nargs = '*', default = None,
                           help = "Organism filtering regular expressions. Clusters for these organisms are included.")
    args_filt.add_argument("--scaffolds", dest = "scaffolds", type = str, nargs = '*', default = None,
                           help = "Clusters on these scaffolds are included.")
    args_filt.add_argument('-mc', '--max-clusters', dest = 'max_clusters', type = int, default = None,
                           help = "The maximum number of clusters extracted regardless of filters.")  
    
    return parser


def setup_logging(verbosity: int) -> None:
    """
    Set up the root logger if it has not been set up yet.
    
    Args:
        verbosity (int): Verbosity level (choices: 0,1,2,3,4).
        
    Returns:
        None
    """
    root_logger = logging.getLogger()
    if root_logger.handlers:
        return None
    
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
        handlers = [logging.StreamHandler(sys.stdout),],
        force = True
        )
    
    return None


def parse_and_validate_arguments(args: argparse.Namespace) -> dict:
    """
    This function parses and validates the arguments received through the command line.
    
    Args:
        args (argparse.Namespace): A Namespace holder object with the parsed argument values
    
    Returns:
        parsed_args (dict): A dictionary holding the parsed and validated argument values.
        
    Raises:
        ValueError: if an invalid argument value was given.
    """
    # Validate arguments
    try:
        if not(args.cores > 0):
            raise ValueError('Number of workers must be a strictly positive integer.')
        if args.max_clusters and not(int(args.max_clusters) > 0):
            raise ValueError('Maximum number of clusters must be a strictly positive integer.')
        if args.score_threshold and not(float(args.score_threshold) > 0):
            raise ValueError('Score threshold must be a strictly positive integer.')
        if args.output_dir.is_dir():
            if args.force:
                LOG.warning("Output folder already exists, but it will be overwritten.")
            else:
                raise FileExistsError("Output folder already exists! Rerun with -f to overwrite it.")
        else:
            args.output_dir.mkdir(parents = True, exist_ok = True)
        if not args.session.is_file():
            raise IOError('Session file does not exist.')
        if not args.nucl_fastas_path.is_dir():
            raise IOError('Nucleotide fasta folder does not exist.')
        if not args.prot_fastas_path.is_dir():
            raise IOError('Protein fasta folder does not exist.')
        
        filelabels_nucl = {'.'.join(p.with_suffix('.gz').stem.split('.')[:-1]) for p in args.nucl_fastas_path.iterdir()}
        filelabels_prot = {p.stem for p in args.prot_fastas_path.iterdir()}
        # At least one full overlap is sufficient
        if not (filelabels_nucl <= filelabels_prot or filelabels_nucl >= filelabels_prot):
            raise IOError('Filelabels in nucleotide and protein fasta folders do not fully overlap.')
            
    except (ValueError, IOError) as err:
        LOG.critical(err)
        raise err
    
    # Convert validated arguments to dictionary
    parsed_args = vars(args)
    # Group cluster filter arguments
    filter_names = ['cluster_numbers', 'score_threshold', 'organisms', 'scaffolds', 'max_clusters']
    cluster_filters = {filt: parsed_args[filt] for filt in filter_names}
    parsed_args = {k:v for k,v in parsed_args.items() if k not in filter_names}
    parsed_args['cluster_filters'] = cluster_filters
    
    return parsed_args


def locate_nucleotide_sequences(scaffolds: list[Scaffold]) -> pl.DataFrame:
    """
    Locate all nucleotide sequences to be fetched from session information.
    
    Deduces where to get the nucleotides sequences for every cluster in the session file
    from location information saved in there.
    
    Args:
        scaffolds (list): List of cblaster Scaffold objects sourced from the session file
        
    Returns:
        df (polars.DataFrame): DataFrame holding the file location and genomic range
            of every cluster in the session file.
    """
    nuc_locations = []
    for scaff in scaffolds:
        for cluster in scaff.clusters:
            record = {'scaffold': scaff.accession,
                      'cluster_number': cluster.number,
                      'start': cluster.start,
                      'end': cluster.end,
                      'assembly': cluster.subjects[0].sequence}
            nuc_locations.append(record)
    
    return pl.from_dicts(nuc_locations)


def locate_protein_sequences(scaffolds: list[Scaffold]) -> pl.DataFrame:
    """
    Locate all protein sequences to be fetched from session information.
    
    Deduces where to get the protein sequences for every cluster in the session file
    from location information saved in there.
    
    Args:
        scaffolds (list): List of cblaster Scaffold objects sourced from the session file
        
    Returns:
        df (polars.DataFrame): DataFrame holding the file locations and protein IDs
            of every protein of each cluster in the session file.
    """
    prot_locations = []
    for scaff in scaffolds:
        for cluster in scaff.clusters:
            for protein in cluster.subjects:
                record = {'cluster_number': cluster.number,
                          'id': protein.name,
                          'assembly': protein.sequence}
                prot_locations.append(record)
    
    return pl.from_dicts(prot_locations)


def _write_one_cluster_genbank(enumerated_scaffold: tuple[int, Scaffold], assemblies: list[str], 
                               prefix: str, flavour: str,
                               nucl_locations: pl.DataFrame, prot_locations: pl.DataFrame,
                               output_dir: Path, nucl_fastas_path: Path, prot_fastas_path: Path,
                               required_genes: list[str]) -> None:
    """
    Write a Genbank cluster file for each cluster on one scaffold.
    
    Fetches the nucleotide and protein sequences for each cluster from the earlier
        located files, and writes them away in a new Genbank file.
        
    Args:
        enumerated_scaffold (tuple): tuple of an integer index and a cblaster Scaffold object
            sourced from the session file. The parent function `write_cluster_genbanks` requires
            the index to keep track of its progress while parallising.
        assemblies (list): List of assembly filenames, sourced from the session.
        prefix (str): String to start the file name of each cluster with.
        flavour (str): Requested flavour of the resulting genbank file. Either regular genbank ('genbank'), or
            a genbank with BigScape-specific fields ('bigscape').
        nucl_locations (polars.DataFrame): DataFrame with the file location of all
            nucleotide sequences to be fetched, as determined by `locate_nucleotide_sequences`.
        prot_locations (polars.DataFrame): DataFrame with the file location of all
            protein sequences to be fetched, as determined by `locate_protein_sequences`.
        output_dir (pathlib.Path): Path of the output folder.
        nucl_fastas_path (pathlib.Path): Path of the nucleotide fastas folder.
        prot_fastas_path (pathlib.Path): Path of the protein fastas folder.
        required_genes (list): List of strings corresponding with the query genes
            that were marked required during the search. Sourced from the session.
    
    Returns:
        None
    """

    idx = enumerated_scaffold[0]
    scaffold = enumerated_scaffold[1]
    
    for cluster in scaffold.clusters:
        # Pinpoint exact nucleotide location
        nuc_location = nucl_locations.filter(pl.col('cluster_number') == cluster.number).to_dicts()[0]
        nuc_file = list(nucl_fastas_path.glob(f'{nuc_location["assembly"]}*'))
        if len(nuc_file) > 1:
            msg = f'Found more than one genome file with this filelabel: {nuc_file}'
            LOG.error(msg)
            raise RuntimeError(msg)
        nuc_file = nuc_file[0]
        
        # Fetch cluster nucleotide sequence from genome fasta
        if '.gz' in nuc_file.suffixes:
            handle = gzip.open(nuc_file, 'rt')
        else:
            handle = open(nuc_file, 'r')
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id == nuc_location['scaffold']:
                nuc_seq = str(record[nuc_location['start']-1 : nuc_location['end']].seq)
                break
        handle.close()
                
        # Pinpoint exact protein location
        prot_location = prot_locations.filter(pl.col('cluster_number') == cluster.number).to_dicts()
        prot_file = list(prot_fastas_path.glob(f'{prot_location[0]["assembly"]}*'))
        prot_ids = {prot['id'] for prot in prot_location}
        if len(prot_file) > 1:
            msg = f'Found more than one protein file with this filelabel: {prot_file}'
            LOG.error(msg)
            raise RuntimeError(msg)
        if len(prot_file) == 0:
            LOG.warning(f'No protein file found for {prot_location[0]["assembly"]}! Skipping...')
            continue
        prot_file = prot_file[0]
        
        # Fetch cluster protein sequences from protein fasta
        with open(prot_file, 'r') as handle:
            cluster_prot_sequences = {record.id: str(record.seq) 
                                      for record in SeqIO.parse(handle, 'fasta') 
                                      if record.id in prot_ids}
        
        # Define sequence Record
        seq_record = cluster_to_record(cluster = cluster,
                                       cluster_prot_sequences = cluster_prot_sequences,
                                       cluster_nuc_sequence = nuc_seq,
                                       organism_name = assemblies[idx],
                                       scaffold_accession = scaffold.accession,
                                       format_ = flavour,
                                       required_genes = required_genes)
        
        # Write Genbank file
        output_file = output_dir / f"{prefix}cluster{cluster.number}.gbk"
        LOG.debug(f'Writing Genbank file {output_file.name}')
        with open(output_file, 'w') as handle:
            SeqIO.write(seq_record, handle, 'genbank')
                
    return None


def write_cluster_genbanks(scaffolds: list[Scaffold], assemblies: list[str], 
                           prefix: str, flavour: str, required_genes: list[str],
                           nucl_locations: pl.DataFrame, prot_locations: pl.DataFrame,
                           output_dir: Path, nucl_fastas_path: Path, prot_fastas_path: Path,
                           n_workers: int = 1, no_progress: bool = False) -> None:
    """
    Write Genbank cluster files for each cluster on all scaffolds in the session.
    
    Fetches the nucleotide and protein sequences for each cluster from the earlier
        located files, and writes them away in new Genbank files. Supports
        parallellisation.
        
    Args:
        scaffolds (list): List of cblaster Scaffold objects sourced from the session file.
        assemblies (list): List of assembly filenames, sourced from the session.
        prefix (str): String to start the file name of each cluster with.
        flavour (str): Requested flavour of the resulting genbank file. Either regular genbank ('genbank'), or
            a genbank with BigScape-specific fields ('bigscape').
        required_genes (list): List of strings corresponding with the query genes
            that were marked required during the search. Sourced from the session.
        nucl_locations (polars.DataFrame): DataFrame with the file location of all
            nucleotide sequences to be fetched, as determined by `locate_nucleotide_sequences`.
        prot_locations (polars.DataFrame): DataFrame with the file location of all
            protein sequences to be fetched, as determined by `locate_protein_sequences`.
        output_dir (pathlib.Path): Path of the output folder.
        nucl_fastas_path (pathlib.Path): Path of the nucleotide fastas folder.
        prot_fastas_path (pathlib.Path): Path of the protein fastas folder.
        n_workers (int): Number of parallel worker threads. Defaults to 1.
        no_progress (bool): Flag to disable showing the progress bar. Defaults to False.
    
    Returns:
        None
    """
    
    with logging_redirect_tqdm(loggers = [LOG]):
        thread_map(lambda x: _write_one_cluster_genbank(enumerated_scaffold = x,
                                                        assemblies = assemblies,
                                                        prefix = prefix,
                                                        flavour = flavour,
                                                        nucl_locations = nucl_locations, 
                                                        prot_locations = prot_locations,
                                                        output_dir = output_dir, 
                                                        nucl_fastas_path = nucl_fastas_path, 
                                                        prot_fastas_path = prot_fastas_path,
                                                        required_genes = required_genes),
                   list(enumerate(scaffolds)),
                   max_workers = n_workers,
                   leave = False,
                   disable = no_progress)
    
    return None


def select_clusters(session: Session, filters: dict) -> tuple:
    """
    Selects the clusters to process based on a selection of filtering parameters.
    
    Returns the cblaster cluster hierarchy filtered by a maximum number, cluster number, score, organism or scaffold.
    
    Args:
        session (cblaster.Session): The cblaster session to get the hierarchy from.
        filters (dict): dictionary of filtering parameter values.
        
    Returns:
        cluster_hierarchy (tuple): Tuple of three checked out groups of objects: clusters, scaffolds and assemblies.
        
    Note:
        This is a wrapper for cblaster's `get_sorted_cluster_hierarchies`.
    """
    cluster_hierarchy = get_sorted_cluster_hierarchies(session = session, **filters)
    
    return cluster_hierarchy


def convert_session(session: Session) -> Session:
    """
    Convert a cblaster session in-memory to a cfoldseeker session.
    
    Adds file labels at the sequence fields of each subject in the session, which
    is necessary to retrieve the correct local genome and proteome files.
    
    Args:
        session (cblaster.Session): The cblaster session to convert.
        
    Returns:
        session (cblaster.Session): The converted session.
        
    Mutates:
        session (cblaster.Session): The converted session.
    """
    for organism in session.organisms:
        for scaff in organism.scaffolds.values():
            for subject in scaff.subjects:
                subject.sequence = organism.name
                subject.id = None
                
    return session


def run_workflow(parsed_args: dict) -> None:
    """
    Execute the sequence export workflow.
    
    Loads the session, locates the nucleotide and protein sequences in the fasta folders,
    and writes the Genbank files.
    
    Args:
        parsed_args (dict): A dictionary holding the parsed and validated argument values.
        
    Returns:
        None
    """
    session = Session.from_file(parsed_args['session'])
    
    ## Convert a cblaster session to a cfoldseeker session
    # Recognise which tool generated this session
    # cblaster leaves the sequence attribute of local searches empty; cfoldseeker fills it with the local filelabel
    first_scaffold = list(session.organisms[0].scaffolds.values())[0]
    match first_scaffold.subjects[0].sequence:
        # cblaster
        case None:
            LOG.info("Detected cblaster session. Converting to cfoldseeker session.")
            session = convert_session(session)
        # cfoldseeker
        case str():
            pass
        case _:
            raise ValueError('Could not determine session type!')
    
    LOG.info('Selecting clusters')
    clusters, scaffolds, assemblies = zip(*select_clusters(session, parsed_args['cluster_filters']))
    
    LOG.info("Locating sequences")
    nucl_locations = locate_nucleotide_sequences(scaffolds)
    prot_locations = locate_protein_sequences(scaffolds)
    
    LOG.info("Fetching sequences")
    write_cluster_genbanks(scaffolds = scaffolds,
                           assemblies = assemblies,
                           prefix = parsed_args['prefix'],
                           flavour = parsed_args['flavour'],
                           required_genes = session.params['require'],
                           nucl_locations = nucl_locations,
                           prot_locations = prot_locations,
                           output_dir = parsed_args['output_dir'],
                           nucl_fastas_path = parsed_args['nucl_fastas_path'],
                           prot_fastas_path = parsed_args['prot_fastas_path'],
                           n_workers = parsed_args['cores'],
                           no_progress = parsed_args['no_progress'])
    
    return None


def main():
    """
    Main entry point for the sequence export tool.
    
    Oversees the complete workflow: parses command-line arguments, and calls
    the workflow.
    """
    # Process arguments
    parser = create_parser()
    # Parse arguments
    args = parser.parse_args()
    # Validate arguments
    parsed_args = parse_and_validate_arguments(args)

    # Set up logging
    setup_logging(args.verbosity)
    
    # Run workflow
    run_workflow(parsed_args)
    
    LOG.info('DONE')
    

if __name__ == "__main__":
    main()
    
    
