#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import gzip
import logging
import tempfile
import polars as pl
from Bio import SeqIO, Entrez
from pathlib import Path
from itertools import batched, chain
from csv import DictWriter
from functools import partial
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm


LOG = logging.getLogger(__name__)

logging.basicConfig(
    level = logging.INFO,
    format = "[%(asctime)s] %(levelname)s [%(filename)s: %(funcName)s] - %(message)s",
    datefmt="%H:%M:%S",
    handlers = [logging.StreamHandler(sys.stdout)]
    )


ALLOWED_EXTENSIONS = ['.gb', '.gbk', '.gbff', '.gb.gz', '.gbk.gz', '.gbff.gz']


def create_parser() -> argparse.ArgumentParser:
    """
    This function creates a parser object that will collect the arguments given through the command line.
    
    Args:
        None
    
    Returns:
        parser (ArgumentParser): An ArgumentParser object holding the CLI ready to collect the arguments when called
    """
    
    parser = argparse.ArgumentParser(
        prog = 'cfoldseeker-cds',
                epilog = 
                """
                Lucas De Vrieze
                (c) 2026 Masschelein lab, VIB
                """,
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = 
                """
                Helper tool to construct a CDS coordinates DB for cfoldseeker
                """,
                add_help = False
                )
    
    parser.add_argument('-i', '--input', dest = 'input', type = Path, default = Path('.'), 
                        help = "Path to folder holding the input files or NCBI package (default: current directory)")
    parser.add_argument('-m', '--mode', dest = 'mode', type = str, required = True, 
                        choices = ['ncbi-gbff', 'ncbi-package', 'bakta-gbff', 'tsv'],
                        help = 'File parsing mode (choices: ncbi-gbff, ncbi-package, bakta-gbff, tsv).')
    parser.add_argument('-o', '--output', dest = 'output', type = Path, default = Path('local_db'), 
                        help = "Filepath to save CDS coordinate DB (default: local_db).")
    parser.add_argument('-t', '--temp', dest = "temp", type = Path, default = Path(tempfile.gettempdir()),
                         help = "Path to store temporary files (default: your OS's default temporary directory).")
    parser.add_argument('-gz', '--gzip', dest = 'gzip', default = False, action = 'store_true', 
                        help = "Gzip output (default: False).")
    taxon_names_fetch = parser.add_mutually_exclusive_group()
    taxon_names_fetch.add_argument('-tna', '--taxon-names-auto', dest = 'fetch_taxa_auto', default = False, action = 'store_true', 
                                   help = "Automatically adds taxon names to use as taxon labels instead of filenames (default: False). For NCBI files, this will fetch them from NCBI. For Bakta files, this will generate a generic taxon name locally.")
    taxon_names_fetch.add_argument('-tnf', '--taxon-names-file', dest = 'fetch_taxa_file', default = None, type = Path, 
                                   help = "File to fetch taxon names from to use as taxon labels instead of filenames (default: None).")
    parser.add_argument('-f', '--force', dest = 'force', default = False, action = 'store_true', 
                        help = "Force overwriting output (default: false).")
    parser.add_argument('-vv', '--verbosity', dest = 'verbosity', default = 3, type = int, choices = [0,1,2,3,4], 
                              help = "Console verbosity level (default: 3 (info))")
    parser.add_argument('-np', '--no-progress', dest = 'no_progress', default = False, action = "store_true", 
                        help = "Don't show progress bar (default: False).")
    parser.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit")      
    
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
        # Input folder
        if not args.input.is_dir():
            raise ValueError('Input folder does not exist.')
            
        # Output folder
        if args.output.exists():
            if args.force:
                LOG.warning("Output already exists, but it will be overwritten.")
            else:
                raise FileExistsError("Output already exists! Rerun with -f to overwrite it.")
        else:
            args.output.parent.mkdir(parents = True, exist_ok = True)
            
        # Temp folder will always be unique
        args.temp.mkdir(parents = True, exist_ok = True)
        args.temp = Path(tempfile.mkdtemp(dir = args.temp))
    
        match args.mode:
            case 'ncbi-gbff' | 'bakta-gbff':
                genbank_files = any(chain(*[args.input.glob(f'*{ext}') 
                                            for ext in ALLOWED_EXTENSIONS]))
                if not(genbank_files):
                    raise ValueError("Input folder does not contain any Genbank file (allowed extensions: .gbk, .gbff, .gb; can be gzipped).")
            case 'ncbi-package':
                genbank_files = any(chain(*[args.input.glob(f'ncbi_dataset/data/*/genomic{ext}')
                                            for ext in ALLOWED_EXTENSIONS]))
                if not(genbank_files):
                    raise ValueError("NCBI package does not contain any Genbank file (mind the .gbff extension; can be gzipped).")
            case 'tsv':
                if not(any(args.input.glob('*.tsv'))):
                    raise ValueError("Input folder does not contain TSV files (mind the .tsv extension).")
                    
        if args.fetch_taxa_file:
            if not args.fetch_taxa_file.is_file():
                raise ValueError("Taxon name file does not exist.")
                    
    except (ValueError, FileExistsError) as err:
        LOG.critical(err)
        raise err
    
    # Convert validated arguments to dictionary
    parsed_args = vars(args)
    
    # Further specify the gzip option
    if args.gzip:
        parsed_args['gzip'] = "gzip"
    else:
        parsed_args['gzip'] = "uncompressed"
    
    return parsed_args


def read_genome(file: str | Path):
    """
    Open the appropriate file handle for a genome file.
    
    Automatically distinguishes between compressed and uncompressed files based on the file extension.
    
    Args:
        file (str | Path): genome file to open
        
    Returns:
        handle: A file handle to open the genome file
    """
    if '.gz' in Path(file).suffixes:
        handle = gzip.open(file, mode = 'rt')
    else:
        handle = open(file, mode = 'r')
        
    return handle


def _add_one_ncbi_genbank_to_db(numbered_filepath: tuple, writer: DictWriter, in_package: bool = False) -> None:
    """
    Parses all CDS records in a single NCBI Genbank file and adds them to a temporary TSV file.
    
    Extracts CDS features and their associated metadata from an NCBI-formatted Genbank file,
    including taxonomic ID, gene tags, product names, genomic coordinates, and strand
    information. Joins exon coordinates for multi-exon CDS records.
    
    Args:
        numbered_filepath (tuple): A tuple containing (index, Path) where index is the file
            number and Path is the file path to the Genbank file.
        writer (csv.DictWriter): Preinitialised DictWriter object to write to the temp TSV file.
        in_package (bool): If True, extracts filename from parent directory (for NCBI package
            structure). If False, uses file stem as filename. Defaults to False.
    
    Returns:
        None
    """
    file = numbered_filepath[1]
    cds_records = []
    
    with read_genome(file) as handle:
        # First parse the Genbank
        records = list(SeqIO.parse(handle, 'genbank'))
        
        # Get filename
        if in_package:
            filelabel = file.parent.name
        else:
            filelabel = file.with_suffix('.gz').with_suffix('').with_suffix('').name
        
        # Get taxon id
        first_source_feature = [feat for feat in records[0].features if feat.type == 'source'][0]
        first_source_feature_dbxrefs = first_source_feature.qualifiers['db_xref']
        taxon_id = [xref for xref in first_source_feature_dbxrefs if 'taxon:' in xref][0]
        taxon_id = int(taxon_id.split(':')[1])
        
        # Parse all CDS features
        for record in records:
            cds_features = [feat for feat in record.features if feat.type == 'CDS']
            for feature in cds_features:
                
                # Try to fetch all necessary data for a CDS record
                try:
                    coord_intervals = ['..'.join([str(part.start), str(part.end)]) 
                                       for part in feature.location.parts]
                    cds_record = {
                        'gene_tag': feature.qualifiers['protein_id'][0],
                        'name': feature.qualifiers['product'][0],
                        'contig': record.id,
                        'coords': ','.join(coord_intervals),
                        'strand': '{0:+}'.format(feature.location.strand)[0],
                        'taxon_id': taxon_id,
                        'filelabel': filelabel,
                        }
                    
                    cds_records.append(cds_record)
                    
                # If some data is lacking, ignore this entry
                except KeyError:
                    continue
    
    # Write away the CDS records for this file
    writer.writerows(cds_records)
    
    return None


def _add_one_bakta_genbank_to_db(numbered_filepath: tuple, writer: DictWriter) -> None:
    """
    Parses all CDS records in a single Bakta Genbank file and adds them to a temporary TSV file.
    
    Extracts CDS features and their associated metadata from a Bakta-formatted Genbank file,
    including taxonomic ID, gene tags, product names, genomic coordinates, and strand
    information. Joins exon coordinates for multi-exon CDS records.
    
    Args:
        numbered_filepath (tuple): A tuple containing (index, Path) where index is the file
            number and Path is the file path to the Genbank file.
        writer (csv.DictWriter): Preinitialised DictWriter object to write to the temp TSV file.
    
    Returns:
        None
    """
    file = numbered_filepath[1]
    cds_records = []
    
    with read_genome(file) as handle:
        # First parse the Genbank
        records = list(SeqIO.parse(handle, 'genbank'))
        
        # Get filelabel
        filelabel = file.with_suffix('.gz').with_suffix('').with_suffix('').name
        
        # Assign generic taxon id
        taxon_id = numbered_filepath[0]
        
        # Parse all CDS features
        for record in records:
            cds_features = [feat for feat in record.features if feat.type == 'CDS']
            for feature in cds_features:
                
                # Try to fetch all necessary data for a CDS record
                try:
                    coord_intervals = ['..'.join([str(part.start), str(part.end)]) 
                                       for part in feature.location.parts]
                    cds_record = {
                        'gene_tag': feature.qualifiers['locus_tag'][0],
                        'name': feature.qualifiers['product'][0],
                        'contig': record.id,
                        'coords': ','.join(coord_intervals),
                        'strand': '{0:+}'.format(feature.location.strand)[0],
                        'taxon_id': taxon_id,
                        'filelabel': filelabel,
                        }
                    
                    cds_records.append(cds_record)
                    
                # If some data is lacking, ignore this entry
                except KeyError:
                    continue
    
    # Write away the CDS records for this file
    writer.writerows(cds_records)
    
    return None


def _add_one_tsv_to_db(numbered_filepath: tuple, writer: DictWriter) -> None:
    """
    Reads all CDS records from a custom TSV file and adds them to a temporary TSV file.
    
    Reads a TSV file with CDS coordinate data. Expected header:
    gene_tag, name, contig, coords, strand, taxon_id, filelabel.
    
    Args:
        numbered_filepath (tuple): A tuple containing (index, Path) where index is the file
            number and Path is the file path to the TSV file.
        writer (csv.DictWriter): Preinitialised DictWriter object to write to the temp TSV file.
    
    Returns:
        None
        
    Note:
        Genomic coordinates (coords) are expected to already have been formatted as joined exons. Example: "1..100,105..200".
    """
    # Envisioned header: ['gene_tag', 'name', 'contig', 'coords', 'strand', 'taxon_id', 'filelabel']
    df = pl.scan_csv(numbered_filepath[1], separator = "\t", has_header = True, comment_prefix = '#')
    
    # Check out the CDS records
    cds_records = df.collect().to_dicts()
    
    # Write away the CDS records for this file
    writer.writerows(cds_records)
    
    return None
    

def parse_files(input_path: Path, parsing_mode: str, temp_cds_db_path: Path,
                no_progress: bool = False) -> pl.LazyFrame:
    """
    Parses all input files and constructs a temporary CDS coordinates database.
    
    Dispatches file parsing based on the specified parsing mode, and returns
    a lazy entry point to the temporary TSV file for further processing.
    
    Args:
        input_path (Path): Path to the folder containing input files.
        parsing_mode (str): File format mode - one of: 'ncbi-gbff', 'ncbi-package',
            'bakta-gbff', or 'tsv'.
        no_progress (bool): If True, suppresses the progress bar during parsing.
            Defaults to False.
    
    Returns:
        A Polars LazyFrame containing concatenated CDS records from all input files
        with columns: gene_tag, name, contig, coords, strand, taxon_id, and filename
        (or taxon_name for TSV formats).
    """
    # Select the right workflows
    match parsing_mode:
        case 'ncbi-gbff':
            LOG.info('Parsing all input files as NCBI Genbanks')
            parser = _add_one_ncbi_genbank_to_db
            files = chain(*[input_path.glob(f'*{ext}') for ext in ALLOWED_EXTENSIONS])
        case 'ncbi-package':
            LOG.info('Parsing all input files as an NCBI Genbank package')
            parser = partial(_add_one_ncbi_genbank_to_db, in_package = True)
            files = chain(*[input_path.glob(f'ncbi_dataset/data/*/genomic{ext}') for ext in ALLOWED_EXTENSIONS])
        case 'bakta-gbff':
            LOG.info('Parsing all input files as Bakta Genbanks')
            parser = _add_one_bakta_genbank_to_db
            files = chain(*[input_path.glob(f'*{ext}') for ext in ALLOWED_EXTENSIONS])
        case 'tsv':
            LOG.info('Parsing all input files as TSVs')
            parser = _add_one_tsv_to_db
            files = input_path.glob('*.tsv')
    inputs = enumerate(files)
            
    # Parse all Genbanks and write them to a temporary TSV.GZ file
    fieldnames = ['gene_tag', 'name', 'contig', 'coords', 'strand', 'taxon_id', 'filelabel']
    with gzip.open(temp_cds_db_path, 'wt') as out_handle:
        writer = DictWriter(out_handle, fieldnames, delimiter = '\t')
        
        for numbered_filepath in tqdm(list(inputs), leave = False, disable = no_progress): 
            parser(numbered_filepath = numbered_filepath, writer = writer)
       
    # Return a LazyFrame entrypoint
    cds_db = pl.scan_csv(temp_cds_db_path, separator = "\t", has_header = False, new_columns = fieldnames)
                
    return cds_db


def check_duplicate_contigs(cds_db: pl.LazyFrame, parsing_mode: str) -> pl.LazyFrame:
    """
    Check for and attempt to fix duplicate contig labels per taxon.
    
    Detects cases where the same contig label appears in multiple taxa. For
    Bakta files, attempts to prepend the existing locus tag prefix to make contigs
    unique. For other formats, exits with an error.
    
    Args:
        cds_db (polars.LazyFrame): A dataframe containing CDS records with 'contig',
            'taxon_id', and 'gene_tag' columns.
        parsing_mode (str)): The format mode used for parsing.
    
    Returns:
        cds_db (polars.LazyFrame): The input DataFrame with modified contig labels if a fix
            was applied (Bakta mode only).
            
    Mutates: 
        cds_db (polars.LazyFrame): The input DataFrame with modified contig labels if a fix
            was applied (Bakta mode only).
                
    Raises:
        RuntimeError: If duplicate contigs are detected and cannot be fixed, or if
            fix attempt fails.
        
    Note:
        This function contains a potential local partial materialisation of the LazyFrame.
        This triggers all files to be parsed.
    """
    LOG.info('Checking for duplicate contig labels')
    
    # Check that no contig label occurs in combination with multiple taxon IDs
    contig_taxa_combs = cds_db.select(['contig', 'taxon_id']).group_by(['contig']).agg(pl.col('taxon_id').unique())
    contig_taxa_combs = contig_taxa_combs.collect() # Local materialisation
    nb_contig_taxa_combs = contig_taxa_combs.with_columns(nb_combs = pl.col('taxon_id').list.len())
    multi_taxa_contig = nb_contig_taxa_combs['nb_combs'] > 1
    if multi_taxa_contig.any():
        LOG.error("Detected duplicate contig labels for a taxon!")
        
        # In case of Bakta files, we may still fix this if Bakta autogenerated a unique locus tax prefix,
        # by prepending that locus tag prefix
        if parsing_mode == 'bakta-gbff':
            LOG.error("Trying to fix this by prepending the locus tag prefix...")
            # Extract the gene tag prefix
            cds_db = cds_db.with_columns(locus_tag_prefix = pl.col('gene_tag').str.split('_').list.reverse().list.slice(1).list.reverse().list.join('_'))
            # Preprend it
            cds_db = cds_db.with_columns(contig = pl.concat_str(['contig', 'locus_tag_prefix'], separator = '_'))
            cds_db = cds_db.drop('locus_tag_prefix')
            
            # Check if fix succeeded
            # => No contig label occurring in combination with multiple taxon IDs
            contig_taxa_combs_fixed = cds_db.group_by(['contig']).agg(pl.col('taxon_id').unique())
            nb_contig_taxa_combs_fixed = contig_taxa_combs_fixed.with_columns(nb_combs = pl.col('taxon_id').list.len())
            
            multi_taxa_contig_fixed = nb_contig_taxa_combs_fixed['nb_combs'] > 1
            contigs_fixed_counts = not(cds_db.select('contig').is_duplicated().any())
            
            # If Not all contig labels are unique despite the fix
            if multi_taxa_contig_fixed.any() and contigs_fixed_counts:
                msg = 'Duplicate gene tag fix failed! Exiting.'
                LOG.critical(msg)
                raise RuntimeError(msg)
            else:
                LOG.warning('Fix succeeded, but I had to change your contig labels!')
        
        # For other input types, we don't provide a potential fix
        else:
            msg = "Duplicate contig labels detected for a taxon! No fix provided for this parsing mode."
            LOG.critical(msg)
            raise RuntimeError(msg)
            
    return cds_db


def set_taxon_labels(cds_db: pl.LazyFrame, fetch_taxa_auto: bool, fetch_taxa_file: Path, 
                     parsing_mode: str, batch_size: int = 250, max_attempts: int = 5) -> pl.LazyFrame:
    """
    Set taxon labels as either scientific names or filenames.
    
    For NCBI files, optionally fetches the scientific names from NCBI Taxonomy
    via BioPython's NCBI Entrez API with retry logic.
    For Bakta files, generates generic labels or uses filenames.
    For TSV files, preserves user-provided annotations.
    
    Args:
        cds_db (polars LazyFrame): Dataframe containing CDS records with 'taxon_id',
            'filename', and optionally 'gene_tag' columns.
        fetch_taxa_ncbi (bool): If True, uses scientific names (NCBI) or generates generic
            names (Bakta) to use as taxon names instead of the filenames. If false,
            the default filenames will be kept, unless a rename file was supplied.
        fetch_taxa_file (Path | None): Path to the rename file with the taxon names to
            replace the current ones sourced from the filenames. Defaults to None.
        parsing_mode (str)): The format mode used for parsing.
        batch_size (int): Number of taxon names to fetch in one batch. Defaults to 250.
        max_attempts (int): Maximum numbers of times to attempt fetching the taxon names
            using Entrez. Defaults to 5.
        
    Returns:
        cds_db (polars.LazyFrame): The input DataFrame with a new 'taxon_name' column and
            the 'filename' column removed.
            
    Mutates:
        cds_db (polars.LazyFrame): The input DataFrame with a new 'taxon_name' column and
            the 'filename' column removed.
            
    Note:
        This function contains a potential local partial materialisation of the LazyFrame
        when fetching taxon names automatically. This triggers all files to be parsed an additional time.
    """
    match parsing_mode:
        # In case of NCBI files
        case str(parsing_mode) if 'ncbi' in parsing_mode:
            
            # Fetch taxon names on-the-fly from NCBI if requested
            if fetch_taxa_auto:
                LOG.info('Fetching taxon names using NCBI Entrez')
                all_taxon_ids = cds_db.select('taxon_id').unique()
                all_taxon_ids = all_taxon_ids.collect().to_series().to_list() # Local materialisation
                
                # Split the list up in batches, and fetch each batch
                with logging_redirect_tqdm(loggers = [LOG]):
                    all_taxon_names = []
                    for batch_idx, batch_ids in tqdm(list(enumerate(batched(all_taxon_ids, batch_size))),
                                                     leave = False):
                        try:
                            batch_names = fetch_taxon_names(batch_ids, max_attempts = max_attempts)
                            all_taxon_names.append(batch_names)
                        except RuntimeError as err:
                            LOG.error(f"Error fetching taxon names for batch {batch_idx}!")
                            raise err
                
                # Chain all fetched lists of taxon names
                all_taxon_names = list(chain(*all_taxon_names))
                
                # Join with the CDS DB
                LOG.info('Adding taxon name column')
                id_name_map = pl.DataFrame({'taxon_id': all_taxon_ids, 'taxon_name': all_taxon_names})
                id_name_map = dict(zip(id_name_map['taxon_id'], id_name_map['taxon_name']))
                cds_db = cds_db.with_columns(pl.col('taxon_id').replace(id_name_map).alias('taxon_name'))
            
            # Use taxon names from local mapping file if requested
            elif fetch_taxa_file:
                LOG.info('Fetching taxon names from local file')
                id_name_map = pl.read_csv(fetch_taxa_file, separator = "\t", new_columns = ['old_name', 'new_name'])
                id_name_map = dict(zip(id_name_map['old_name'], id_name_map['new_name']))
                
                # Replace taxon labels
                LOG.info('Adding taxon name column')
                cds_db = cds_db.with_columns(pl.col('filelabel').replace(id_name_map).alias('taxon_name'))
        
            # Default behaviour: use filenames
            else:
                LOG.info('Using filenames as taxon names')
                cds_db = cds_db.with_columns(pl.col('filelabel').alias('taxon_name'))
            
        # In case of Bakta files
        case 'bakta-gbff':
            # Generate a generic taxon name if requested
            if fetch_taxa_auto:
                LOG.info('Generating generic taxon names')
                cds_db = cds_db.with_columns("Taxon " + pl.col('taxon_id').alias('taxon_name'))
                
            # Use taxon names from local mapping file if requested
            elif fetch_taxa_file:
                LOG.info('Fetching taxon names from local file')
                replacements = pl.read_csv(fetch_taxa_file, separator = "\t", new_columns = ['old_name', 'new_name'])
                replacements = dict(zip(replacements['old_name'], replacements['new_name']))
                
                # Replace taxon labels
                LOG.info('Adding taxon name column')
                cds_db = cds_db.with_columns(pl.col('filelabel').replace(replacements).alias('taxon_name'))
                
            # Default behaviour: use filenames
            else:
                LOG.info('Using filenames as taxon names')
                cds_db = cds_db.with_columns(pl.col('filelabel').alias('taxon_name'))
               
        # For other parsing modes, keep the user's annotations
        case _:
            pass
            
    return cds_db


def fetch_taxon_names(taxon_ids: list, max_attempts: int = 5) -> list:
    """
    Fetch NCBI taxon names for a given list of NCBI taxon IDs.
    
    Fetches summary objects for a given list of NCBI Taxonomy IDs using BioPython's
    Entrez API with retry logic, and extracts the taxon name from each summary.
    
    Args:
        taxon_ids (list): List of NCBI taxon IDs as strings.
        max_attempts (int): Maximum numbers of attempts to retry fetching in
            case of a failure. Defaults to 5.
            
    Returns:
        taxon_names (list): List of NCBI taxon names
    """
    for attempt in range(max_attempts):
        try:
            with Entrez.esummary(db = 'taxonomy', id = taxon_ids) as handle:
                records = list(Entrez.read(handle))
            taxon_names = [str(i['ScientificName']) for i in records]
            break
        except:
            if attempt+1<max_attempts:
                LOG.warning(f'Failed fetching taxon names in attempt {attempt+1}. Retrying...')
                continue
            else:
                msg = f'Failed fetching taxon names in {max_attempts} attempts. Giving up...'
                LOG.error(msg)
                raise RuntimeError(msg)
    
    return taxon_names


def run_workflow(parsed_args: dict) -> None:
    """
    Execute the complete CDS database construction workflow.
    
    Loads and parses input files, validates contig uniqueness, assigns taxon labels,
    and writes the final CDS coordinates database to disk as a tab-separated file.
    Supports optional gzip compression.
    
    Args:
        parsed_args (dict): A dictionary holding the parsed and validated argument values.
    
    Returns:
        None
    
    Note:
        This workflow uses a lazy parsing method to limit RAM usage. This comes at
        the cost of the files being parsed multiple times, which is slower.
    """
    # Unpack arguments
    (input_path, output_path, 
     temp_path, parsing_mode,
     fetch_taxa_auto, fetch_taxa_file, 
     no_progress, gzip) = map(parsed_args.get, ['input', 'output', 'temp', 'mode', 
                                                'fetch_taxa_auto', 'fetch_taxa_file', 'no_progress', 'gzip'])
    temp_cds_db_path = temp_path / 'temp_cds_db.tsv.gz'
        
    # Parse the input files
    cds_db = parse_files(input_path, parsing_mode, temp_cds_db_path = temp_cds_db_path, no_progress = no_progress)
    
    # Check for duplicate contig labels
    cds_db = check_duplicate_contigs(cds_db, parsing_mode)
    
    # Keep taxon names as assembly labels or use the filenames
    cds_db = set_taxon_labels(cds_db, fetch_taxa_auto, fetch_taxa_file, parsing_mode)
        
    # Write results
    LOG.info('Writing DB to disk')
    cds_db = cds_db.select(['gene_tag', 'name', 'contig', 'strand', 'coords', 'taxon_id', 'taxon_name', 'filelabel'])
    cds_db.sink_csv(output_path, separator = '\t', include_header = False, compression = gzip)
    
    return None


def main():
    """
    Main entry point for the CDS database construction tool.
    
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
    
    
if __name__ == "__main__":
    main()

