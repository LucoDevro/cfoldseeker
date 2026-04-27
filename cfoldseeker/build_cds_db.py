#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import logging
import warnings
warnings.filterwarnings('ignore')
from Bio import Entrez
import polars as pl
from pathlib import Path
from tqdm.contrib.concurrent import thread_map


LOG = logging.getLogger(__name__)
logging.basicConfig(
    level = logging.INFO,
    format = "[%(asctime)s] %(levelname)s [%(filename)s: %(funcName)s] - %(message)s",
    datefmt="%H:%M:%S",
    handlers = [logging.StreamHandler(sys.stdout)]
    )


def parse_arguments() -> argparse.Namespace:
    """
    This function parses the arguments given through the command line.
    
    Args:
        None
    
    Returns:
        A Namespace object holding the parsed arguments
    """
    
    parser = argparse.ArgumentParser(
        prog = 'cfoldseeker-cds',
                epilog = 
                """
                Lucas De Vrieze, Miguel Biltjes
                (c) 2026 Masschelein lab, VIB
                """,
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = 
                """
                Helper tool to construct a CDS coordinates DB for cfoldseeker
                """,
                add_help = False
                )
    
    parser.add_argument('-i', '--input', dest = 'input', type = Path, default = Path('.'), help = "Path to folder holding the input files or NCBI package (default: current directory)")
    parser.add_argument('-m', '--mode', dest = 'mode', type = str, required = True, choices = ['ncbi-gff', 'ncbi-package', 'bakta-gff', 'tsv'],
                        help = 'File parsing mode (choices: ncbi-gff, ncbi-package, bakta-gff, tsv).')
    parser.add_argument('-o', '--output', dest = 'output', type = Path, default = Path('local_db'), help = "Filepath to save CDS coordinate DB (default: local_db).")
    parser.add_argument('-gz', '--gzip', dest = 'gzip', default = False, action = 'store_true', help = "Gzip output (default: False).")
    parser.add_argument('-tn', '--use-taxon-names', dest = 'use_taxa', default = False, action = 'store_true', help = "Use taxon names as labels for all files instead of filenames (default: False).")
    parser.add_argument('-c', '--cores', dest = 'cores', type = int, default = 1, help = "Number of cores available to use (default: 1).")
    parser.add_argument('-f', '--force', dest = 'force', default = False, action = 'store_true', help = "Force overwriting output (default: false).")
    parser.add_argument('-np', '--no-progress', dest = 'no_progress', default = False, action = "store_true", help = "Don't show progress bar (default: False).")
    parser.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit")      

    args = parser.parse_args()
    
    if not args.input.is_dir():
        msg = 'Input folder does not exist.'
        LOG.critical(msg)
        raise argparse.ArgumentError(msg)
    
    match args.mode:
        case 'ncbi-gff' | 'bakta-gff':
            any(args.input.glob('*.gff')), "Input folder does not contain GFF files (mind the .gff extension)."
        case 'ncbi-package':
            any(args.input.glob('ncbi_dataset/data/*/genomic.gff')), "NCBI package does not contain GFF files."
        case 'tsv':
            any(args.input.glob('*.tsv')), "Input folder does not contain TSV files (mind the .tsv extension)."
    
    if args.output.exists():
        if args.force:
            LOG.warning("Output already exists, but it will be overwritten.")
        else:
            msg = "Output already exists! Rerun with -f to overwrite it."
            LOG.error(msg)
            raise IOError(msg)
    else:
        args.output.parent.mkdir(parents = True, exist_ok = True)
    
    return args


def _parse_one_ncbi_gff(numbered_filepath: tuple, in_package: bool = False) -> pl.DataFrame:
    """
    Parses a single NCBI GFF file into a Polars DataFrame.
    
    Extracts CDS features and their associated metadata from an NCBI-formatted GFF file,
    including taxonomic ID, gene tags, product names, genomic coordinates, and strand
    information. Aggregates exon coordinates for multi-exon CDS records.
    
    Args:
        numbered_filepath (tuple): A tuple containing (index, Path) where index is the file
            number and Path is the file path to the GFF file.
        in_package (bool): If True, extracts filename from parent directory (for NCBI package
            structure). If False, uses file stem as filename. Defaults to False.
    
    Returns:
        A Polars DataFrame with columns: gene_tag, name, contig, coords, strand,
        taxon_id, and filename. Exons are aggregated into comma-separated coordinates.
        
    Note:
        The column filename is a temporary column that is removed by a later method in the workflow.
        It is necessary to construct a taxon name column with NCBI Assembly accession IDs.
    """
    df = pl.scan_csv(numbered_filepath[1], separator = "\t", has_header = False, comment_prefix = '#',
                     new_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
                     ).select(['seqid', 'type', 'start', 'end', 'strand', 'attributes'])
    
    df_regions = df.filter(pl.col('type') == 'region')
    df_cds = df.filter(pl.col('type') == 'CDS')
    
    # Extract taxon ID
    taxon_id = df_regions.select(pl.col('attributes').str.extract(r'taxon:([0-9]+)')).head(1)
    taxon_id = taxon_id.collect().rows()[0][0]
    
    # Get filename
    if in_package:
        filename = numbered_filepath[1].parent.name
    else:
        filename = numbered_filepath[1].stem

    # Populate the CDS coords DB
    df_cds = df_cds.drop('type')
    cds_record = df_cds.select([pl.col('attributes').str.extract(r'protein_id=([^;]+)').alias('gene_tag'),
                            pl.col('attributes').str.extract(r'product=([^;]+)').alias('name'),
                            pl.col('seqid').alias('contig'),
                            pl.concat_str(['start', 'end'], separator = '..').alias('coords'),
                            pl.col('strand'),
                            pl.lit(taxon_id).alias('taxon_id'),
                            pl.lit(filename).alias('filename') # temporary filename column to help construct a taxon name column later
                            ])
    cds_record = cds_record.drop_nulls(subset = 'gene_tag')
    
    # Aggregate multiple CDSes (i.e. exons) in one record
    cds_record = cds_record.group_by(pl.all().exclude('coords')).agg(pl.col('coords').str.join(','))
    
    return cds_record.collect()


def _parse_one_bakta_gff(numbered_filepath: tuple) -> pl.DataFrame:
    """
    Parses a single Bakta GFF file into a Polars DataFrame.
    
    Extracts CDS features and metadata from a Bakta-formatted GFF file, including
    locus tags, product names, genomic coordinates, and strand information. Aggregates
    exon coordinates for multi-exon CDS records.
    
    Args:
        numbered_filepath (tuple): A tuple containing (index, Path) where index is used as
            the generic taxon ID and Path is the file path to the GFF file.
    
    Returns:
        A Polars DataFrame with columns: gene_tag, name, contig, coords, strand,
        taxon_id, and filename. Exons are aggregated into comma-separated coordinates.
    """
    df = pl.scan_csv(numbered_filepath[1], separator = "\t", has_header = False, comment_prefix = '#',
                     new_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
                     ).select(['seqid', 'type', 'start', 'end', 'strand', 'attributes'])
    
    df_cds = df.filter(pl.col('type') == 'CDS')
    
    # Assign generic taxon ID
    taxon_id = numbered_filepath[0]
    
    # Get filename
    filename = numbered_filepath[1].stem
    
    # Populate the CDS coords DB
    df_cds = df_cds.drop('type')
    cds_record = df_cds.select([pl.col('attributes').str.extract(r'locus_tag=([^;]+)').alias('gene_tag'),
                            pl.col('attributes').str.extract(r'product=([^;]+)').alias('name'),
                            pl.col('seqid').alias('contig'),
                            pl.concat_str(['start', 'end'], separator = '..').alias('coords'),
                            pl.col('strand'),
                            pl.lit(taxon_id).alias('taxon_id'),
                            pl.lit(filename).alias('filename') # temporary filename column to help construct a taxon name column later
                            ])
    cds_record = cds_record.drop_nulls(subset = 'gene_tag')
    
    # Aggregate multiple CDSes (i.e. exons) in one record
    cds_record = cds_record.group_by(pl.all().exclude('coords')).agg(pl.col('coords').str.join(','))
    
    return cds_record.collect()


def _parse_one_tsv(numbered_filepath: tuple) -> pl.DataFrame:
    """
    Parses a single TSV file into a Polars DataFrame.
    
    Reads a tab-separated values file with CDS coordinate data. Expected header:
    gene_tag, name, contig, start, end, strand, taxon_id, taxon_name. Aggregates
    exon coordinates for multi-exon CDS records.
    
    Args:
        numbered_filepath (tuple): A tuple containing (index, Path) where index is unused
            and Path is the file path to the TSV file.
    
    Returns:
        A Polars DataFrame with columns: gene_tag, name, contig, coords, strand,
        taxon_id, taxon_name. Exons are aggregated into comma-separated coordinates.
    """
    # Envisioned header: ['gene_tag', 'name', 'contig', 'start', 'end', 'strand', 'taxon_id', 'taxon_name']
    # Exons at multiple lines
    df = pl.scan_csv(numbered_filepath[1], separator = "\t", has_header = True, comment_prefix = '#')
    
    # Populate the CDS coords DB
    cds_record = df.with_columns(pl.concat_str(['start', 'end'], separator = '..').alias('coords'))
    cds_record = cds_record.drop(['start', 'end'])
    
    # Aggregate multiple CDSes (i.e. exons) in one record
    cds_record = cds_record.group_by(pl.all().exclude('coords')).agg(pl.col('coords').str.join(','))
    
    return cds_record.collect()


def parse_inputs(input_path: Path, parsing_mode: str, n_workers: int = 1, no_progress: bool = False) -> pl.DataFrame:
    """
    Parses all input files and constructs a draft CDS coordinates database.
    
    Dispatches file parsing based on the specified parsing mode, using
    parallel processing to handle multiple files efficiently. Concatenates all parsed
    results into a single DataFrame.
    
    Args:
        input_path (Path): Path to the folder containing input files.
        parsing_mode (str): File format mode - one of: 'ncbi-gff', 'ncbi-package',
            'bakta-gff', or 'tsv'.
        n_workers (int): Number of worker threads for parallel file parsing.
            Defaults to 1.
        no_progress (bool): If True, suppresses the progress bar during parsing.
            Defaults to False.
    
    Returns:
        A Polars DataFrame containing concatenated CDS records from all input files
        with columns: gene_tag, name, contig, coords, strand, taxon_id, and filename
        (or taxon_name for TSV formats).
    """
    # Parse all GFFs
    match parsing_mode:
        case 'ncbi-gff':
            LOG.info('Parsing all input files as NCBI GFFs')
            parsed_gffs_to_concat = thread_map(_parse_one_ncbi_gff, list(enumerate(input_path.glob('*.gff'))),
                                               max_workers = n_workers,
                                               leave = False,
                                               disable = no_progress)
        case 'ncbi-package':
            LOG.info('Parsing all input files as an NCBI GFF package')
            parsed_gffs_to_concat = thread_map(lambda x: _parse_one_ncbi_gff(x, in_package = True), 
                                               list(enumerate(input_path.glob('ncbi_dataset/data/*/genomic.gff'))),
                                               max_workers = n_workers,
                                               leave = False,
                                               disable = no_progress)
        case 'bakta-gff':
            LOG.info('Parsing all input files as Bakta GFFs')
            parsed_gffs_to_concat = thread_map(_parse_one_bakta_gff, list(enumerate(input_path.glob('*.gff'))),
                                               max_workers = n_workers,
                                               leave = False,
                                               disable = no_progress)
        case 'tsv':
            LOG.info('Parsing all input files as TSVs')
            parsed_gffs_to_concat = thread_map(_parse_one_tsv, list(enumerate(input_path.glob('*.tsv'))),
                                               max_workers = n_workers,
                                               leave = False,
                                               disable = no_progress)
                
    # Concatenate all parsed GFF tables
    LOG.info('Constructing CDS coordinates DB')
    cds_db = pl.concat(parsed_gffs_to_concat)
    
    return cds_db


def check_duplicate_contigs(cds_db: pl.DataFrame, parsing_mode: str) -> pl.DataFrame:
    """
    Check for and attempt to fix duplicate contig labels per taxon.
    
    Detects cases where the same contig label appears in multiple taxa. For
    Bakta GFF files, attempts to prepend the existing locus tag prefix to make contigs
    unique. For other formats, exits with an error.
    
    Args:
        cds_db (polars.DataFrame): A dataframe containing CDS records with 'contig',
            'taxon_id', and 'gene_tag' columns.
        parsing_mode (str)): The format mode used for parsing ('bakta-gff', 'ncbi-gff',
            'ncbi-package', or 'tsv').
    
    Returns:
        cds_db (polars.DataFrame): The input DataFrame with modified contig labels if a fix
            was applied (Bakta mode only).
            
    Mutates: 
        cds_db (polars.DataFrame): The input DataFrame with modified contig labels if a fix
            was applied (Bakta mode only).
                
    Raises:
        RuntimeError: If duplicate contigs are detected and cannot be fixed, or if
            fix attempt fails.    
    """
    # Check that no contig label occurs in combination with multiple taxon IDs
    contig_taxa_combs = cds_db.group_by(['contig']).agg(pl.col('taxon_id').unique())
    nb_contig_taxa_combs = contig_taxa_combs.with_columns(nb_combs = pl.col('taxon_id').list.len())
    multi_taxa_contig = nb_contig_taxa_combs['nb_combs'] > 1
    if multi_taxa_contig.any():
        LOG.error("Detected duplicate contig labels for a taxon!")
        
        # In case of Bakta GFFs, we may still fix this if Bakta autogenerated a unique locus tax prefix,
        # by prepending that locus tag prefix
        if parsing_mode == 'bakta-gff':
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


def set_taxon_labels(cds_db: pl.DataFrame, use_taxa: bool, parsing_mode: str, max_attempts: int = 3) -> pl.DataFrame:
    """
    Set taxon labels as either scientific names or filenames.
    
    For NCBI files, optionally fetches the scientific names from NCBI Taxonomy
    via BioPython's NCBI Entrez API with retry logic.
    For Bakta GFF files, generates generic labels or uses filenames.
    For TSV files, preserves user-provided annotations.
    
    Args:
        cds_db (polars DataFrame): Dataframe containing CDS records with 'taxon_id',
            'filename', and optionally 'gene_tag' columns.
        use_taxa (bool): If True, uses scientific names (NCBI) or generates generic
            names (Bakta). If False, uses filenames as taxon labels.
        parsing_mode (str)): The format mode used for parsing ('ncbi-gff', 'ncbi-package',
            'bakta-gff', or 'tsv').
        max_attempts (int): Maximum numbers of times to attempt fetching the taxon names
            using Entrez. Defaults to 3.
        
    Returns:
        cds_db (polars.DataFrame): The input DataFrame with a new 'taxon_name' column and
            the 'filename' column removed.
            
    Mutates:
        cds_db (polars.DataFrame): The input DataFrame with a new 'taxon_name' column and
            the 'filename' column removed.
            
    Note:
        This function removes the temporary column 'filename' if present, as it may
        have been introduced when parsing NCBI GFF files.
    """
    # In case of NCBI files
    if 'ncbi' in parsing_mode:
        # Fetch all taxon names if requested
        if use_taxa:
            LOG.info('Fetching taxon names using NCBI Entrez')
            all_taxon_ids = cds_db.select('taxon_id').unique().to_series().to_list()
            for attempt in range(max_attempts):
                try:
                    with Entrez.esummary(db = 'taxonomy', id = all_taxon_ids) as handle:
                        records = list(Entrez.read(handle))
                    all_taxon_names = [str(i['ScientificName']) for i in records]
                    break
                except:
                    if attempt+1<max_attempts:
                        LOG.warning(f'Failed fetching taxon names in attempt {attempt+1}. Retrying...')
                        continue
                    else:
                        msg = f'Failed fetching taxon names in {max_attempts} attempts. Giving up...'
                        LOG.error(msg)
                        raise RuntimeError(msg)
            
            # Join with the CDS DB
            LOG.info('Adding taxon name column')
            id_name_map = pl.DataFrame({'taxon_id': all_taxon_ids, 'taxon_name': all_taxon_names})
            cds_db = cds_db.join(id_name_map, on = 'taxon_id', how = 'left', maintain_order = "left")
        
        # If not requested, use filenames
        else:
            cds_db = cds_db.with_columns(pl.col('filename').alias('taxon_name'))
            
    # In case of Bakta GFF files
    elif parsing_mode == 'bakta-gff':
        # Generate a generic taxon name if requested
        if use_taxa:
            LOG.info('Generating generic taxon names')
            cds_db = cds_db.with_columns("Taxon " + pl.col('taxon_id').alias('taxon_name'))
        # Otherwise, use the filename
        else:
            cds_db = cds_db.with_columns(pl.col('filename').alias('taxon_name'))
            
    # For other parsing modes, keep the user's annotations
            
    # Drop the temporary filename column
    cds_db = cds_db.drop('filename', strict = False)
    
    return cds_db


def main():
    """
    Main entry point for the CDS database construction tool.
    
    Orchestrates the complete workflow: parses command-line arguments, loads
    and parses input files, validates contig uniqueness, assigns taxon labels,
    and writes the final CDS coordinates database to disk as a tab-separated file.
    Supports optional gzip compression.
    """
    # Process arguments
    args = parse_arguments()
    cores = args.cores
    input_path = args.input
    output_path = args.output
    parsing_mode = args.mode
    use_taxa = args.use_taxa
    no_progress = args.no_progress
    if args.gzip:
        gzip = "gzip"
    else:
        gzip = "uncompressed"
        
    # Parse the input files
    cds_db = parse_inputs(input_path, parsing_mode, n_workers = cores, no_progress = no_progress)
    
    # Check for duplicate contig labels
    cds_db = check_duplicate_contigs(cds_db, parsing_mode)
    
    # Keep taxon names as assembly labels or use the filenames
    cds_db = set_taxon_labels(cds_db, use_taxa, parsing_mode)
        
    # Write results
    LOG.info('Writing DB to disk')
    cds_db = cds_db.select(['gene_tag', 'name', 'contig', 'strand', 'coords', 'taxon_id', 'taxon_name'])
    cds_db.write_csv(output_path, separator = '\t', include_header = False, compression = gzip)
    
    
if __name__ == "__main__":
    main()

