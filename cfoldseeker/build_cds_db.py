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
from concurrent.futures import ThreadPoolExecutor

LOG = logging.getLogger(__name__)
logging.basicConfig(
    level = logging.INFO,
    format = "[%(asctime)s] %(levelname)s [%(filename)s: %(funcName)s] - %(message)s",
    datefmt="%H:%M:%S",
    handlers = [logging.StreamHandler(sys.stdout)]
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog = 'build_cds_db.py',
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
    
    parser.add_argument('-i', '--input', dest = 'input', type = Path, default = Path('.'), help = "Path to folder containing the input files (default: current directory)")
    parser.add_argument('-m', '--mode', dest = 'mode', type = str, required = True, choices = ['ncbi-gff', 'ncbi-package', 'bakta-gff', 'tsv', 'excel'],
                        help = 'File parsing mode (choices: ncbi-gff, bakta-gff, tsv, excel).')
    parser.add_argument('-o', '--output', dest = 'output', type = Path, default = Path('local_db'), help = "Filepath to save CDS coordinate DB (default: local_db).")
    parser.add_argument('-gz', '--gzip', dest = 'gzip', default = False, action = 'store_true', help = "Gzip output (default: False).")
    parser.add_argument('-tn', '--use-taxon-names', dest = 'use_taxa', default = False, action = 'store_true', help = "Use taxon names as labels for all files instead of filenames (default: False).")
    parser.add_argument('-c', '--cores', dest = 'cores', type = int, default = 1, help = "Number of cores available to use (default: 1).")
    parser.add_argument('-f', '--force', dest = 'force', default = False, action = 'store_true', help = "Force overwriting output (default: false).")
    parser.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit")      

    args = parser.parse_args()
    
    assert args.input.exists() and args.input.is_dir(), 'Input folder does not exist.'
    match args.mode:
        case 'ncbi-gff' | 'bakta-gff':
            any(args.input.glob('*.gff')), "Input folder does not contain GFF files (mind the .gff extension)."
        case 'ncbi-package':
            any(args.input.glob('ncbi_dataset/data/*/genomic.gff')), "NCBI package does not contain GFF files."
        case 'tsv':
            any(args.input.glob('*.tsv')), "Input folder does not contain TSV files (mind the .tsv extension)."
        case 'excel':
            any(args.input.glob('*.xlsx')), "Input folder does not contain Excel files (mind the .xlsx extension)."
    if args.output.exists():
        if args.force:
            LOG.warning("Output already exists, but it will be overwritten.")
        else:
            LOG.error("Output already exists! Rerun with -f to overwrite it.")
            sys.exit()
    else:
        args.output.parent.mkdir(parents = True, exist_ok = True)
    
    return args


def parse_one_ncbi_gff(numbered_filepath: tuple, in_package: bool = False) -> pl.DataFrame:
    """
    Parse one NCBI GFF file into a Polars DataFrame
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


def parse_one_bakta_gff(numbered_filepath: tuple) -> pl.DataFrame:
    """
    Parse one Bakta GFF file into a Polars DataFrame
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


def parse_one_tsv(numbered_filepath: tuple) -> pl.DataFrame:
    """
    Parse one TSV file into a Polars DataFrame
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


def parse_one_excel(numbered_filepath: tuple) -> pl.DataFrame:
    """
    Parse one Excel file into a Polars DataFrame
    """
    # Envisioned header: ['gene_tag', 'name', 'contig', 'start', 'end', 'strand', 'taxon_id', 'taxon_name']
    # Exons at multiple lines
    df = pl.read_excel(numbered_filepath[1], has_header = True).lazy()
    
    # Populate the CDS coords DB
    cds_record = df.with_columns(pl.concat_str(['start', 'end'], separator = '..').alias('coords'))
    cds_record = cds_record.drop(['start', 'end'])
    
    # Aggregate multiple CDSes (i.e. exons) in one record
    cds_record = cds_record.group_by(pl.all().exclude('coords')).agg(pl.col('coords').str.join(','))
    
    return cds_record.collect()


def parse_inputs(input_path: Path, parsing_mode: str, n_workers: int = 1) -> pl.DataFrame:
    """
    Parses all input files and returns a draft CDS coordinates DB.
    """
    # Parse all GFFs
    with ThreadPoolExecutor(max_workers = n_workers) as executor:
        match parsing_mode:
            case 'ncbi-gff':
                LOG.info('Parsing all input files as NCBI GFFs')
                parsed_gffs_to_concat = executor.map(parse_one_ncbi_gff, enumerate(input_path.glob('*.gff')))
            case 'ncbi-package':
                LOG.info('Parsing all input files as an NCBI GFF package')
                parsed_gffs_to_concat = executor.map(lambda x: parse_one_ncbi_gff(x, in_package = True), 
                                                     enumerate(input_path.glob('ncbi_dataset/data/*/genomic.gff')))
            case 'bakta-gff':
                LOG.info('Parsing all input files as Bakta GFFs')
                parsed_gffs_to_concat = executor.map(parse_one_bakta_gff, enumerate(input_path.glob('*.gff')))
            case 'tsv':
                LOG.info('Parsing all input files as TSVs')
                parsed_gffs_to_concat = executor.map(parse_one_tsv, enumerate(input_path.glob('*.tsv')))
            case 'excel':
                LOG.info('Parsing all input files as Excel files')
                parsed_gffs_to_concat = executor.map(parse_one_excel, enumerate(input_path.glob('*.xlsx')))
                
    # Concatenate all parsed GFF tables
    LOG.info('Constructing CDS coordinates DB')
    cds_db = pl.concat(parsed_gffs_to_concat)
    
    return cds_db


def check_duplicate_contigs(cds_db: pl.DataFrame, parsing_mode: str) -> pl.DataFrame:
    """
    Checks for duplicate contig labels per taxon. Tries to correct these in case of Bakta GFFs.
    """
    # Check that no contig label occurs in combination with multiple taxon IDs
    contig_taxa_combs = cds_db.group_by(['contig']).agg(pl.col('taxon_id').unique())
    nb_contig_taxa_combs = contig_taxa_combs.with_columns(nb_combs = pl.col('taxon_id').list.len())
    multi_taxa_contig = nb_contig_taxa_combs['nb_combs'] > 1
    if multi_taxa_contig.any():
        LOG.error("Detected duplicate contig labels for a taxon!")
        
        # In case of Bakta GFFs, we may still fix this if Bakta autogenerated a unique assembly ID.
        if parsing_mode == 'bakta-gff':
            LOG.error("Trying to fix this by prepending the locus tag prefix...")
            # Prepend gene tag prefix as a potential fix
            cds_db = cds_db.with_columns(locus_tag_prefix = pl.col('gene_tag').str.split('_').list.reverse().list.slice(1).list.reverse().list.join('_'))
            cds_db = cds_db.with_columns(contig = pl.concat_str(['contig', 'locus_tag_prefix'], separator = '_'))
            cds_db = cds_db.drop('locus_tag_prefix')
            
            # Check if fix succeeded
            # No contig label occurring in combination with multiple taxon IDs
            contig_taxa_combs_fixed = cds_db.group_by(['contig']).agg(pl.col('taxon_id').unique())
            nb_contig_taxa_combs_fixed = contig_taxa_combs_fixed.with_columns(nb_combs = pl.col('taxon_id').list.len())
            multi_taxa_contig_fixed = nb_contig_taxa_combs_fixed['nb_combs'] > 1
            # Not all contig labels are unqiue
            contigs_fixed_counts = not(cds_db.select('contig').is_duplicated().any())
            if multi_taxa_contig_fixed.any() and contigs_fixed_counts:
                LOG.critical('Fix failed! Exiting.')
                sys.exit()
            else:
                LOG.error('Fix succeeded, but I had to change your contig labels!')
        
        # For other input types, we don't provide a potential fix
        else:
            LOG.critical("Please make sure your contig labels are unique for each taxon! Exiting.")
            sys.exit()
            
    return cds_db


def set_taxon_labels(cds_db: pl.DataFrame, use_taxa: bool, parsing_mode: str) -> pl.DataFrame:
    """
    Sets the taxon labels. Taxon names are more human-friendly, while filenames are better for downstream data joining purposes
    """
    # In case of NCBI files
    if 'ncbi' in parsing_mode:
        if use_taxa:
            # Fetch all taxon names
            LOG.info('Fetching taxon names using NCBI Entrez')
            all_taxon_ids = cds_db.select('taxon_id').unique().to_series().to_list()
            with Entrez.esummary(db = 'taxonomy', id = all_taxon_ids) as handle:
                records = list(Entrez.read(handle))
            all_taxon_names = [str(i['ScientificName']) for i in records]
            
            # Join with the CDS DB
            LOG.info('Adding taxon name column')
            id_name_map = pl.DataFrame({'taxon_id': all_taxon_ids, 'taxon_name': all_taxon_names})
            cds_db = cds_db.join(id_name_map, on = 'taxon_id', how = 'left', maintain_order = "left")
        # If not, use filenames
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
    # Process arguments
    args = parse_arguments()
    cores = args.cores
    input_path = args.input
    output_path = args.output
    parsing_mode = args.mode
    use_taxa = args.use_taxa
    if args.gzip:
        gzip = "gzip"
    else:
        gzip = "uncompressed"
        
    # Parse the input files
    cds_db = parse_inputs(input_path, parsing_mode, cores)
    
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

