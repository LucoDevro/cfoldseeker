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
        prog = 'build_db.py',
                epilog = 
                """
                Lucas De Vrieze, Miguel Biltjes
                (c) 2026 Masschelein lab, VIB
                """,
                formatter_class = argparse.RawDescriptionHelpFormatter,
                description = 
                """
                Helper tool to construct a CDS coordinates DB directly from a folder of NCBI GFF files.
                """,
                add_help = False
                )
    
    parser.add_argument('-i', '--input', dest = 'input', type = Path, default = Path('.'), help = "Path to folder containing the NCBI GFF files (default: current directory)")
    parser.add_argument('-c', '--cores', dest = 'cores', type = int, default = 1, help = "Number of cores available to use (default: 1).")
    parser.add_argument('-o', '--output', dest = 'output', type = Path, default = Path('local_db'), help = "Filepath to save CDS coordinate DB (default: local_db).")
    parser.add_argument('-gz', '--gzip', dest = 'gzip', default = False, action = 'store_true', help = "Gzip output (default: False).")
    parser.add_argument('-f', '--force', dest = 'force', default = False, action = 'store_true', help = "Force overwriting output (default: false).")
    parser.add_argument('-h', '--help', action = 'help', help = "Show this help message and exit")      

    args = parser.parse_args()
    
    assert args.input.exists() and args.input.is_dir() and any(args.input.glob('*.gff')), 'Input folder path does not exist or does not contain GFF files.'
    if args.output.exists():
        if args.force:
            LOG.warning("Output already exists, but it will be overwritten.")
        else:
            LOG.error("Output already exists! Rerun with -f to overwrite it.")
            sys.exit()
    else:
        args.output.parent.mkdir(parents = True, exist_ok = True)
    
    return args


def parse_one_gff(path: Path) -> pl.DataFrame:
    """
    Parse one NCBI GFF table into a Polars DataFrame
    """
    df = pl.scan_csv(path, separator = "\t", has_header = False, comment_prefix = '#',
                     new_columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
                     ).select(['seqid', 'type', 'start', 'end', 'strand', 'attributes'])
    
    df_regions = df.filter(pl.col('type') == 'region')
    df_cds = df.filter(pl.col('type') == 'CDS')
    
    # Extract taxon ID
    taxon_id = df_regions.select(pl.col('attributes').str.extract(r'taxon:([0-9]+)')).head(1)
    taxon_id = taxon_id.collect().rows()[0][0]
    
    # Populate the CDS DB
    df_cds = df_cds.drop('type')
    cds_record = df_cds.select([pl.col('attributes').str.extract(r'protein_id=([^;]+)').alias('gene_tag'),
                            pl.col('attributes').str.extract(r'product=([^;]+)').alias('name'),
                            pl.col('seqid').alias('contig'),
                            pl.concat_str(['start', 'end'], separator = '..').alias('coords'),
                            pl.col('strand'),
                            pl.lit(taxon_id).alias('taxon_id')
                            ])
    cds_record = cds_record.drop_nulls(subset = 'gene_tag')
    
    # Aggregate multiple CDSes (i.e. exons) in one record
    cds_record = cds_record.group_by(pl.all().exclude('coords')).agg(pl.col('coords').str.join(','))
    
    return cds_record.collect()


def main():
    # Process arguments
    args = parse_arguments()
    cores = args.cores
    input_path = args.input
    output_path = args.output
    if args.gzip:
        gzip = "gzip"
    else:
        gzip = "uncompressed"
    
    # Parse all GFFs
    LOG.info('Parsing all GFF files')
    with ThreadPoolExecutor(max_workers = cores) as executor:
        parsed_gffs_to_concat = executor.map(parse_one_gff, input_path.glob('*.gff'))
        
    LOG.info('Construct CDS coordinates DB')
    cds_db = pl.concat(parsed_gffs_to_concat)
    
    # Fetch all taxon names
    LOG.info('Fetch taxon names using NCBI Entrez')
    all_taxon_ids = cds_db.select('taxon_id').unique().to_series().to_list()
    with Entrez.esummary(db = 'taxonomy', id = all_taxon_ids) as handle:
        records = list(Entrez.read(handle))
    all_taxon_names = [str(i['ScientificName']) for i in records]
    
    # Join with the CDS DB
    LOG.info('Add taxon name column')
    id_name_map = pl.DataFrame({'taxon_id': all_taxon_ids, 'taxon_name': all_taxon_names})
    cds_db = cds_db.join(id_name_map, on = 'taxon_id', how = 'left', maintain_order = "left")
    
    # Write results
    LOG.info('Write DB to disk')
    cds_db.write_csv(output_path, separator = '\t', include_header = True, compression = gzip)
    
    
if __name__ == "__main__":
    main()

