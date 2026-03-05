#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import io
import logging
from Bio import SeqIO

LOG = logging.getLogger(__name__)


def extract_genomic_information_kegg(gene_entry: str) -> dict:
    """
    Extracts the genomic information from a pulled KEGG Gene record.
    """
    # Genomic positions are at the POSITION line
    position_line = [line for line in gene_entry.split('\n') if 'POSITION' in line]
    
    position_info = {}
    if len(position_line) == 0:
        return position_info
    else:
        position_line = position_line[0]
        # If there is a scaffold mentioned, get it
        if ':' in position_line:
            internal_scaffold_id, coords = position_line.split(':')
            internal_scaffold_id = internal_scaffold_id[12:]
        # If not, leave empty. The downstream processing will handle this
        else:
            internal_scaffold_id, coords = '', position_line
        # Extract the coordinates of all exons, ignoring indefinite boundaries
        coords = coords.translate(str.maketrans('', '', '<>'))
        coord_groups = re.findall(r'\d+\.\.\d+', coords)
        # If no coordinate information, return the response dictionary empty
        if len(coord_groups) == 0:
            return position_info
        # Else, parse it
        coord_groups = [i.split('..') for i in coord_groups]
        coord_groups = [[int(j) for j in i] for i in coord_groups]
        if 'complement' in coords:
            strand = "-"
        else:
            strand = "+"
        
        # Gather the results
        position_info['scaffold'] = internal_scaffold_id
        position_info['coords'] = coord_groups
        position_info['strand'] = strand
        
        return position_info
    

def extract_scaffold_mapping_kegg(genome_entry: str) -> dict:
    """
    Maps all KEGG scaffold IDs for a Genome entry to the associated GenBank/RefSeq IDs.
    """
    lines = genome_entry.split('\n')
    ## First the CHROMOSOME field
    # Find the start
    start_chromosome = [idx for idx,line in enumerate(lines) if 'CHROMOSOME' in line]
    if len(start_chromosome) == 0:
        mapping_scaffolds = {}
    else:
        # Expand it
        index = start_chromosome[0]
        chromosome_field = [lines[index][12:]]
        while index < len(lines) and lines[index+1].startswith(' '):
            index += 1
            chromosome_field.append(lines[index][12:])
        # Parse it
        internal_scaffold_ids = [re.split(r'[;\s]', l)[0] for l in chromosome_field]
        internal_scaffold_ids = ['' if l == 'Circular' else l for l in internal_scaffold_ids]
        scaffolds = [l.split(':')[1][:-1] for l in chromosome_field]
        mapping_scaffolds = dict(zip(internal_scaffold_ids, scaffolds))
    
    ## Then the PLASMIDS field
    # Find the start
    start_plasmid = [idx for idx,line in enumerate(lines) if 'PLASMID' in line]
    if len(start_plasmid) == 0:
        mapping_plasmids = {}
    else:
        # Expand it
        index = start_plasmid[0]
        plasmid_field = [lines[index][12:]]
        while index < len(lines) and lines[index+1].startswith(' '):
            index += 1
            plasmid_field.append(lines[index][12:])
        # Parse it
        internal_plasmid_ids = [re.split(r'[;\s]', l)[0] for l in plasmid_field]
        plasmids = [l.split(':')[1][:-1] for l in plasmid_field]
        mapping_plasmids = dict(zip(internal_plasmid_ids, plasmids))
    
    ## Wrap it in a dictionary
    mapping = mapping_scaffolds | mapping_plasmids
    
    return mapping


def extract_genomic_information_ena(record: str) -> dict:
    """
    Extracts the genomic information from a pulled ENA GenPept record.
    """
    position_info = {}
    # catch for empty or bad results
    if record == None:
        return None
    embl = io.StringIO(record)
    seq_record = list(SeqIO.parse(embl, format = 'embl'))[0]
    cds = [f for f in seq_record.features if f.type == 'CDS'][0]
    # Genome coordinates
    parts = cds.location.parts
    coord_groups = [[int(p.start)+1, int(p.end)] for p in parts] # start coordinate is one off in BioPython parsing
    # Scaffold
    scaffold = list({p.ref for p in parts})[0]
    # Strand
    strand = cds.location.strand
    if strand == 1:
        strand = '+'
    elif strand == -1:
        strand = '-'
    
    # collect
    position_info['coords'] = coord_groups
    position_info['strand'] = strand
    position_info['scaffold'] = scaffold
    
    return position_info

