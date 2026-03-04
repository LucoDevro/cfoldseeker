#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import time
import sys
import requests
import json
import io
import logging
import polars as pl
import itertools as it
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from kegg_pull.pull import MultiProcessMultiplePull
from Bio import SeqIO

from cfoldseeker.classes import Hit, Search, _sanitise_hit_attr

LOG = logging.getLogger(__name__)


class RemoteSearch(Search):
    
    def __init__(self, query, mapping_table_path, params = {}, hits = [], clusters = [], 
                 output_folder = Path('.'), temp_folder = Path('.')):
        
        super().__init__(query, params, hits, clusters, output_folder, temp_folder)
        
        LOG.debug(f'Scanning ID mapping table from {str(mapping_table_path)}')
        self.mapping_table: pl.LazyFrame = pl.scan_csv(mapping_table_path, has_header = False, separator = "\t",
                                                       new_columns = ['Uniprot', 'DB', 'ID'])
                
        return None
    
    def __repr__(self) -> str:
        return f"Remote Search of {','.join(list(self.query.keys()))} in {self.params['db']} with {len(self.clusters)} clusters identified"
    
    
    def run_foldseek(self) -> None:
        """
        Submits queries to the FoldSeek webserver and collects the results.
        """

        FOLDSEEK_SUBMISSION_URL = "https://search.foldseek.com/api/ticket"
        FOLDSEEK_RESULTS_URL = "https://search.foldseek.com/api/result"
        
        """
        Submits one structure file to the FoldSeek API and returns the submission ticket in dictionary form.
        """
        def submit_foldseek_query(query_path: Path, dbs: list, taxfilters: list) -> None:
            with open(query_path, "rb") as f:
                files = {"q": f}
                data = [("mode", "3diaa")]
                for db in dbs:
                    data.append(("database[]", db))
                for taxfilt in taxfilters:
                    data.append(('taxfilter', taxfilt))
                LOG.debug(f"Posting request on FoldSeek webserver {FOLDSEEK_SUBMISSION_URL}")
                LOG.debug(f"using the following parameters: {data}")
                response = requests.post(FOLDSEEK_SUBMISSION_URL, files=files, data=data)
                if response.status_code == 200:
                    LOG.info(f'Query {query_path} successfully submitted!')
                    return response.json()
                else:
                    LOG.exception(f"Error submitting query {query_path}!")
                    
            return None
        
        """
        Checks the status of the job associated with the provided job ID. Returns the status flag.
        """
        def check_query_status(job_id: str) -> str:
            url = f"{FOLDSEEK_SUBMISSION_URL}/{job_id}"
            LOG.info(f'Checking out status of job {job_id}')
            results = requests.get(url).json()
            status = results['status']
            
            return status
        
        """
        Waits for the FoldSeek job associated with the provided job ID to complete, and retrieves its result.
        Returns the parsed results in a dictionary.
        """
        def retrieve_foldseek_results(job_id: str) -> dict:
            while True:
                status = check_query_status(job_id)
                if status == "COMPLETE":
                    LOG.info(f"Job {job_id} has completed! Downloading results...")
                    entry = 0
                    url = f"{FOLDSEEK_RESULTS_URL}/{job_id}/{entry}"
                    results = requests.get(url).json()
                    break
                else:
                    LOG.info(f"Job {job_id} has not completed yet. Waiting another 10 seconds...")
                    time.sleep(10)
                    
            return results
            
        
        # First submit all query proteins to FoldSeek
        LOG.info('Posting queries at FoldSeek webserver')
        query_foldseek = lambda x: submit_foldseek_query(x, self.params['db'], self.params['taxfilters'])
        with ThreadPoolExecutor(max_workers = self.params['max_workers']) as executor:
            tickets = dict(zip(self.query.keys(), executor.map(query_foldseek, self.query.values())))
                
        # Then wait for the results and retrieve them automatically when completed
        all_job_ids = [ticket['id'] for ticket in tickets.values()]
        with ThreadPoolExecutor(max_workers = self.params['max_workers']) as executor:
            all_results = dict(zip(self.query.keys(), executor.map(retrieve_foldseek_results, all_job_ids)))
        
        LOG.info('FINISHED PART 1')
        LOG.info('All queries have been processed and downloaded!')
        self.hits = all_results
        
        LOG.info('Saving results in temporary files')
        for query,result in all_results.items():
            with open((self.TEMP_DIR / f"foldseek_result_{query}").with_suffix('.json'), "w") as handle:
                json.dump(result, handle)
        
        return None
    
    
    def parse_foldseek_results(self) -> None:
        """
        Parses the FoldSeek results and creates Hit objects for the passing hits
        """
        
        ## Load the thresholds from params
        max_eval = self.params['max_eval']
        min_score = self.params['min_score']
        min_seqid = self.params['min_seqid']
        min_qcov = self.params['min_qcov']
        min_tcov = self.params['min_tcov']
        
        LOG.debug('Applying the following hit filtering thresholds:')
        LOG.debug(f'bitscore >= {min_score}')
        LOG.debug(f'query coverage >= {min_qcov}')
        LOG.debug(f'target coverage >= {min_tcov}')
        LOG.debug(f'minimum sequence identity threshold >= {min_seqid}')
        LOG.debug(f'maximum evalue threshold <= {max_eval}')
        
        ## Parse the FoldSeek results
        all_hits = []
        for query, task in self.hits.items():
            LOG.debug(f"Parsing FoldSeek results for {query}")
            for results_by_db in task['results']:
                db = results_by_db['db']
                for hit_entry in results_by_db['alignments'][0]:
                    target = hit_entry['target'] # uniprot ID and name
                    if 'afdb' in db:
                        db_id = target.split('-')[1]
                    elif db == 'pdb100':
                        db_id = target[:4].upper() + '.' + target[22]
                    else:
                        LOG.error(f"Unsupported DB ({db}) for ID {target}")
                        continue
                    name = ' '.join(target.split(' ')[1:])
                    taxon_name = hit_entry['taxName'] # taxon name
                    taxon_id = hit_entry['taxId'] # taxon ID
                    evalue = float(hit_entry['eval']) # FoldSeek e-value
                    score = int(hit_entry['score']) # FoldSeek hit score
                    seqid = float(hit_entry['seqId']) # Sequence identity with the query protein
                    qcov = (int(hit_entry['qEndPos']) - int(hit_entry['qStartPos'])) / int(hit_entry['qLen']) * 100 # Query coverage
                    tcov = (int(hit_entry['dbEndPos']) - int(hit_entry['dbStartPos']) / int(hit_entry['dbLen'])) * 100 # Target coverage
                    
                    # Create Hit object and collect it if it passes all thresholds
                    if all((evalue <= max_eval, 
                            score >= min_score, 
                            seqid >= min_seqid, 
                            qcov >= min_qcov, 
                            tcov >= min_tcov)):
                        hit = Hit(db_id, query, name = name, taxon_name = taxon_name, taxon_id = taxon_id, db = db,
                                  evalue = evalue, score = score, seqid = seqid, qcov = qcov, tcov = tcov)
                        all_hits.append(hit)
                        
        if len(all_hits) == 0:
            LOG.error("No hits identified by FoldSeek!")
            sys.exit()
            
        ## Filter out redundant hit instances of hits found in multiple databases
        # Group by Uniprot ID and DB to find potentially redundant hits
        LOG.debug('Removing redundant hits found in multiple DBs')
        hit_db_groups = {}
        for hit in all_hits:
            if (hit.db_id, hit.db) in hit_db_groups.keys():
                hit_db_groups[(hit.db_id, hit.db)].append(hit)
            else:
                hit_db_groups[(hit.db_id, hit.db)] = [hit]
        # Keep the first DB hit instance for every Uniprot-DB group
        filtered_hits = []
        all_hit_ids = {h.db_id for h in all_hits}
        for hit_id in all_hit_ids:
            hits_this_hit_id = [g for g in hit_db_groups.keys() if g[0] == hit_id]
            db_hits_to_keep = hits_this_hit_id[0]
            filtered_hits.append(hit_db_groups[db_hits_to_keep])
                
        # Flatten out and save
        filtered_hits = list(it.chain(*filtered_hits))
        self.hits = filtered_hits
        
        LOG.info('FINISHED PART 2')
        LOG.info(f'Found {len(filtered_hits)} gene hits.')
        
        return None
            
                
    def crossref_afdb(self) -> list:
        """
        Crossrefs all hits in the various supported AFDB flavours where possible and
        retrieves the genomic neighbourhood information (scaffold IDs and coordinates).
        """
        
        def prepare_mapping_dict(all_uniprot_ids: list, db: str) -> dict:
            """
            Effectively extracts a mapping dictionary from the full Uniprot crossref LazyFrame, linking the Uniprot ID (keys)
            with a list of one or more multiple crossrefs in the provided target database (values).
            """
            LOG.debug(f"Preparing a mapping dictionary from UniProt to {db}")
            
            res = self.mapping_table.filter(pl.col('Uniprot').is_in(all_uniprot_ids)
                                            ).filter(pl.col('DB') == db).drop('DB')
            res = res.group_by('Uniprot').all()
            res = res.collect() # Materialise the LazyFrame into a DataFrame
            res = dict(zip(res['Uniprot'], res['ID'].to_list()))
            
            return res
        

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
                        

        def pull_from_ena(entry: str) -> None:
            """
            Pulls a GenPept record from the ENA Browser API.
            """
            ENA_BROWSER_URL = "https://www.ebi.ac.uk/ena/browser/api/embl"
            
            url = f"{ENA_BROWSER_URL}/{entry}"
            LOG.debug(f'Going to pull GenPept record from {url}')
            response = requests.get(url)
            if response.status_code == 200:
                return response.text
            else:
                LOG.warning(f'Error pulling GenPept record {entry}')
            
            return None
            
               
        ### Identify the AFDB hits. These will be crossreffed by this function
        all_afdb_hits = [h for h in self.hits if 'afdb' in h.db]
            
        ### Get the crossreffing IDs            
        ## First, try to get KEGG IDs
        LOG.info('Searching crossrefs in KEGG')
        all_uniprot_ids = list({h.db_id for h in all_afdb_hits})
        hits_failed_kegg = []
        # Extract a mapping table for KEGG IDs from the Uniprot crossref mapping table
        all_uniprot_kegg = prepare_mapping_dict(all_uniprot_ids, 'KEGG')
        # Fill all KEGG IDs
        for h in all_afdb_hits:
            if h.db_id in all_uniprot_kegg.keys():
                h.crossref_id = all_uniprot_kegg[h.db_id]
                h.crossref_method = "KEGG"
            else:
                hits_failed_kegg.append(h)
                continue
        LOG.info(f'{len(hits_failed_kegg)} hits have no crossref in KEGG.')
        
        ## Then, try to get GenPept IDs for the hits that did not get a KEGG ID.
        LOG.info('Searching crossrefs in GenPept')
        uniprot_ids_failed_kegg = list({h.db_id for h in hits_failed_kegg})
        hits_failed_genpept = []
        all_uniprot_genpept = prepare_mapping_dict(uniprot_ids_failed_kegg, 'EMBL-CDS')
        # Fill the GenPept IDs
        for h in all_afdb_hits:
            if h.db_id in all_uniprot_genpept.keys():
                genpept_id = all_uniprot_genpept[h.db_id]
                genpept_id = [i for i in genpept_id if i != '-'] # Filter out empty crossrefs (often mRNA records)
                h.crossref_id = genpept_id
                h.crossref_method = "GenPept"
            else:
                hits_failed_genpept.append(h)
                continue
        LOG.info(f'..., of which {len(hits_failed_genpept)} hits have no crossref in GenPept.')
        print([h.db_id for h in hits_failed_genpept])
            
        ## Split records with multiple crossrefs and discard the ones without crossref
        LOG.debug("Sanitising the hit list")
        all_afdb_hits = _sanitise_hit_attr(all_afdb_hits, 'crossref_id')
        
        ### Do the actual crossreffing
        LOG.info('Navigating the crossrefs')
        processed_hits = []
        ## Get the scaffold and genome positions for the KEGG crossreffing method
        ## Pull all KEGG Gene and Genome records
        LOG.info("Pulling KEGG crossrefs")
        pull = MultiProcessMultiplePull(n_workers = self.params['max_workers'])
        
        # KEGG Gene IDs
        all_kegg_gene_ids = list({h.crossref_id for h in all_afdb_hits if h.crossref_method == 'KEGG'})
        LOG.debug(f'Going to pull {len(all_kegg_gene_ids)} KEGG Gene entries')
        _, gene_records = pull.pull_dict(all_kegg_gene_ids)
        # KEGG Genome IDs
        all_kegg_genome_ids = list({'genome:' + i.split(':')[0] for i in all_kegg_gene_ids})
        LOG.debug(f'Going to pull {len(all_kegg_genome_ids)} KEGG Genome entries')
        _, genome_records = pull.pull_dict(all_kegg_genome_ids)
        genome_records = {k.split(':')[-1]: v for k,v in genome_records.items()}
        
        ## Extract scaffold and coordinate information from the pulled KEGG records
        LOG.debug('Parsing downloaded KEGG entries')
        genomic_info = {k: extract_genomic_information_kegg(v) for k,v in gene_records.items()}
        scaffold_mappings = {k: extract_scaffold_mapping_kegg(v) for k,v in genome_records.items()}
        # Pass the extracted information on to the Hit objects
        LOG.debug('Extract CDS coordinates')
        for h in all_afdb_hits:
            try:
                position_info = genomic_info[h.crossref_id]
                scaffold_info = scaffold_mappings[h.crossref_id.split(':')[0]]
            except KeyError:
                continue
            if len(position_info) == 0 or len(scaffold_info) == 0:
                continue
            h.coords = position_info['coords']
            h.strand = position_info['strand']
            h.scaff = scaffold_info[position_info['scaffold']]
            
            processed_hits.append(h)
                
        ## Get the scaffold and genome positions for the GenPept crossreffing method
        ## Pull all GenPept records from ENA in EMBL format
        LOG.info("Pulling GenPept crossrefs from ENA")
        all_genpept_gene_ids = list({h.crossref_id for h in all_afdb_hits if h.crossref_method == "GenPept"})
        LOG.debug(f'Going to pull {len(all_genpept_gene_ids)} ENA GenPept entries')
        with ThreadPoolExecutor(max_workers = self.params['max_workers']) as executor:
            ena_records = dict(zip(all_genpept_gene_ids, executor.map(pull_from_ena, all_genpept_gene_ids)))
            
        ## Extract scaffold and coordinate information from the pulled ENA records
        LOG.debug("Parsing downloaded GenPept entries")
        genomic_info = {k: extract_genomic_information_ena(v) for k,v in ena_records.items()}
        # Pass the extracted information on the Hit objects
        LOG.debug('Extract CDS coordinates')
        for h in all_afdb_hits:
            try:
                position_info = genomic_info[h.crossref_id]
                if position_info == None:
                    continue
                h.scaff = position_info['scaffold']
                h.coords = position_info['coords']
                h.strand = position_info['strand']
            except KeyError:
                continue
            
            processed_hits.append(h)
        
        ### Update the hit set
        all_afdb_hits = processed_hits
        LOG.info("FINISHED PART 2B")
        LOG.info(f'{len(processed_hits)} hits have been processed.')
        
        return all_afdb_hits
    

    def run(self) -> None:
        """
        Complete remote search workflow run
        """
        
        LOG.info("STARTING PART 1: Executing FoldSeek search")
        self.run_foldseek()
        
        LOG.info("STARTING PART 2: Parsing FoldSeek results")
        self.parse_foldseek_results()
        
        LOG.info('STARTING PART 2B: Fetching CDS coordinates via crossreffing')
        afdb_hits = self.crossref_afdb()
        self.hits = afdb_hits
        
        LOG.info('STARTING PART 3: Identifying gene clusters')
        self.identify_clusters()
        
        return None
        
    