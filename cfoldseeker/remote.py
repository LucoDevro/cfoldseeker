#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import json
import logging
import warnings
import polars as pl
import itertools as it
from pathlib import Path
from copy import deepcopy
from tqdm.contrib.concurrent import thread_map
from tqdm.contrib.logging import logging_redirect_tqdm
from kegg_pull.pull import MultiProcessMultiplePull
from Bio import Entrez

from cfoldseeker.classes import Hit, Search
from cfoldseeker.communication import (submit_foldseek_query, retrieve_foldseek_results, 
                                       pull_dict_from_unisave, pull_from_ena)
from cfoldseeker.remote_parsers import (extract_genomic_information_ena, 
                                        extract_genomic_information_kegg, 
                                        extract_scaffold_mapping_kegg)


LOG = logging.getLogger(__name__)


class RemoteSearch(Search):
    """
    Subclass executing the workflow for gene cluster identification from remote protein searches.
    
    Extends the Search base class to perform FoldSeek-based searches against remote
    databases, parse results, and cross-reference hits to genomic data.
    Uses a local copy of the UniProt mapping table to retrieve genomic coordinates for each protein
    with fallback strategies using KEGG and ENA.
    
    New attributes:
        mapping_table: Polars LazyFrame containing UniProt ID cross-references to
            various databases (KEGG, EMBL-CDS, etc.).
            
    Inherits from:
        Search: Base class providing the cluster identification and output generation capabilities.
        
    See also:
        Local: Sister class providing the search and parsing capabilities for local database searches
    """
    
    def __init__(self, query, mapping_table_path, params = {}, hits = [], clusters = [], 
                 output_folder = Path('.'), temp_folder = Path('.')):
        """
        Initialise a RemoteSearch instance.
        
        Sets up the RemoteSearch object by calling the parent Search class
        initialiser and scanning the UniProt ID mapping table from a TSV file.
        
        Args:
            query: Dictionary mapping query protein names to file paths.
            mapping_table_path: Path to a TSV file containing UniProt cross-references.
                Expected format: three columns (Uniprot, DB, ID) without header.
            params: Dictionary of search parameters and configuration options.
                Defaults to an empty dictionary.
            hits: List of Hit objects from previous searches. Defaults to empty list.
            clusters: List of identified gene clusters. Defaults to empty list.
            output_folder: Path object for final output directory. Defaults to current
                directory ('.').
            temp_folder: Path object for temporary file directory. Defaults to current
                directory ('.').
        
        Returns:
            None
        """
        
        super().__init__(query, params, hits, clusters, output_folder, temp_folder)
        
        LOG.debug(f'Scanning ID mapping table from {str(mapping_table_path)}')
        self.mapping_table: pl.LazyFrame = pl.scan_csv(mapping_table_path, has_header = False, separator = "\t",
                                                       new_columns = ['Uniprot', 'DB', 'ID'])
                
        return None
    
    
    def __repr__(self) -> str:
        return f"Remote Search of {','.join(list(self.query.keys()))} in {self.params['db']} with {len(self.clusters)} clusters identified"
    
    
    def _sanitise_hit_attr(self, hits: list[Hit], attr: str) -> list[Hit]:
        """
        Sanitise the hit list for the provided attribute.
        
        Sanitising means standardising the hit table. In case of multiple values for a given attribute,
        it splits a hit instance into multiple hit instances that are identical except for that exact
        crossref attribute. The split hits are first introduced in a sublist replacing the attribute value
        in the old hit instance. Then the full hit list is flattened. In case of an absent attribute, the hit is discarded.
        
        Args:
            hits (list[Hit]): hits that must be sanitised.
            attr (str): name of the attribute that must be sanitised in each hit
            
        Mutates:
            Hit objects in list: split or unpack attribute values in case of multiple ones
            hits: empty attribute values make an instance be discarded from the list
        
        Returns:
            hits (list[Hit]): sanitised hits
        """
        sanitised_hits = []
        split_counter = 0
        discard_counter = 0
        untouched_counter = 0
        
        # Loop over all hits, but keep track of the index
        for h in hits:
            # Determine whether we should split a record (# attribute values > 1)
            nb_copies = len(getattr(h, attr))
            # If we need to split a record, make as much copies as there are attribute values 
            # and distribute these over each split instance's attribute value
            if nb_copies > 1:
                split_counter += 1
                split_records = []
                copy_counter = 0
                while copy_counter < nb_copies:
                    new_record = deepcopy(h)
                    setattr(new_record, attr, getattr(h, attr)[copy_counter])
                    split_records.append(new_record)
                    copy_counter += 1
                    
                # Insert the list of the split hit instances at the index of the old record
                sanitised_hits.append(split_records)
                
            # If we don't need to split the record, just unpack the singleton list
            elif nb_copies == 1:
                untouched_counter += 1
                setattr(h, attr, getattr(h, attr)[0])
                
                # Insert the unpacked list at the index of the old record
                sanitised_hits.append([h])
                
            # If the attribute is empty, don't include it in the new hit set
            else:
                discard_counter += 1
                continue
        
        # Flatten the full hit list again
        hits = list(it.chain(*sanitised_hits))
        
        LOG.debug(f"{split_counter} hit records have been exploded while sanitising.")
        LOG.debug(f"{discard_counter} hit records have been discarded while sanitising.")
        LOG.debug(f"{untouched_counter} hit records have been left untouched while sanitising.")
        
        return hits
    
    
    def run_foldseek(self) -> None:
        """
        Submit protein queries to FoldSeek and retrieve the results.
        
        Sends all query structures to the FoldSeek webserver in parallel, collects
        submission tickets, monitors job status, and downloads completed results.
        Stores raw JSON results in the temporary folder.
        
        Returns:
            None
        """

        # First submit all query proteins to FoldSeek
        LOG.info('Submitting queries to FoldSeek webserver')
        query_foldseek = lambda x: submit_foldseek_query(x, self.params['db'], self.params['taxfilters'])
        query_paths = list(self.query.values())
        query_labels = list(self.query.keys())
        
        with logging_redirect_tqdm(loggers = [LOG]):
            runs = thread_map(query_foldseek, query_paths, 
                              max_workers = self.params['max_workers'],
                              leave = False,
                              disable = self.params['no_progress'])
            tickets = dict(zip(query_labels, runs))
        all_job_ids = [ticket['id'] for ticket in tickets.values()]
                
        # Then wait for the results and retrieve them automatically when completed
        LOG.info('Checking status and retrieving results when ready')
        with logging_redirect_tqdm(loggers = [LOG]):
            runs_results = thread_map(retrieve_foldseek_results, all_job_ids,
                                      max_workers = self.params['max_workers'],
                                      leave = False,
                                      disable = self.params['no_progress'])
            all_results = dict(zip(self.query.keys(), runs_results))
        
        LOG.info('All queries have been processed and downloaded!')
        self.hits = all_results
        
        LOG.info('Saving results in temporary folder')
        for query,result in all_results.items():
            with open((self.TEMP_DIR / f"foldseek_result_{query}").with_suffix('.json'), "w") as handle:
                json.dump(result, handle)
        
        return None
    
    
    def passes_criteria(self, hit: Hit):
        """
        Check if a hit passes the criteria set for this search.
        
        Returns:
            (bool): Dit the hit pass all criteria?
        """
        criteria_passes = (hit.evalue <= self.params['max_eval'],
                           hit.score >= self.params['min_score'],
                           hit.seqid >= self.params['min_seqid'],
                           hit.qcov >= self.params['min_qcov'],
                           hit.tcov >= self.params['min_tcov'])
        
        return all(criteria_passes)
    
    
    def parse_foldseek_results(self) -> None:
        """
        Parse FoldSeek results and create Hit objects for hits passing the predefined criteria.
        
        Extracts hit information from the raw FoldSeek results, and applies filtering
        thresholds based on e-value, bit score, sequence identity, and coverage.
        Creates Hit objects for hits meeting all criteria and removes redundant
        hits found in multiple databases, keeping only the first instance per
        UniProt ID.
        
        Returns:
            None
            
        Raises:
            RuntimeError: If the hit list is empty after applying the criteria
            
        Note:
            Logs a warning when there are no hits for a certain query and DB pair.
        """
        
        LOG.debug('Applying the following hit filtering thresholds:')
        LOG.debug(f'bitscore >= {self.params["min_score"]}')
        LOG.debug(f'query coverage >= {self.params["min_qcov"]}')
        LOG.debug(f'target coverage >= {self.params["min_tcov"]}')
        LOG.debug(f'minimum sequence identity threshold >= {self.params["min_seqid"]}')
        LOG.debug(f'maximum evalue threshold <= {self.params["max_eval"]}')
        
        ## Parse the FoldSeek results
        all_hits = []
        # Each query requires parsing a separate result
        for query, task in self.hits.items():
            LOG.debug(f"Parsing FoldSeek results for {query}")
            # In each result, multiple databases may have been queried
            for results_by_db in task['results']:
                db = results_by_db['db']
                # Since every query has been submitted separately there is only one set of alignments to parse per result set
                alignments = results_by_db['alignments']
                
                # Warn if no hit has been found in this DB
                if len(alignments) == 0:
                    LOG.warning(f'No hit has been found for query {query} in database {db}')
                    continue
                
                for hit_entry in alignments[0]:
                    target = hit_entry['target'] # uniprot ID and name
                    if 'afdb' in db:
                        db_id = target.split('-')[1]
                    elif db == 'pdb100':
                        db_id = target[:4].upper() + '.' + target[22]
                    else:
                        LOG.error(f"Unsupported DB ({db}) for ID {target}")
                        continue
                    
                    # Gather the information for each hit
                    name = ' '.join(target.split(' ')[1:])
                    taxon_name = hit_entry['taxName'] # taxon name
                    taxon_id = hit_entry['taxId'] # taxon ID
                    evalue = float(hit_entry['eval']) # FoldSeek e-value
                    score = int(hit_entry['score']) # FoldSeek hit score
                    seqid = float(hit_entry['seqId']) # Sequence identity with the query protein
                    qcov = (int(hit_entry['qEndPos']) - int(hit_entry['qStartPos'])) / int(hit_entry['qLen']) * 100 # Query coverage (in %)
                    tcov = (int(hit_entry['dbEndPos']) - int(hit_entry['dbStartPos'])) / int(hit_entry['dbLen']) * 100 # Target coverage (in %)
                    
                    # Create Hit object and add it if it passes all criteria
                    hit = Hit(db_id, query, name = name, taxon_name = taxon_name, taxon_id = taxon_id, db = db,
                                  evalue = evalue, score = score, seqid = seqid, qcov = qcov, tcov = tcov)
                    if self.passes_criteria(hit):
                        all_hits.append(hit)
        
        ## Stop if hit list is empty
        if len(all_hits) == 0:
            msg = "No hits pass the criteria!"
            LOG.error(msg)
            raise RuntimeError(msg)
            
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
            hits_this_hit_id = [g for g in hit_db_groups.keys() if g[0] == hit_id] # all hits with this hit ID
            db_hits_to_keep = hits_this_hit_id[0] # Keep the first one
            filtered_hits.append(hit_db_groups[db_hits_to_keep])
                
        # Flatten out and save
        filtered_hits = list(it.chain(*filtered_hits))
        self.hits = filtered_hits
        
        LOG.info(f'Found {len(filtered_hits)} gene hits.')
        
        return None
    
    
    def prepare_mapping_dict(self, all_uniprot_ids: list, db: str) -> dict:
        """
        Extract a mapping dictionary from the local UniProt mapping table.
        
        Filters the full UniProt cross-reference LazyFrame to extract IDs
        for a specific target database, grouping multiple cross-references
        per UniProt ID into a dictionary of lists.
        
        Args:
            all_uniprot_ids: List of UniProt accession numbers to extract.
            db: Target database name (e.g., 'KEGG', 'EMBL-CDS').
        
        Returns:
            A dictionary mapping UniProt IDs (keys) to lists of cross-reference
            IDs in the target database (values).
        """
        LOG.debug(f"Preparing a mapping dictionary from UniProt to {db}")
        
        mapping = self.mapping_table.filter(pl.col('Uniprot').is_in(all_uniprot_ids)
                                        ).filter(pl.col('DB') == db).drop('DB')
        mapping = mapping.group_by('Uniprot').all()
        mapping = mapping.collect() # Materialise the LazyFrame into a DataFrame
        mapping = dict(zip(mapping['Uniprot'], mapping['ID'])) # Create a mapping dictionary
        
        return mapping
    
    
    def crossref_afdb_via_kegg(self, afdb_hits: list[Hit]) -> list[Hit]:
        """
        Cross-reference AFDB hits to KEGG IDs using the UniProt mapping table
        
        Retrieves the KEGG IDs for all hits and updates the the crossref ID and
        method if any.
        
        Args:
            afdb_hits: List of Hit objects without genomic information and cross-references
            
        Mutates:
            Hit objects in afdb_hits: Fills cross-reference ID and method attributes
        
        Returns:
            hits_failed_legg (list[Hit]): list of Hits that have no cross-reference to KEGG in 
            the UniProt mapping table. These need to be processed by a different cross-referencing
            method.
        """
        all_uniprot_ids = list({h.db_id for h in afdb_hits})
        hits_failed_kegg = []
        
        # Extract a mapping table for KEGG IDs from the Uniprot crossref mapping table
        all_uniprot_kegg = self.prepare_mapping_dict(all_uniprot_ids, 'KEGG')
        
        # Fill all KEGG IDs
        for h in afdb_hits:
            if h.db_id in all_uniprot_kegg.keys():
                h.crossref_id = all_uniprot_kegg[h.db_id]
                h.crossref_method = "KEGG"
            else:
                hits_failed_kegg.append(h)
        
        return hits_failed_kegg
    
    
    def crossref_afdb_via_genpept(self, afdb_hits: list[Hit]) -> list[Hit]:
        """
        Cross-reference AFDB hits to GenPept IDs in ENA using the UniProt mapping table
        
        Retrieves the GenPept IDs in ENA for all hits and updates the the crossref ID and
        method if any.
        
        Args:
            afdb_hits: List of Hit objects without genomic information and cross-references
            
        Mutates:
            Hit objects in afdb_hits: Fills cross-reference ID and method attributes
        
        Returns:
            hits_failed_legg (list[Hit]): list of Hits that have no cross-reference to GenPept in 
            the UniProt mapping table. These need to be processed by a different cross-referencing
            method.
        """
        all_uniprot_ids = list({h.db_id for h in afdb_hits})
        hits_failed_genpept = []
        
        # Extract a mapping table for GenPept IDs from the Uniprot crossref mapping table
        all_uniprot_genpept = self.prepare_mapping_dict(all_uniprot_ids, 'EMBL-CDS')
        
        # Fill the GenPept IDs
        for h in afdb_hits:
            if h.db_id in all_uniprot_genpept.keys():
                genpept_id = all_uniprot_genpept[h.db_id]
                
                # Filter out empty crossrefs (often mRNA records)
                genpept_id = [i for i in genpept_id if i != '-']
                if len(genpept_id) == 0:
                    hits_failed_genpept.append(h)
                    continue
                
                h.crossref_id = genpept_id
                h.crossref_method = "GenPept"
            else:
                hits_failed_genpept.append(h)
                
        return hits_failed_genpept
    
    
    def crossref_afdb_via_wgs_genpept(self, afdb_hits: list[Hit]) -> list[Hit]:
        """
        Cross-reference AFDB hits to WGS-GenPept IDs in UniSave using the UniSave API.
        
        Retrieves the WGS-GenPept IDs in UniSave for all hits and updates the the crossref ID and
        method if any.
        
        Args:
            afdb_hits: List of Hit objects without genomic information and cross-references
            
        Mutates:
            Hit objects in afdb_hits: Fills cross-reference ID and method attributes
        
        Returns:
            hits_failed_legg (list[Hit]): list of Hits that have no cross-reference to WGS-GenPept in UniSave.
            These need to be processed by a different cross-referencing method.
        """
        all_uniprot_ids = list({h.db_id for h in afdb_hits})
        hits_failed_wgs_genpept = []
        
        # Pull all UniSave records
        all_unisaves = pull_dict_from_unisave(all_uniprot_ids,
                                              max_workers = self.params['max_workers'],
                                              no_progress = self.params['no_progress'])
        
        # Fill the WGS-GenPept IDs
        for h in afdb_hits:
            unisave_record = all_unisaves[h.db_id]
            wgs_genpept_lines = [l for l in unisave_record.split('\n') 
                                 if 'DR   EMBL' in l 
                                 and 'NOT_ANNOTATED_CDS' not in l]
            
            # Skip if there is no cross-reference information
            if len(wgs_genpept_lines) == 0:
                hits_failed_wgs_genpept.append(h)
                continue
            
            # Extract what should be the WGS-GenPept ID
            wgs_genpept_id = wgs_genpept_lines[0].split(';')[2].strip()
            
            # Check if it has the right format
            if len(re.findall(r'[A-Z0-9]+\.[1-9]', wgs_genpept_id)) == 0:
                hits_failed_wgs_genpept.append(h)
                continue
            
            h.crossref_id = [wgs_genpept_id]
            h.crossref_method = "WGS-GenPept"
            
        return hits_failed_wgs_genpept
    
    
    def pull_and_parse_kegg_records(self, afdb_hits: list[Hit]) -> list[Hit]:
        """
        Pull and parse the KEGG records associated with an AFDB hit. Update the hit attributes.
        
        Pulls the KEGG records associated with an AFDB hit. It immediately extracts the scaffold and strand
        information, and the genomic coordinates from the record and updates the hit's attributes accordingly.
        
        Args:
            afdb_hits: List of Hit objects with cross-references
        
        Mutates:
            Hit objects in afdb_hits: Fills scaffold, strand and coordinates attributes
            
        Returns:
            processed_hits (list[Hit]): Hit objects with updated genomic location attributes retrieved from
            this cross-referencing method.
        """
        processed_hits = []
        pull = MultiProcessMultiplePull(n_workers = self.params['max_workers'])
        
        # KEGG Gene IDs
        all_kegg_gene_ids = list({h.crossref_id for h in afdb_hits if h.crossref_method == 'KEGG'})
        LOG.info(f'Going to pull {len(all_kegg_gene_ids)} KEGG Gene entries')
        _, gene_records = pull.pull_dict(all_kegg_gene_ids)
            
        # KEGG Genome IDs
        all_kegg_genome_ids = list({'genome:' + i.split(':')[0] for i in all_kegg_gene_ids})
        LOG.info(f'Going to pull {len(all_kegg_genome_ids)} KEGG Genome entries')
        _, genome_records = pull.pull_dict(all_kegg_genome_ids)
        genome_records = {k.split(':')[-1]: v for k,v in genome_records.items()}
        
        # Extract scaffold and coordinate information from the pulled KEGG records
        LOG.debug('Parsing downloaded KEGG entries')
        genomic_info = {k: extract_genomic_information_kegg(v) for k,v in gene_records.items()}
        scaffold_mappings = {k: extract_scaffold_mapping_kegg(v) for k,v in genome_records.items()}
        
        # Pass the extracted information on to the Hit objects
        LOG.debug('Extract CDS coordinates')
        for h in afdb_hits:
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
            
        return processed_hits
    
    
    def pull_and_parse_genpept_records(self, afdb_hits: list[Hit]) -> list[Hit]:
        """
        Pull and parse the (WGS-)GenPept records associated with an AFDB hit. Update the hit attributes.
        
        Pulls the (WGS-)GenPept records associated with an AFDB hit. It immediately extracts the scaffold and strand
        information, and the genomic coordinates from the record and updates the hit's attributes accordingly.
        
        Args:
            afdb_hits: List of Hit objects with cross-references
        
        Mutates:
            Hit objects in afdb_hits: Fills scaffold, strand and coordinates attributes
            
        Returns:
            processed_hits (list[Hit]): Hit objects with updated genomic location attributes retrieved from
            this cross-referencing method.
        """
        processed_hits = []
        
        all_genpept_gene_ids = list({h.crossref_id for h in afdb_hits if "GenPept" in h.crossref_method})
        
        # Pull GenPept ENA records
        LOG.info(f'Going to pull {len(all_genpept_gene_ids)} ENA GenPept entries')
        ena_records_pulled = thread_map(pull_from_ena, all_genpept_gene_ids,
                                        max_workers = self.params['max_workers'],
                                        leave = False,
                                        disable = self.params['no_progress'])
        ena_records = dict(zip(all_genpept_gene_ids, ena_records_pulled))
            
        # Extract scaffold and coordinate information from the pulled ENA records
        LOG.debug("Parsing downloaded GenPept entries")
        genomic_info = {k: extract_genomic_information_ena(v) for k,v in ena_records.items()}
        # Pass the extracted information on the Hit objects
        LOG.debug('Extract CDS coordinates')
        for h in afdb_hits:
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
            
        return processed_hits
    
    
    def update_version_digits(self, processed_hits: list[Hit], max_attempts: int = 3) -> list[Hit]:
        """
        Update the version digits of the scaffold IDs of every hit.
        
        Adds or updates the version digit of the scaffold ID of every hit. Retrieves the
        latest scaffold ID from the NCBI Entrez API with retry logic, and updates the hit's scaffold
        IDs accordingly.
        
        Args:
            processed_hits (list[Hit]): hits with filled cross-reference and genomic location attributes
            max_attempts (int): Maximum number of times to try getting the most recent version digits from Entrez.
        
        Mutates:
            Hit objects in processed_hits: Updates the scaffold attribute with the most recent version digit.
        
        Returns:
            processed_hits (list[Hit]): hits with an updated version digit in the scaffold attribute of every hit.            
        """
        # Get scaffold IDs with version digits using the NCBI Entrez API.
        all_old_scaffold_ids = [h.scaff for h in processed_hits]
        with warnings.catch_warnings(action = "ignore"):
            for attempt in range(max_attempts):
                try:
                    with Entrez.efetch(db = 'nucleotide',
                                       id = all_old_scaffold_ids,
                                       rettype = 'acc',
                                       retmode = 'text') as handle:
                        all_new_scaffold_ids = [l.rstrip() for l in handle.readlines()]
                    break
                except:
                    if attempt < max_attempts:
                        LOG.warning(f'Failed getting version digits from NCBI in attempt {attempt+1}. Max. attempts: {max_attempts}')
                    else:
                        msg = f'Failed getting version digits from NCBI in {max_attempts} attempts!'
                        LOG.error(msg)
                        raise RuntimeError(msg)
        
        ## Make mapping between old scaffold IDs and new ones
        version_mappings = {}
        for id_without in all_old_scaffold_ids:
            mapped_id = next((id_with for id_with in all_new_scaffold_ids if id_with.startswith(id_without)), None)
            if mapped_id:
                version_mappings[id_without] = mapped_id
        
        ## Replace the scaffold IDs
        for h in processed_hits:
            if h.scaff in version_mappings.keys():
                h.scaff = version_mappings[h.scaff]
                
        return processed_hits
    
                
    def crossref_afdb(self) -> list:
        """
        Cross-reference AFDB hits and retrieve genomic neighbourhood information.
        
        Systematically cross-references AlphaFold DB (AFDB) hits to genomic data
        through three methods: KEGG IDs from KEGG, GenPept IDs from ENA, and WGS-GenPept IDs
        from UniSave. For each hit, extracts scaffold IDs, CDS coordinates, and strand information.
        Updates scaffold IDs to their latest versions using NCBI Entrez.
        
        Returns:
            A list of Hit objects with populated genomic information (scaffold,
            coordinates, and strand data).
        """
        
        ### Identify the AFDB hits
        all_afdb_hits = [h for h in self.hits if 'afdb' in h.db]
            
        ### Get the crossreffing IDs            
        ## Method 1: try to get KEGG IDs
        LOG.info('Crossref method 1: Searching crossref IDs in KEGG')
        hits_failed_kegg = self.crossref_afdb_via_kegg(all_afdb_hits)
        LOG.info(f'{len(hits_failed_kegg)} hits have no crossref in KEGG.')
        
        ## Method 2: try to get GenPept IDs for the hits that did not get a KEGG ID.
        LOG.info('Crossref method 2: Searching crossref IDs in GenPept')
        hits_failed_genpept = self.crossref_afdb_via_genpept(hits_failed_kegg)
        LOG.info(f'..., of which {len(hits_failed_genpept)} hits have no crossref in GenPept.')
        
        ## Method 3: try to get WGS GenPept IDs for the hits that failed the previous two methods
        LOG.info('Crossref method 3: Searching crossref IDs in WGS-GenPept')
        hits_failed_wgs_genpept = self.crossref_afdb_via_wgs_genpept(hits_failed_genpept)
        LOG.info(f'..., of which {len(hits_failed_wgs_genpept)} hits have no crossref in WGS-GenPept.')
            
        ## Split records with multiple crossrefs and discard the ones without crossref
        LOG.debug("Sanitising the hit list")
        all_afdb_hits = self._sanitise_hit_attr(all_afdb_hits, 'crossref_id')
        
        ### Pulling the crossreffed records
        ## Get the scaffold and genome positions for the KEGG crossreffing method
        ## Pull all KEGG Gene and Genome records
        LOG.info("Crossref method 1: Pulling and parsing crossreffed KEGG records")
        processed_hits_kegg = self.pull_and_parse_kegg_records(all_afdb_hits)
                
        ## Get the scaffold and genome positions for the GenPept crossreffing method
        ## Pull all GenPept records from ENA in EMBL format
        LOG.info("Crossref method 2 & 3: Pulling and parsing crossreffed GenPept crossrefs from ENA")
        processed_hits_genpept = self.pull_and_parse_genpept_records(all_afdb_hits)
        
        ## Gather everything
        processed_hits = processed_hits_kegg + processed_hits_genpept
        
        ### Get the latest version digit of all accession codes
        LOG.info('Adding latest version digits to accession codes')
        processed_hits = self.update_version_digits(processed_hits)
        
        ### Final update of the hit set
        all_afdb_hits = processed_hits
        LOG.info(f'{len(all_afdb_hits)} hits have been processed.')
        
        return all_afdb_hits
    

    def run(self) -> None:
        """
        Execute the complete remote search workflow.
        
        Orchestrates all processing steps in sequence: running FoldSeek remotely via the webserver,
        parsing the results, cross-referencing using different sources, pulling the cross-referenced
        records, parsing the genomic neighbourhood coordinates from these records, and identifying 
        the gene clusters.
        
        Returns:
            None
        """
        
        LOG.info("STARTING PART 1: Executing FoldSeek search")
        self.run_foldseek()
        LOG.info('FINISHED PART 1')
        
        LOG.info("STARTING PART 2: Parsing FoldSeek results")
        self.parse_foldseek_results()
        LOG.info('FINISHED PART 2')
        
        LOG.info('STARTING PART 2B: Fetching CDS coordinates via crossreffing')
        afdb_hits = self.crossref_afdb()
        self.hits = afdb_hits
        LOG.info("FINISHED PART 2B")
        
        LOG.info('STARTING PART 3: Identifying gene clusters')
        self.identify_clusters()
        LOG.info('FINISHED PART 3')
        
        return None
        
    