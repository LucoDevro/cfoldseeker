#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import time
import requests
import polars as pl
import networkx as nx
import itertools as it
from pathlib import Path
from copy import deepcopy
from concurrent.futures import ThreadPoolExecutor
from kegg_pull.pull import MultiProcessMultiplePull
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from cblaster.classes import Session

DEFAULTS = {'mode': 'local',
            'db': ('afdb-proteome', 'afdb-swissprot'),
            'max_eval': 1,
            'min_prob': 0,
            'min_score': 0,
            'min_seqid': 25,
            'min_qcov': 70,
            'gap': 1000,
            'min_hits': 2,
            'min_cov_qrs': 2,
            'require': []
            }


class Hit:
    
    def __init__(self, uniprot, query, kegg_id = [], name = '', taxon_name = '', taxon_id = 0,
                 db = "", evalue = 1, prob = 1, score = 0, seqid = 0, qcov = 0,
                 scaff = '', coords = [], strand = ''):
        
        self.query: str = query #ID of the homologous query protein
        
        # ID attributes
        self.uniprot: str = uniprot #Uniprot ID
        self.kegg_id: list = kegg_id #KEGG Gene ID
        self.scaff: str = scaff #RefSeq ID of the scaffold encoding the hit
        
        # FoldSeek hit properties
        self.name: list = name #name of the hit in the KEGG entry
        self.taxon_name: str = taxon_name #name of the taxon in which this hit was found
        self.taxon_id: int = taxon_id
        self.db: str = db #Structure database the hit was found in
        self.evalue: float = evalue #evalue of the FoldSeek hit
        self.prob: float = prob #FoldSeek hit probability score
        self.score: int = score #FoldSeek score
        self.seqid: float = seqid #Sequence identity with the query protein
        self.qcov: float = qcov #Query coverage with the query protein
        
        # Genomic properties
        self.coords: list = coords #list of genomic coordinates of the encoding gene's exons
        self.strand: str = strand #DNA strand the encoding gene is part from
    
    def __repr__(self):
        return f"{self.query} Hit {self.uniprot}\t {self.scaff} {self.start()}-{self.end()} ({self.strand})"
    
    # Returns start coordinate of the first exon
    def start(self):
        try:
            return min(it.chain(*self.coords))
        except ValueError:
            return None
    
    # Returns end coordinate of the last exon
    def end(self):
        try:
            return max(it.chain(*self.coords))
        except ValueError:
            return None
    
    # Returns the sum of the exon lengths
    def length(self):
        return sum([abs(c[1] - c[0] + 1) for c in self.coords])
    
    # Returns the smallest distance between gene extremities, regardless of overlaps
    def distance(self, other_hit):
        all_dist = [abs(self.start() - other_hit.start()),
                    abs(self.start() - other_hit.end()),
                    abs(self.end() - other_hit.start()),
                    abs(self.end() - other_hit.end())]
        return min(all_dist)
    
    # Determines whether two hits are quasi-identical, i.e. there is at least a 90% overlap in the coordinate intervals
    def quasi_identical(self, other_hit, quasi_identity_threshold = 0.9):
        max_start = max(self.start(), other_hit.start())
        min_end = min(self.end(), other_hit.end())
        overlap_length = min_end - max_start
        x_length = abs(self.end() - self.start())
        y_length = abs(other_hit.end() - other_hit.start())
        x_overlap = overlap_length / x_length
        y_overlap = overlap_length / y_length
        if min(x_overlap, y_overlap) >= quasi_identity_threshold:
            return True
        else:
            return False
        
    # Determines whether one hit overlaps another, i.e. the first's coordinates' interval overlaps at least 90% with the second one's
    def part_of(self, other_hit, overlap_threshold = 0.9):
        max_start = max(self.start(), other_hit.start())
        min_end = min(self.end(), other_hit.end())
        overlap_length = min_end - max_start
        x_length = abs(self.end() - self.start())
        x_overlap = overlap_length / x_length
        if x_overlap >= overlap_threshold:
            return True
        else:
            return False
        

class Cluster:
    
    def __init__(self, hits, number = 0):
        
        self.hits: list(Hit) = hits
        self.number: int = number
        
        # Calculate cluster scores by summing hit scores
        self.score: int = sum([h.score for h in self.hits])
        
        # Cluster coordinates are defined as the most extreme coordinates
        self.start: int = min([h.start() for h in self.hits])
        self.end: int = max([h.end() for h in self.hits])
        
        # Cluster length is defined by the sum of hits' exons
        self.length: int = sum([h.length() for h in self.hits])
        
        # Take over cluster strand from the first hit. Throw a warning if strands do not match across the hits in the cluster
        self.strand = self.hits[0].strand
        common_strand = {h.strand for h in self.hits}
        if len(common_strand) == 0:
            print("Warning! Different coding strands found in your cluster.")
        
        # Same for scaffold ID
        self.scaff: str = self.hits[0].scaff
        common_scaff = {h.scaff for h in self.hits}
        if len(common_scaff) == 0:
            print("Warning! Different scaffolds found in your cluster.")
        
        # Same for taxon ID
        self.taxon_id: str = self.hits[0].taxon_id
        common_taxon_id = {h.taxon_id for h in self.hits}
        if len(common_taxon_id) == 0:
            print("Warning! Different taxon IDs found in your cluster.")
            
        self.taxon_name: str = self.hits[0].taxon_name
        common_taxon_name = {h.taxon_name for h in self.hits}
        if len(common_taxon_name) == 0:
            print("Warning! Different taxon names found in your cluster.")
        
        return None
    
    def __repr__(self):
        
        return f"Cluster {self.number}: {len(self.hits)} proteins from {self.scaff} ({self.start} - {self.end}), ({self.strand})\tScore: {self.score}"
    

class Search:
    
    def __init__(self, query, mapping_table_path, params = DEFAULTS,
                 hits = [], clusters = []):
        
        self.query: list = query # dictionary of query names as keys and structure filepaths as values
        self.params: dict = params # dictionary containing the search configuration
        self.hits: list = hits # list of Hit objects
        self.clusters: list = clusters # list of Cluster objects
        # The full Uniprot crossref mapping table casted in a polars DataFrame
        self.mapping_table: dict = pl.read_csv("uniprot_kegg_nucl.gz", has_header = False, separator = "\t",
                                               new_columns = ['Uniprot_ID', 'Target_DB', 'ID'])
                
        return None
    
    def __repr__(self):
        return f"Search of {','.join(list(self.query.keys()))} with {len(self.clusters)} clusters identified"
    
    def run_foldseek(self):
        """
        Submits queries to the FoldSeek webserver and collects the results.
        """

        FOLDSEEK_SUBMISSION_URL = "https://search.foldseek.com/api/ticket"
        FOLDSEEK_RESULTS_URL = "https://search.foldseek.com/api/result"
        
        """
        Submits one structure file to the FoldSeek API and returns the submission ticket in dictionary form.
        """
        def submit_foldseek_query(query_path, dbs):
            with open(query_path, "rb") as f:
                files = {"q": f}
                data = [("mode", "3diaa")]
                for db in dbs:
                    data.append(("database[]", db))
                response = requests.post(FOLDSEEK_SUBMISSION_URL, files=files, data=data)
                if response.status_code == 200:
                    return response.json()
                else:
                    raise Exception(f"Error submitting query {query_path}!")
                    
            return None
        
        """
        Checks the status of the job associated with the provided job ID. Returns the status flag.
        """
        def check_query_status(job_id):
            url = f"{FOLDSEEK_SUBMISSION_URL}/{job_id}"
            results = requests.get(url).json()
            status = results['status']
            
            return status
        
        """
        Waits for the FoldSeek job associated with the provided job ID to complete, and retrieves its result.
        Returns the parsed results in a dictionary.
        """
        def retrieve_foldseek_results(job_id):
            while True:
                status = check_query_status(job_id)
                if status == "COMPLETE":
                    entry = 0
                    url = f"{FOLDSEEK_RESULTS_URL}/{job_id}/{entry}"
                    results = requests.get(url).json()
                    break
                else:
                    time.sleep(2)
                    
            return results
            
        
        # First submit all query proteins to FoldSeek
        query_foldseek = lambda x: submit_foldseek_query(x, self.params['db'])
        with ThreadPoolExecutor(max_workers = 5) as executor:
            tickets = dict(zip(self.query.keys(), executor.map(query_foldseek, self.query.values())))
                
        # Then wait for the results and retrieve them automatically when completed
        all_job_ids = [ticket['id'] for ticket in tickets.values()]
        with ThreadPoolExecutor(max_workers = 5) as executor:
            all_results = dict(zip(self.query.keys(), executor.map(retrieve_foldseek_results, all_job_ids)))
        
        self.hits = all_results
        
        return None
    
    
    def parse_foldseek_results(self):
        """
        Parses the FoldSeek results and creates Hit objects for the passing hits
        """
        
        ## Load the thresholds from params
        max_eval = self.params['max_eval']
        min_prob = self.params['min_prob']
        min_score = self.params['min_score']
        min_seqid = self.params['min_seqid']
        min_qcov = self.params['min_qcov']
        
        ## Parse the FoldSeek results
        all_hits = []
        for query, task in self.hits.items():
            for results_by_db in task['results']:
                db = results_by_db['db']
                for hit_entry in results_by_db['alignments'][0]:
                    target = hit_entry['target'] # uniprot ID and name
                    uniprot = target.split('-')[1]
                    name = ' '.join(target.split(' ')[1:])
                    taxon_name = hit_entry['taxName'] # taxon name
                    taxon_id = hit_entry['taxId'] # taxon ID
                    evalue = float(hit_entry['eval']) # FoldSeek e-value
                    prob = float(hit_entry['prob']) # FoldSeek probability score
                    score = int(hit_entry['score']) # FoldSeek hit score
                    seqid = float(hit_entry['seqId']) # Sequence identity with the query protein
                    qcov = (int(hit_entry['qEndPos']) - int(hit_entry['qStartPos']))/int(hit_entry['qLen'])*100 # Sequence covergage
                    
                    # Create Hit object and collect it if it passes all thresholds
                    if evalue <= max_eval and prob >= min_prob and score >= min_score and seqid >= min_seqid and qcov >= min_qcov:
                        hit = Hit(uniprot, query, name = name, taxon_name = taxon_name, taxon_id = taxon_id,
                                  db = db, evalue = evalue, prob = prob, score = score, seqid = seqid, qcov = qcov)
                        all_hits.append(hit)
            
        ## Filter out redundant hits from multiple databases
        # Group by uniprot ID to find potentially redundant hits
        filtered_hits = {}
        for hit in all_hits:
            try:
                filtered_hits[hit.uniprot].append(hit)
            except KeyError:
                filtered_hits[hit.uniprot] = [hit]
                
        # Extract the redundant hits from the grouped list
        redundant_hits = [filtered_hits.pop(k) for k,v in list(filtered_hits.items()) if len(v) > 1]
        
        # Take the first hit to resolve any redundancy
        redundant_hits = [[v[0]] for v in redundant_hits]
        
        # Put everything back together and sort by Uniprot ID
        filtered_hits = list(filtered_hits.values()) + redundant_hits
        filtered_hits = sorted(list(it.chain(*filtered_hits)), key = lambda x: x.uniprot)
        
        self.hits = filtered_hits
        
        return None
            
                
    def crossref(self):
        """
        Crossrefs all hits where possible and retrieves the genomic neighbourhood information (scaffold IDs and coordinates)
        """
        def extract_positional_information(entry: str) -> str:
            """
            Extracts the genomic positional information from a full pulled KEGG record.
            """
            try:
                return [line for line in entry.split('\n') if 'POSITION' in line][0][12:]
            except IndexError:
                return None
        
        def prepare_mapping_dict(crossref_df: pl.DataFrame, target_db: str) -> dict:
            """
            Extracts a mapping dictionary from the full Uniprot crossref mapping table, linking the Uniprot ID (keys) with a list of
            one or more multiple crossrefs in the provided target database (values).
            """
            res = crossref_df.filter(crossref_df['Target_DB'] == target_db).drop('Target_DB')
            res = res.group_by('Uniprot_ID').all()
            res = dict(zip(res['Uniprot_ID'], res['ID'].to_list()))
            
            return res
        
        def sanitise_attr(self, attr):
            """
            Auxiliary function that standardises the hits. In case of multiple values for a given attribute, it splits it into 
            multiple hits that are identical except for that exact crossref attribute. The split hits are first introduced instead of the old hit before the full
            hit list is being flattened.
            In case of an absent or default attribute value, the hit is discarded.
            
            Mutates:
                self.hits: list: the hit list containing Hit objects that represent proteins that are structurally homologous to
                                 one of the query proteins.
            """
            sanitised_hits = []
            # Loop over all hits, but keep track of the index
            for h in self.hits:
                # Determine whether we should split a record
                nb_copies = len(getattr(h, attr))
                # If we need to split a record, make as much copies as there are crossrefs and adjust each crossref to one of these
                if nb_copies > 1:
                    split_records = []
                    copy_counter = 0
                    while copy_counter < nb_copies:
                        new_record = deepcopy(h)
                        setattr(new_record, attr, getattr(h, attr)[copy_counter])
                        split_records.append(new_record)
                        copy_counter += 1
                    # Insert at the index of the old record
                    sanitised_hits.append(split_records)
                # If we don't need to split the record, just unlist it
                elif nb_copies == 1:
                    setattr(h, attr, getattr(h, attr)[0])
                    # Insert at the index of the old record
                    sanitised_hits.append([h])
                # If it's empty, don't include it in the new hit set
                else:
                    continue
            
            # Flatten the full hit list again
            self.hits = list(it.chain(*sanitised_hits))
            
            return None
        
        
        ## First, fill all KEGG IDs
        # Extract a mapping table for KEGG IDs from the Uniprot crossref mapping table
        all_uniprot_ids = [h.uniprot for h in self.hits]
        all_cross_refs = self.mapping_table.filter(self.mapping_table['Uniprot_ID'].is_in(all_uniprot_ids))
        all_uniprot_kegg = prepare_mapping_dict(all_cross_refs, 'KEGG')
        # Fill the KEGG IDs, if possible
        for h in self.hits:
            try:
                h.kegg_id = all_uniprot_kegg[h.uniprot]
            except KeyError:
                continue
        # Split records with multiple crossrefs and discard the ones without crossref
        sanitise_attr(self, 'kegg_id')
            
        ## Then, pull all KEGG records
        all_kegg_ids = [h.kegg_id for h in self.hits]
        pull = MultiProcessMultiplePull(n_workers = 2)
        _, records = pull.pull_dict(all_kegg_ids)
        
        ## Extract genomic information from the pulled KEGG entries
        positions = {k: extract_positional_information(v) for k,v in records.items()}
        # Extract a mapping table for RefSeq Nucleotide IDs from the Uniprot crossref mapping table
        all_uniprot_refseq = prepare_mapping_dict(all_cross_refs, 'RefSeq_NT')
        # Start sorting out the genomic information
        processed_hits = []
        for h in self.hits:
            position_info = positions[h.kegg_id]
            # If there is not even position info, there is no point of keeping the hit
            if position_info == None:
                continue
            # If there are no numeric characters, there is no positional information
            elif not(any(char.isdigit() for char in position_info)):
                continue
            # If there is no colon, then there probably is no scaffold information,
            # but we have a second data source in our Uniprot crossref mapping table
            elif ':' not in position_info:
                try:
                    refseq_id = all_uniprot_refseq[h.uniprot] # This might result in a double crossref, which will be fixed later on
                    coords = position_info
                # Bad luck, there's no RefSeq crossref in Uniprot either, discard the hit
                except KeyError:
                    continue
            # Check whether it has the format of a complete genomic position field as we would like to have it from KEGG
            elif len(re.findall(r'.+:.*<?[0-9]+\.\.>?[0-9]+.*', position_info)) > 0:
                refseq_id, coords = position_info.split(':')
                refseq_id = [refseq_id] # wrapped in a list for consistency with the output of the case above
            # Something else not encountered so far
            else:
                print(f'Invalid positional record found for {h.kegg_id}! Skipping...')
                continue
                        
            # Save the scaffold and the genomic coordinates
            h.scaff = refseq_id

            # Parse the genomic coordinates. The fact that we got here in the loop implies that there is parsable information.
            coords = coords.translate(str.maketrans('', '', '<>'))
            coord_groups = re.findall(r'\d+\.\.\d+', coords)
            coord_groups = [i.split('..') for i in coord_groups]
            coord_groups = [[int(j) for j in i] for i in coord_groups]
            h.coords = coord_groups
            
            # Need to do some further processing in case of a reverse strand gene
            if 'complement' in coords:
                h.strand = "-"
            else:
                h.strand = '+'
                
            # Stack the processed hits
            processed_hits.append(h)
                
        ## Update the hit set and sanitise the scaff attribute
        self.hits = processed_hits
        sanitise_attr(self, 'scaff')
        
        return None
        
    
    def identify_clusters(self):
        """
        Identifies the gene clusters among the hits.
        """
        
        ### Load the requirements from params
        gap = self.params['gap']
        min_hits = self.params['min_hits']
        min_covered_queries = self.params['min_cov_qrs']
        require = self.params['require']
        
        ### Cluster identification
        ## First, make groups by scaffold and taxon ID
        grouped_scaffs = {}
        for h in self.hits:
            try:
                grouped_scaffs[(h.scaff, h.taxon_id)].append(h)
            except KeyError:
                grouped_scaffs[(h.scaff, h.taxon_id)] = [h]
        
        ## Then, calculate the intergenic distance between all connections and filter out the ones failing the intergenic threshold
        close_groups = []
        for group, hits in grouped_scaffs.items():
            dists = {pair: Hit.distance(*pair) for pair in it.combinations(hits, 2)}
            dists = {k:v for k,v in dists.items() if v <= gap}
            close_groups.append(list(dists.keys()))
        
        ## Then identify the clusters by identifying series of connected Hit objects using an undirected network graph
        clusters = []
        for cg in close_groups:
            if len(cg) == 0:
                continue
            G = nx.Graph()
            G.add_edges_from(cg)
            cluster = list(nx.connected_components(G))
            clusters.append(cluster)
        clusters = list(it.chain(*clusters))
        
        ### Reducing redundancy introduced by crossreffing
        ## Collapse quasi-identical hits within clusters
        collapsed_clusters = []
        for cluster in clusters:
            # First identify the quasi_identical hit pairs
            quasi_identical = {pair: Hit.quasi_identical(*pair) for pair in it.combinations(cluster, 2)}
            quasi_identical = [k for k,v in quasi_identical.items() if v]
            if len(quasi_identical) == 0:
                collapsed_clusters.append(cluster)
                continue
            # Then identify groups of quasi-identical pairs using an undirected network graph
            G = nx.Graph()
            G.add_edges_from(quasi_identical)
            quasi_identical_groups = list(nx.connected_components(G))
            # Collapse the groups by taking the first item of each group
            collapsed_groups = [{list(qig)[0]} for qig in quasi_identical_groups]
            hits_to_remove = set.union(*quasi_identical_groups) # the complete groups
            hits_to_readd = set.union(*collapsed_groups) # the selected singletons
            collapsed_clusters.append(cluster - hits_to_remove | hits_to_readd)
        
        ## Collapse partial overlap hits within clusters (by keeping the smaller parts and discarding apparent fusion proteins)
        collapsed_clusters_no_fusions = []
        for cluster in collapsed_clusters:
            if len(cluster) == 1:
                collapsed_clusters_no_fusions.append(cluster)
                continue
            # First identify hits that are part of another hit
            master_and_parts = {pair: Hit.part_of(*pair) for pair in it.permutations(cluster, 2)}
            master_and_parts = [k for k,v in master_and_parts.items() if v]
            if len(master_and_parts) == 0:
                collapsed_clusters_no_fusions.append(cluster)
                continue
            # Then identify the master (or sink) hits (i.e. the fusion hits) using a directed network graph
            G = nx.DiGraph()
            G.add_edges_from(master_and_parts)
            masters = {node for node in G.nodes() if G.out_degree(node) == 0}
            # Remove these master hits from the cluster
            collapsed_clusters_no_fusions.append(cluster - masters)
            
        ### Applying cluster requirements
        ## Select for requested cluster configurations
        # Minimum number of hits in a cluster
        resulting_set = [cl for cl in collapsed_clusters_no_fusions if len(cl) >= min_hits]
        # Minimum number of covered queries and required queries
        covered_queries = [{h.query for h in cl} for cl in resulting_set]
        resulting_set = [cl for cl,qrs in zip(resulting_set, covered_queries)
                         if len(qrs) >= min_covered_queries and set(require) <= qrs]
        
        ### Preparing results
        ## Reorder the clusters by genomic coordinates
        results = [sorted(cl, key = lambda h: h.start()) for cl in resulting_set]
        
        ## Create the Cluster objects from the processed hit clusters
        res_objects = [Cluster(cl, number = idx) for idx,cl in enumerate(results)]
        
        ## Rank by cluster score and renumber
        res_objects.sort(key = lambda x: x.score, reverse = True)
        for idx,cl in enumerate(res_objects):
            cl.number = idx+1
        
        self.clusters = res_objects
        
        return None
   
    
    def generate_cblaster_session(self):
        """
        Generates a cblaster session object
        """
        def get_sequence_length_from_cif(file):
            """
            Determines the sequence length of the structure in the input CIF file.
            """
            structure = MMCIF2Dict(file)
            res_ids = [int(i) for i in structure['_entity_poly_seq.num']]
            return max(res_ids)
        
        def get_clusters_by_id(self, nbs):
            """
            Returns the cluster object that has the given cluster number.
            """
            return [cl for cl in self.clusters if cl.number in nbs]
        
        def get_taxon_name_from_taxon_id(self, txid):
            """
            Returns the taxon name associated with a given taxon ID.
            """
            return [cl.taxon_name for cl in self.clusters if cl.taxon_id == txid][0]
        
        session_dict = {}
        
        ### Queries field
        session_dict['queries'] = list(self.query.keys())
        
        ### Query field
        cblaster_query = {}
        cblaster_query['indices'] = []
        ## Query Subjects field
        cblaster_query_subjects = []
        previous_end = -500
        for qry,file in self.query.items():
            this_subject = {}
            this_subject['id'] = None
            this_subject['hits'] = []
            this_subject['name'] = qry
            this_subject['ipg'] = None
            
            # Sequence length is determined from the number of modelled residues in the CIF file
            # cblaster session files take a margin of 500 aa between the hits
            length = get_sequence_length_from_cif(file)
            this_subject['start'] = previous_end + 500
            this_subject['end'] = this_subject['start'] + length
            previous_end = this_subject['end']
            
            this_subject['strand'] = 1
            this_subject['sequence'] = ""
            
            cblaster_query_subjects.append(this_subject)
            
        cblaster_query['subjects'] = cblaster_query_subjects
        cblaster_query['intermediate_genes'] = []
        cblaster_query['score'] = len(self.query.keys())
        cblaster_query['start'] = 0
        cblaster_query['end'] = cblaster_query_subjects[-1]['end']
        cblaster_query['number'] = 0
        
        session_dict['query'] = cblaster_query
        
        ### Params field
        cblaster_params = {}
        cblaster_params['mode'] = self.params['mode']
        cblaster_params['database'] = list(self.params['db'])
        cblaster_params['min_identity'] = self.params['min_seqid']
        cblaster_params['min_coverage'] = self.params['min_qcov']
        cblaster_params['max_evalue'] = self.params['max_eval']
        cblaster_params['require'] = self.params['require']
        cblaster_params['query_file'] = None
        cblaster_params['rid'] = None
        cblaster_params['entrez_query'] = ""
        
        session_dict['params'] = cblaster_params
        
        ### Organisms field
        cblaster_organisms = {}
        
        ## Group the session's cluster records in the same way as they will be structured in the session file,
        ## which will make processing much easier
        
        # Prepare a polars dataframe of taxon IDs (organisms), scaffold IDs and cluster numbers (clusters)
        taxon_ids = [cl.taxon_id for cl in self.clusters]
        scaffs = [cl.scaff for cl in self.clusters]
        cl_nbs = [cl.number for cl in self.clusters]
        cl_df = pl.DataFrame({'taxon_ids': taxon_ids, 'scaffolds': scaffs, 'cluster_number': cl_nbs})
        
        # Group by taxon ID and scaffold ID, and cast into a dictionary
        grouped_cl_df = cl_df.group_by(['taxon_ids', 'scaffolds']).all()
        grouped_cl_dict = {(row[0], row[1]): get_clusters_by_id(self, row[2]) for row in grouped_cl_df.iter_rows()}
        
        ## Make the cblaster session fields inside out, i.e. populate organisms first with scaffolds (and other attributes),
        ## then populate the scaffolds with clusters and subjects, then populate the clusters with links to the subjects.
        for (txid, scaff), clusters in grouped_cl_dict.items():
            # Create a new organism instance of there's no one for this taxon ID
            if txid in cblaster_organisms.keys():
                this_organism = cblaster_organisms[txid]
            else:
                this_organism = {'name': get_taxon_name_from_taxon_id(self, txid),
                                 'strain': f"(txid: {txid})",
                                 'scaffolds': {}}
                cblaster_organisms[txid] = this_organism
            
            # Create a new scaffold instance of there is no one for this scaffold ID in this taxon
            if scaff in this_organism['scaffolds'].keys():
                this_scaffold = this_organism['scaffolds'][scaff]
            else:
                this_scaffold = {'accession': scaff,
                                 'subjects': [],
                                 'clusters': []}
                this_organism['scaffolds'][scaff] = this_scaffold
                
            # Populate this scaffold with cluster for each cluster identified on this scaffold
            for cl in clusters:
                cblaster_this_cluster = {}
                cblaster_this_cluster['indices'] = list(range(len(cl.hits)))
                cblaster_this_cluster['subjects'] = []
                cblaster_this_cluster['intermediate_genes'] = []
                cblaster_this_cluster['score'] = cl.score
                cblaster_this_cluster['start'] = cl.start
                cblaster_this_cluster['end'] = cl.end
                cblaster_this_cluster['number'] = cl.number
                this_scaffold['clusters'].append(cblaster_this_cluster)
                
                # Populate this scaffold with subjects (hits) for each hit identified on this scaffold
                # and add the link with the cluster
                these_hits = cl.hits
                for hit in these_hits:
                    cblaster_this_subject = {}
                    cblaster_this_subject['id'] = None
                    cblaster_this_subject['hits'] = [{'query': hit.query,
                                                      'subject': hit.uniprot,
                                                      'identity': hit.seqid,
                                                      'coverage': hit.qcov,
                                                      'evalue': hit.evalue,
                                                      'bitscore': hit.score}]
                    cblaster_this_subject['name'] = hit.uniprot
                    cblaster_this_subject['ipg'] = hit.kegg_id
                    cblaster_this_subject['start'] = hit.start()
                    cblaster_this_subject['end'] = hit.end()
                    cblaster_this_subject['strand'] = int(f"{hit.strand}1")
                    cblaster_this_subject['sequence'] = None
                    this_scaffold['subjects'].append(cblaster_this_subject)
        
        ## Discard the taxon ID and scaffold ID indexing to get the cblaster session formatting
        cblaster_organisms = list(cblaster_organisms.values())
        cblaster_organisms_new = []
        for organism in cblaster_organisms:
            organism['scaffolds'] = list(organism['scaffolds'].values())
            cblaster_organisms_new.append(organism)
            
        session_dict['organisms'] = cblaster_organisms_new
        
        ### Construct the cblaster Session
        session = Session.from_dict(session_dict)
        
        return session
        
    
# Queries
query = {'query1': '/home/lucas/bin/cfoldseeker_old/AF3_models/fold_sco5088_model_0.cif',
         'query2': '/home/lucas/bin/cfoldseeker_old/AF3_models/fold_sco5089_model_0.cif'
         }
    
# Search parameters
params = {'mode': 'local',
          'db': ('afdb-proteome', 'afdb-swissprot'),
          'max_eval': 1,
          'min_prob': 0,
          'min_score': 0,
          'min_seqid': 25,
          'min_qcov': 70,
          'gap': 1000,
          'min_hits': 2,
          'min_cov_qrs': 2,
          'require': []
        }
    
s = Search(query, "uniprot_kegg_nucl.gz", params = params)
s.run_foldseek()
s.parse_foldseek_results()
s.crossref()
s.identify_clusters()
session = s.generate_cblaster_session()
