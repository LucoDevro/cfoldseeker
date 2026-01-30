#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import time
import sys
import operator
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
            'max_length': 20000,
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
        self.mapping_table: dict = pl.read_csv("uniprot_kegg.gz", has_header = False, separator = "\t",
                                               new_columns = ['Uniprot_ID', 'KEGG_ID'])
                
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
                        
        if len(all_hits) == 0:
            print("No hits identified by FoldSeek!")
            sys.exit()
            
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
        filtered_hits = sorted(list(it.chain(*filtered_hits)), key = operator.attrgetter('uniprot'))
        
        self.hits = filtered_hits
        
        return None
            
                
    def crossref(self):
        """
        Crossrefs all hits where possible and retrieves the genomic neighbourhood information (scaffold IDs and coordinates)
        """
        def extract_positional_information(gene_entry: str) -> dict:
            """
            Extracts the genomic positional information from a pulled KEGG Gene record.
            """
            position_line = [line for line in gene_entry.split('\n') if 'POSITION' in line]
            position_info = {}
            if len(position_line) == 0:
                return position_info
            else:
                position_line = position_line[0]
                if ':' in position_line:
                    internal_scaffold_id, coords = position_line.split(':')
                    internal_scaffold_id = internal_scaffold_id[12:]
                else:
                    internal_scaffold_id, coords = '', position_line
                coords = coords.translate(str.maketrans('', '', '<>'))
                coord_groups = re.findall(r'\d+\.\.\d+', coords)
                if len(coord_groups) == 0:
                    return position_info
                coord_groups = [i.split('..') for i in coord_groups]
                coord_groups = [[int(j) for j in i] for i in coord_groups]
                if 'completement' in coords:
                    strand = "-"
                else:
                    strand = "+"
                
                position_info['scaffold'] = internal_scaffold_id
                position_info['coords'] = coord_groups
                position_info['strand'] = strand
                
                return position_info

        def extract_scaffold_information(genome_entry: str) -> dict:
            """
            Extracts the chromosome information from a pulled KEGG Genome record.
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
                        
        def prepare_mapping_dict(crossref_df: pl.DataFrame) -> dict:
            """
            Extracts a mapping dictionary from the full Uniprot crossref mapping table, linking the Uniprot ID (keys) with a list of
            one or more multiple crossrefs in the provided target database (values).
            """
            res = crossref_df.group_by('Uniprot_ID').all()
            res = dict(zip(res['Uniprot_ID'], res['KEGG_ID'].to_list()))
            
            return res
        
        def sanitise_hit_attr(hits, attr) -> list:
            """
            Auxiliary function that sanitises the hit list. In case of multiple values for a given attribute, it splits it into 
            multiple hit instances that are identical except for that exact crossref attribute. The split hits are first introduced
            in a list replacing the old hit instance before the full hit list is being flattened.
            In case of an absent attribute, the hit is discarded.
            """
            sanitised_hits = []
            # Loop over all hits, but keep track of the index
            for h in hits:
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
            hits = list(it.chain(*sanitised_hits))
            
            return hits
        
        
        ## First, fill all KEGG IDs
        # Extract a mapping table for KEGG IDs from the Uniprot crossref mapping table
        all_uniprot_ids = [h.uniprot for h in self.hits]
        all_cross_refs = self.mapping_table.filter(self.mapping_table['Uniprot_ID'].is_in(all_uniprot_ids))
        all_uniprot_kegg = prepare_mapping_dict(all_cross_refs)
        # Fill the KEGG IDs, if possible
        for h in self.hits:
            if h.uniprot in all_uniprot_kegg.keys():
                h.kegg_id = all_uniprot_kegg[h.uniprot]
            else:
                continue
        # Split records with multiple crossrefs and discard the ones without crossref
        self.hits = sanitise_hit_attr(self.hits, 'kegg_id')
            
        ## Pull all KEGG Gene and Genome records
        pull = MultiProcessMultiplePull(n_workers = 10)
        # KEGG Gene IDs
        all_kegg_gene_ids = [h.kegg_id for h in self.hits]
        _, gene_records = pull.pull_dict(all_kegg_gene_ids)
        # KEGG Genome IDs
        all_kegg_genome_ids = list({'genome:' + i.split(':')[0] for i in all_kegg_gene_ids})
        _, genome_records = pull.pull_dict(all_kegg_genome_ids)
        genome_records = {k.split(':')[-1]: v for k,v in genome_records.items()}
        
        ## Extract scaffold and coordinate information
        positions = {k: extract_positional_information(v) for k,v in gene_records.items()}
        scaffolds = {k: extract_scaffold_information(v) for k,v in genome_records.items()}
        processed_hits = []
        for h in self.hits:
            try:
                position_info = positions[h.kegg_id]
                scaffold_info = scaffolds[h.kegg_id.split(':')[0]]
            except KeyError:
                continue
            if len(position_info) == 0 or len(scaffold_info) == 0:
                continue
            h.coords = position_info['coords']
            h.strand = position_info['strand']
            h.scaff = scaffold_info[position_info['scaffold']]
            
            processed_hits.append(h)
                
        ## Update the hit set and sanitise the scaff attribute
        self.hits = processed_hits
        
        return None
        
    
    def identify_clusters(self):
        """
        Identifies the gene clusters among the hits.
        """
        
        ### Load the requirements from params
        gap = self.params['gap']
        max_length = self.params['max_length']
        min_hits = self.params['min_hits']
        min_covered_queries = self.params['min_cov_qrs']
        require = self.params['require']
        
        ### Cluster identification
        ## First, make groups by scaffold and taxon ID
        grouped_scaffs = {}
        for h in self.hits:
            if (h.scaff, h.taxon_id) in grouped_scaffs.keys():
                grouped_scaffs[(h.scaff, h.taxon_id)].append(h)
            else:
                grouped_scaffs[(h.scaff, h.taxon_id)] = [h]
        
        ## Then, calculate the intergenic distance between all connections and filter out the ones failing the intergenic threshold
        close_groups = []
        for group, hits in grouped_scaffs.items():
            dists = {pair: Hit.distance(*pair) for pair in it.combinations(hits, 2)}
            dists = {k:v for k,v in dists.items() if v <= gap}
            close_groups.append(list(dists.keys()))
        if len(close_groups) == 0:
            print("No cluster could be identified!")
            sys.exit()
        
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
        
        ### Apply intra-cluster requirements
        # Minimum number of hits in a cluster
        clusters_filt = [cl for cl in clusters if len(cl) >= min_hits]
        # Minimum number of covered queries and required queries
        covered_queries = [{h.query for h in cl} for cl in clusters_filt]
        clusters_filt = [cl for cl,qrs in zip(clusters_filt, covered_queries)
                        if len(qrs) >= min_covered_queries and set(require) <= qrs]
        
        ### Create the Cluster objects from the filtered hit clusters
        res_objects = [Cluster(list(cl), number = idx) for idx,cl in enumerate(clusters_filt)]
        
        ### Filter for maximum cluster length
        res_objects_filt = [cl for cl in res_objects if cl.length <= max_length]
        
        ## Rank by cluster score and renumber
        res_objects_filt.sort(key = operator.attrgetter('score'), reverse = True)
        for idx,cl in enumerate(res_objects_filt):
            cl.number = idx+1
        
        self.clusters = res_objects_filt
        
        return None
   
    
    def generate_cblaster_session(self):
        """
        Generates a cblaster session object
        """
        def get_sequence_length_from_cif(file):
            """
            Determines the CDS sequence length of the protein structure in the input CIF file.
            """
            structure = MMCIF2Dict(file)
            res_ids = [int(i) for i in structure['_entity_poly_seq.num']]
            return max(res_ids)*3
        
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
            # Keep track of the number of hits covered on this scaffold for proper references in output files
            nb_hits_covered = 0
            for cl in clusters:
                cblaster_this_cluster = {}
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
                hits_covered_this_cluster = len(these_hits)
                cblaster_this_cluster['indices'] = list(range(nb_hits_covered, nb_hits_covered + hits_covered_this_cluster))
                nb_hits_covered += hits_covered_this_cluster
                for hit in these_hits:
                    cblaster_this_subject = {}
                    cblaster_this_subject['id'] = None
                    cblaster_this_subject['hits'] = [{'query': hit.query,
                                                      'subject': hit.kegg_id,
                                                      'identity': hit.seqid,
                                                      'coverage': hit.qcov,
                                                      'evalue': hit.evalue,
                                                      'bitscore': hit.score}]
                    cblaster_this_subject['name'] = hit.kegg_id
                    cblaster_this_subject['ipg'] = hit.uniprot
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
          'max_length': 20000,
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

with open("test_session.json", "w") as handle:
    session.to_json(fp = handle)
    
with open('test_summary', 'w') as handle:
    session.format(form = "summary", fp = handle)
    
with open('test_binary', 'w') as handle:
    session.format(form = "binary", fp = handle)
