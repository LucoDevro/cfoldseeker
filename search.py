#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from classes import Hit, Cluster

import sys
import operator
import polars as pl
import networkx as nx
import itertools as it
from abc import ABC
from copy import deepcopy
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from cblaster.classes import Session


DEFAULTS = {'mode': 'local',
            'db': ('afdb-proteome', 'afdb-swissprot'),
            'max_eval': 1,
            'min_prob': 0,
            'min_score': 0,
            'min_seqid': 25,
            'min_qcov': 70,
            'min_tcov': 70,
            'max_gap': 1000,
            'max_length': 20000,
            'min_hits': 2,
            'min_cov_qrs': 2,
            'require': []
            }


def _sanitise_hit_attr(hits: list, attr: str) -> list:
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


class Search(ABC):
    
    def __init__(self, query, mapping_table_path, params = DEFAULTS,
                  hits = [], clusters = []):
        
        self.query: list = query # dictionary of query names as keys and structure filepaths as values
        self.params: dict = params # dictionary containing the search configuration
        self.hits: list = hits # list of Hit objects
        self.clusters: list = clusters # list of Cluster objects
                
        return None
    
    def __repr__(self):
        return f"Search of {','.join(list(self.query.keys()))} with {len(self.clusters)} clusters identified"
    
    
    def identify_clusters(self):
        """
        Identifies the gene clusters among the hits.
        """
        
        ### Load the requirements from params
        max_gap = self.params['max_gap']
        max_length = self.params['max_length']
        min_hits = self.params['min_hits']
        min_covered_queries = self.params['min_cov_qrs']
        require = self.params['require']
        
        ### Cluster identification
        ## First, make groups by scaffold
        scaff_groups = {}
        for h in self.hits:
            if h.scaff in scaff_groups.keys():
                scaff_groups[h.scaff].append(h)
            else:
                scaff_groups[h.scaff] = [h]
        
        ## Then, calculate the intergenic distance between all connections on the same scaffold
        ## and filter out the ones failing the intergenic threshold
        ## and filter out self-hits as these are not genuine collocalised genes
        close_groups = []
        for _, hits in scaff_groups.items():
            # Calculate the intergenetic distances and find the self-hits
            pairs_to_test = list(it.combinations(hits, 2))
            self_hits = {pair: Hit.same_spot(*pair) for pair in pairs_to_test}
            dists = {pair: Hit.intergenic_distance(*pair) for pair in pairs_to_test}
            
            # Apply the filtering
            dists = {k:v for k,v in dists.items() if v <= max_gap} # apply max gap criterium
            dists = {k:v for k,v in dists.items() if not self_hits[k]} # filter out self-hits
            
            # Collect
            if len(dists) > 0:
                close_groups.append(list(dists.keys()))
        if len(close_groups) == 0:
            print("No cluster could be identified!")
            sys.exit()
        
        ## Then identify the clusters by chained hits on the same scaffold using a directed network graph
        ## Account for multi-hits and -crossrefs by generating all possible hit chains when encountering hits on the same genomic spot
        clusters = []
        for cg in close_groups:
            # Order every hit pair so from up- to downstream
            reordered_cg = [sorted(pair, key = operator.methodcaller('start')) for pair in cg]
            
            # Identify the hit chains
            G = nx.DiGraph()
            G.add_edges_from(reordered_cg)
            chains = list(nx.weakly_connected_components(G))
            
            # Then, identify all possible chains by generating chains for all multi-hit or -crossref combinations
            all_paths = []
            for chain in chains:
                subG = G.subgraph(chain)
                
                # Identify all possible hits to start a chain and to end a chain
                min_start = min([h.start() for h in chain])
                max_start = max([h.start() for h in chain])
                firsts = [h for h in chain if h.start() == min_start]
                lasts = [h for h in chain if h.start() == max_start]
                
                # Generate all possible hit chains
                all_paths_this_chain = [list(nx.all_simple_paths(subG, first, last)) for first in firsts for last in lasts]
                all_paths_this_chain = list(it.chain(*all_paths_this_chain))
                # Filter out paths containing shortcuts, i.e. keep only longest paths
                max_path_length = max([len(p) for p in all_paths_this_chain])
                all_paths_this_chain = [p for p in all_paths_this_chain if len(p) == max_path_length]
                all_paths.append(all_paths_this_chain)
            
            # Save the chains for this neighbour group
            all_paths = list(it.chain(*all_paths))
            clusters.append(all_paths)
        
        # Flatten out all results
        clusters = list(it.chain(*clusters))
        
        ### Apply intra-cluster filtering requirements
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
                                                      'subject': hit.crossref_id,
                                                      'identity': hit.seqid,
                                                      'coverage': hit.qcov,
                                                      'evalue': hit.evalue,
                                                      'bitscore': hit.score}]
                    cblaster_this_subject['name'] = hit.crossref_id
                    cblaster_this_subject['ipg'] = hit.db_id
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
    
    