#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import operator
import itertools as it
import sys
import logging
import polars as pl
import networkx as nx
from abc import ABC, abstractmethod
from copy import deepcopy
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from cblaster.classes import Session


LOG = logging.getLogger(__name__)


class Hit:
    
    def __init__(self, db_id, query, crossref_id = [], crossref_method = '', name = '', 
                 taxon_name = '', taxon_id = 0, db = "", evalue = 1, score = 0, 
                 seqid = 0, qcov = 0, tcov = 0, scaff = '', coords = [], strand = ''):
        
        self.query: str = query #ID of the homologous query protein
        
        # ID attributes
        self.db_id: str = db_id #ID of the hit in its DB
        self.db: str = db #Structure database the hit was found in
        self.crossref_id: list = crossref_id #ID used for crossreffing (either ID from KEGG or GenPept)
        self.crossref_method: str = crossref_method #Method used for crossreffing (either KEGG or GenPept)
        
        # FoldSeek hit properties
        self.name: list = name #annotation
        self.taxon_name: str = taxon_name #name of the taxon in which this hit was found
        self.taxon_id: int = taxon_id
        self.evalue: float = evalue #evalue of the FoldSeek hit
        self.score: int = score #FoldSeek score
        self.seqid: float = seqid #Sequence identity with the query protein
        self.qcov: float = qcov #Query coverage
        self.tcov: float = tcov #Target coverage
        
        # Genomic properties
        self.scaff: str = scaff #RefSeq or GenBank ID of the scaffold encoding the hit
        self.coords: list = coords #list of genomic coordinates of the encoding gene's exons
        self.strand: str = strand #DNA strand the encoding gene is part from
    
    def __repr__(self) -> str:
        return f"{self.query} Hit {self.db_id}\t {self.scaff} {self.start()}-{self.end()} ({self.strand})"
    
    def as_dict(self) -> dict:
        return {'query': self.query,
                'db_id': self.db_id,
                'db': self.db,
                'crossref_id': self.crossref_id,
                'crossref_method': self.crossref_method,
                'name': self.name,
                'taxon_name': self.taxon_name,
                'taxon_id': self.taxon_id,
                'evalue': self.evalue,
                'score': self.score,
                'seqid': self.seqid,
                'qcov': self.qcov,
                'tcov': self.tcov,
                'scaff': self.scaff,
                'coords': ','.join(['..'.join([str(i) for i in seq]) for seq in self.coords]),
                'strand': self.strand}
    
    # Returns start coordinate of the first exon
    def start(self) -> int | None:
        try:
            return min(it.chain(*self.coords))
        except ValueError:
            return None
    
    # Returns end coordinate of the last exon
    def end(self) -> int | None:
        try:
            return max(it.chain(*self.coords))
        except ValueError:
            return None
    
    # Returns the sum of the exon lengths
    def length(self) -> int:
        return sum([abs(c[1] - c[0] + 1) for c in self.coords])
    
    # Returns the intergenic distance between two genes. Negative if they overlap
    def intergenic_distance(self, other_hit: 'Hit') -> int:
        first = min([self, other_hit], key = operator.methodcaller('start'))
        first_start = first.start()
        first_end = first.end()
        last = max([self, other_hit], key = operator.methodcaller('start'))
        last_start = last.start()
        last_end = last.end()
        
        # In case of full overlap, return the negative of the length of the smaller gene
        if last_start >= first_start and last_end <= first_end:
            return -last.length()
        # In case of at most a partial overlap, return the intergenic distance.
        else:
            return last_start - first_end 
        
    # Checks whether two hits are at exactly the same genomic coordinates
    def same_location(self, other_hit: 'Hit') -> bool:
        first = min([self, other_hit], key = operator.methodcaller('start'))
        last = max([self, other_hit], key = operator.methodcaller('start'))
        return last.start() >= first.start() and last.end() <= first.end() and first.scaff == last.scaff
    

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
        self.strand: str = self.hits[0].strand
        common_strand = {h.strand for h in self.hits}
        if len(common_strand) == 0:
            LOG.warning(f'Different coding strands found among the gene hits in cluster {number}!')
        
        # Same for scaffold ID
        self.scaff: str = self.hits[0].scaff
        common_scaff = {h.scaff for h in self.hits}
        if len(common_scaff) == 0:
            LOG.warning(f"Different scaffolds found among the gene hits in cluster {number}!")
        
        # Same for taxon ID
        self.taxon_id: str = self.hits[0].taxon_id
        common_taxon_id = {h.taxon_id for h in self.hits}
        if len(common_taxon_id) == 0:
            LOG.warning(f"Different taxon IDs found among the gene hits in cluster {number}.")
            
        self.taxon_name: str = self.hits[0].taxon_name
        common_taxon_name = {h.taxon_name for h in self.hits}
        if len(common_taxon_name) == 0:
            LOG.warning(f"ifferent taxon names found ammong the gene hits in cluster {number}.")
        
        return None
    
    def __repr__(self) -> str:
        return f"Cluster {self.number}: {len(self.hits)} proteins from {self.scaff} ({self.start} - {self.end}), ({self.strand})\tScore: {self.score}"
    
    def as_dict(self) -> dict:
        return {'hits': ','.join([h.db_id for h in self.hits]),
                'number': self.number,
                'score': self.score,
                'start': self.start,
                'end': self.end,
                'length': self.length,
                'strand': self.strand,
                'scaff': self.scaff,
                'taxon_id': self.taxon_id,
                'taxon_name': self.taxon_name}
    

def _sanitise_hit_attr(hits: list, attr: str) -> list:
    """
    Auxiliary function that sanitises the hit list. In case of multiple values for a given attribute, it splits it into 
    multiple hit instances that are identical except for that exact crossref attribute. The split hits are first introduced
    in a list replacing the old hit instance before the full hit list is being flattened.
    In case of an absent attribute, the hit is discarded.
    """
    sanitised_hits = []
    split_counter = 0
    discard_counter = 0
    untouched_counter = 0
    # Loop over all hits, but keep track of the index
    for h in hits:
        # Determine whether we should split a record
        nb_copies = len(getattr(h, attr))
        # If we need to split a record, make as much copies as there are crossrefs and adjust each crossref to one of these
        if nb_copies > 1:
            split_counter += 1
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
            untouched_counter += 1
            setattr(h, attr, getattr(h, attr)[0])
            # Insert at the index of the old record
            sanitised_hits.append([h])
        # If it's empty, don't include it in the new hit set
        else:
            discard_counter += 1
            continue
    
    # Flatten the full hit list again
    hits = list(it.chain(*sanitised_hits))
    
    LOG.debug(f"{split_counter} hit records have been exploded while sanitising.")
    LOG.debug(f"{discard_counter} hit records have been discarded while sanitising.")
    LOG.debug(f"{untouched_counter} hit records have been left untouched while sanitising.")
    
    return hits


class Search(ABC):
    
    def __init__(self, query, params = {}, hits = [], clusters = [], 
                output_folder = Path('.'), temp_folder = Path('.')):
        
        self.query: list = query # list of query filepaths
        self.params: dict = params # dictionary containing the search configuration
        self.hits: list = hits # list of Hit objects
        self.clusters: list = clusters # list of Cluster objects
        
        self.OUTPUT_DIR: Path = output_folder
        self.TEMP_DIR: Path = temp_folder
        LOG.debug(f'Created temporary folder at {self.TEMP_DIR}.')
        
        return None
    
    
    def __repr__(self) -> str:
        return f"Search of {','.join(list(self.query.keys()))} with {len(self.clusters)} clusters identified"
    
    @abstractmethod
    def run(self):
        pass
    
    @abstractmethod
    def run_foldseek(self):
        pass
    
    @abstractmethod
    def parse_foldseek_results(self):
        pass
    
    
    def identify_clusters(self) -> None:
        """
        Identifies the gene clusters among the hits.
        """
        
        ### Load the requirements from params
        max_gap: int = self.params['max_gap']
        max_length: int = self.params['max_length']
        min_hits: int = self.params['min_hits']
        min_covered_queries: int = self.params['min_cov_qrs']
        require: list = self.params['require']
        
        LOG.debug('Applying the following cluster identification criteria:')
        LOG.debug(f'maximum intergenic gap >= {max_gap}')
        LOG.debug(f'maximum cluster length >= {max_length}')
        LOG.debug(f'minimum number of hits in a cluster >= {min_hits}')
        LOG.debug(f'minimum number of queries covered by a cluster >= {min_covered_queries}')
        LOG.debug(f'queries required to present in a cluster <= {require}')
        
        ### Cluster identification
        ## First, make groups by scaffold
        LOG.info('Grouping hits by scaffold')
        scaff_groups = {}
        for h in self.hits:
            if h.scaff in scaff_groups.keys():
                scaff_groups[h.scaff].append(h)
            else:
                scaff_groups[h.scaff] = [h]
        
        ## Then, calculate the intergenic distance between all connections on the same scaffold
        ## and filter out the ones failing the intergenic threshold
        ## and filter out self-hits as these are not genuine collocalised genes
        LOG.info("Calculating intergenic distances")
        close_groups = []
        for _, hits in scaff_groups.items():
            # Calculate the intergenetic distances and find the self-hits
            pairs_to_test = list(it.combinations(hits, 2))
            self_hits = {pair: Hit.same_location(*pair) for pair in pairs_to_test}
            dists = {pair: Hit.intergenic_distance(*pair) for pair in pairs_to_test}
            
            # Apply the filtering
            dists = {k:v for k,v in dists.items() if v <= max_gap} # apply max gap criterium
            dists = {k:v for k,v in dists.items() if not self_hits[k]} # filter out self-hits
            
            # Collect
            if len(dists) > 0:
                close_groups.append(list(dists.keys()))
                
        if len(close_groups) == 0:
            LOG.error("No hit groups passed the distance criteria!")
            sys.exit()
        
        ## Then identify the clusters by finding chains of distance pairs on the same scaffold using a directed network graph
        ## Account for multi-hits and -crossrefs by generating all possible hit chains when encountering pairs on the same genomic location
        LOG.info("Identifying gene clusters from chains of distance pairs")
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
                # Keep only the longest paths (discard the paths with shortcuts skipping a gene)
                max_path_length = max([len(p) for p in all_paths_this_chain])
                all_paths_this_chain = [p for p in all_paths_this_chain if len(p) == max_path_length]
                all_paths.append(all_paths_this_chain)
            
            # Save the chains for this neighbour group
            all_paths = list(it.chain(*all_paths))
            clusters.append(all_paths)
        
        # Flatten out all results
        clusters = list(it.chain(*clusters))
        
        ### Apply intra-cluster filtering requirements
        LOG.info("Filtering identified clusters")
        # Minimum number of hits in a cluster
        clusters_filt = [cl for cl in clusters if len(cl) >= min_hits]
        # Minimum number of covered queries and required queries
        covered_queries = [{h.query for h in cl} for cl in clusters_filt]
        clusters_filt = [cl for cl,qrs in zip(clusters_filt, covered_queries)
                         if len(qrs) >= min_covered_queries and set(require) <= qrs]
        
        ### Create the Cluster objects from the filtered hit clusters
        res_objects = [Cluster(cl, number = idx) for idx,cl in enumerate(clusters_filt)]
        
        ### Filter for maximum cluster length
        res_objects_filt = [cl for cl in res_objects if cl.length <= max_length]
        
        ## Rank by cluster score and renumber
        LOG.info('Sorting and renumbering by cluster score')
        res_objects_filt.sort(key = operator.attrgetter('score'), reverse = True)
        for idx,cl in enumerate(res_objects_filt):
            cl.number = idx+1
        
        ### Save
        self.clusters = res_objects_filt
        LOG.info(f"Identified {len(res_objects_filt)} gene clusters passing the criteria")
        
        ### Update the hits attribute after filtering at cluster level
        LOG.debug('Discarding hits not present in the identified gene clusters')
        self.hits = [h for cl in self.clusters for h in cl.hits]
        
        return None
   
    
    def generate_tables(self, output_folder: Path) -> None:
        """
        Saves the hit and cluster lists in separate overview tables.
        """
        # First the hits
        LOG.debug('Generating hit table')
        all_hit_data = [h.as_dict() for h in self.hits]
        all_hit_data_df = pl.DataFrame(all_hit_data, schema = ['db_id', 'query', 'scaff', 'strand', 'coords', 'db', 'crossref_id',
                                                               'crossref_method', 'name', 'taxon_name', 'taxon_id', 'evalue', 
                                                               'score', 'seqid', 'qcov', 'tcov'])
        all_hit_data_df.write_csv(output_folder / 'hits.tsv', include_header = True, separator = "\t")
        
        # Then the clusters
        LOG.debug('Writing gene cluster table')
        all_cluster_data = [cl.as_dict() for cl in self.clusters]
        all_cluster_data_df = pl.DataFrame(all_cluster_data, schema = ['number', 'hits', 'start', 'end', 'length', 'score', 'scaff',
                                                                       'strand', 'taxon_name', 'taxon_id', 'hits'])
        all_cluster_data_df.write_csv(output_folder / 'clusters.tsv', include_header = True, separator = "\t")
        
        return None
    
    
    def generate_cblaster_session(self) -> Session:
        """
        Generates a cblaster session object
        """
        def get_sequence_length_from_cif(file: Path) -> int:
            """
            Determines the CDS sequence length of the protein structure in the input CIF file.
            """
            structure = MMCIF2Dict(file)
            res_ids = [int(i) for i in structure['_entity_poly_seq.num']]
            return max(res_ids)*3 # codon triplets
        
        def get_clusters_by_id(self, nbs: list) -> list:
            """
            Returns the cluster object that has the given cluster number.
            """
            return [cl for cl in self.clusters if cl.number in nbs]
        
        def get_taxon_name_from_taxon_id(self, txid: int) -> str:
            """
            Returns the taxon name associated with a given taxon ID.
            """
            return [cl.taxon_name for cl in self.clusters if cl.taxon_id == txid][0]
        
        session_dict = {}
        
        ### Queries field
        LOG.debug("Generating Queries field")
        session_dict['queries'] = list(self.query.keys())
        
        ### Query field
        LOG.debug('Generating Query field')
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
        LOG.debug("Generating Params field")
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
        LOG.debug("Generating Organisms field")
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
            # Create a new organism instance if there's no one for this taxon ID
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
                
            # Populate this scaffold with clusters for each cluster identified on this scaffold
            # Keep track of the number of hits covered on this scaffold for proper references in the subject links
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
                                                      'coverage': hit.tcov,
                                                      'evalue': hit.evalue,
                                                      'bitscore': hit.score}]
                    cblaster_this_subject['name'] = hit.crossref_id
                    cblaster_this_subject['ipg'] = hit.db_id
                    cblaster_this_subject['start'] = hit.start()
                    cblaster_this_subject['end'] = hit.end()
                    cblaster_this_subject['strand'] = int(f"{hit.strand}1")
                    cblaster_this_subject['sequence'] = None
                    this_scaffold['subjects'].append(cblaster_this_subject)
        
        ## Discard the taxon ID and scaffold ID nested index to get the cblaster session format
        cblaster_organisms = list(cblaster_organisms.values())
        cblaster_organisms_new = []
        for organism in cblaster_organisms:
            organism['scaffolds'] = list(organism['scaffolds'].values())
            cblaster_organisms_new.append(organism)
            
        session_dict['organisms'] = cblaster_organisms_new
        
        ### Construct the cblaster Session
        LOG.debug('Constructing the cblaster session')
        session = Session.from_dict(session_dict)
        
        return session
    
    