#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import operator
import itertools as it
import logging
import polars as pl
import networkx as nx
from abc import ABC, abstractmethod
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from cblaster.classes import Session


LOG = logging.getLogger(__name__)


class Hit:
    """
    Represents a single protein structure hit from a FoldSeek search.
    
    This class encapsulates information about a homologous protein structure match,
    including its database identifiers, sequence similarity metrics, genomic location,
    and taxonomic data.
    
    Attributes:
        query (str): ID of the homologous query protein.
        db_id (str): ID of the hit in its structure database.
        db (str): Structure database the hit was found in.
        crossref_id (list): ID used for cross-referencing (either KEGG or GenPept ID).
        crossref_method (str): Method used for cross-referencing (either KEGG, GenPept, WGS-GenPept, or local).
        name (str): Annotation or description of the hit.
        taxon_name (str): Name of the taxon in which this hit was found.
        taxon_id (int): Identifier of the taxon in which this hit was found.
        evalue (float): E-value of the FoldSeek hit.
        score (int): FoldSeek alignment score.
        seqid (float): Sequence identity percentage with the query protein.
        qcov (float): Query coverage percentage.
        tcov (float): Target coverage percentage.
        scaff (str): RefSeq or GenBank ID of the scaffold encoding the hit.
        coords (list): Genomic coordinates of the encoding gene's exons.
        strand (str): DNA strand the encoding gene is located on ('+' or '-').
    """
    
    def __init__(self, db_id, query, crossref_id = [], crossref_method = '', name = '', 
                 taxon_name = '', taxon_id = 0, db = "", evalue = 1, score = 0, 
                 seqid = 0, qcov = 0, tcov = 0, scaff = '', coords = [], strand = ''):
        """
        Initialise a Hit object with FoldSeek search results and genomic information.
        
        Args:
            db_id (str): ID of the hit in the structure database.
            query (str): ID of the homologous query protein.
            crossref_id (list, optional): Cross-reference ID. Defaults to [].
            crossref_method (str, optional): Cross-referencing method. Defaults to ''.
            name (str, optional): Hit annotation. Defaults to ''.
            taxon_name (str, optional): Taxon name. Defaults to ''.
            taxon_id (int, optional): Taxonomic ID. Defaults to 0.
            db (str, optional): Structure database name. Defaults to "".
            evalue (float, optional): E-value threshold. Defaults to 1.
            score (int, optional): FoldSeek score. Defaults to 0.
            seqid (float, optional): Sequence identity percentage. Defaults to 0.
            qcov (float, optional): Query coverage percentage. Defaults to 0.
            tcov (float, optional): Target coverage percentage. Defaults to 0.
            scaff (str, optional): Scaffold ID. Defaults to ''.
            coords (list, optional): Exon coordinates. Defaults to [].
            strand (str, optional): DNA strand. Defaults to ''.
        """
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
        """
        Convert the Hit object to a dictionary.
        
        Returns:
            dict: Dictionary with all Hit attributes; coordinates are formatted as
                  double-dot-separated exon pairs joined by commas.
        """
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
        """
        Return the start coordinate of the first exon.
        
        Returns:
            int | None: Minimum genomic coordinate across all exons, or None if
                       no coordinates are defined.
        """
        try:
            return min(it.chain(*self.coords))
        except ValueError:
            return None
        
    
    # Returns end coordinate of the last exon
    def end(self) -> int | None:
        """
        Return the end coordinate of the last exon.
        
        Returns:
            int | None: Maximum genomic coordinate across all exons, or None if
                       no coordinates are defined.
        """
        try:
            return max(it.chain(*self.coords))
        except ValueError:
            return None
    
    
    # Returns the sum of the exon lengths
    def length(self) -> int:
        """
        Return the total length in base pairs of all exons.
        
        Returns:
            int: Sum of lengths across all exons, calculated as
                 |end - start + 1| for each exon.
        """
        return sum([abs(c[1] - c[0] + 1) for c in self.coords])
    
    
    # Returns the intergenic distance between two genes. Negative if they overlap
    def intergenic_distance(self, other_hit: 'Hit') -> int:
        """
        Calculate the intergenic distance between this hit and another hit.
        
        For genes on the same scaffold, computes the distance between the end of
        the upstream gene and the start of the downstream gene. If genes overlap,
        returns the negative of the length of the overlapping gene.
        
        Args:
            other_hit (Hit): The other Hit object to measure distance to.
            
        Returns:
            int: Intergenic distance in base pairs (positive for gaps, negative
                 for overlaps). Returns the negative of the length of the smaller
                 gene in case of a full overlap.
        """
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
        """
        Check if two hits are at exactly the same genomic coordinates.
        
        Args:
            other_hit (Hit): The other Hit object to compare.
            
        Returns:
            bool: True if both hits are on the same scaffold and their genomic
                  coordinates completely overlap, False otherwise.
        """
        first = min([self, other_hit], key = operator.methodcaller('start'))
        last = max([self, other_hit], key = operator.methodcaller('start'))
        return last.start() == first.start() and last.end() == first.end() and first.scaff == last.scaff
    

class Cluster:
    """
    Represents a gene cluster containing one or more protein hits.
    
    A cluster groups proximal hits on the same genomic scaffold that meet
    specified clustering criteria. All hits in a cluster are expected to share
    the same scaffold and taxon.
    
    Attributes:
        hits (list[Hit]): List of Hit objects in the cluster.
        number (int): Cluster identifier/rank number.
        score (int): Cumulative score of all hits in the cluster.
        start (int): Minimum genomic coordinate across all hits.
        end (int): Maximum genomic coordinate across all hits.
        length (int): Total length in base pairs of all exons across hits.
        scaff (str): Scaffold/contig ID (taken from first hit).
        taxon_id (str): Taxonomic ID (taken from first hit).
        taxon_name (str): Taxonomic name (taken from first hit).
    """
    def __init__(self, hits, number = 0):
        """
        Initialise a Cluster object with a list of hits.
        
        Args:
            hits (list[Hit]): List of Hit objects to cluster together.
            number (int, optional): Cluster identifier number. Defaults to 0.
            
        Raises:
            ValueError: If hits have non-unique scaffold, taxon ID, or taxon name attributes.
        """
        self.hits: list[Hit] = hits
        self.number: int = number
        
        # Calculate cluster scores by summing hit scores
        self.score: int = sum([h.score for h in self.hits])
        
        # Cluster coordinates are defined as the most extreme coordinates
        self.start: int = min([h.start() for h in self.hits])
        self.end: int = max([h.end() for h in self.hits])
        
        # Cluster length is defined by the sum of hits' exons
        self.length: int = sum([h.length() for h in self.hits])
        
        # Take over the shared scaffold ID if it's unique among the hits
        common_scaff = {h.scaff for h in self.hits}
        if len(common_scaff) > 1:
            msg = f"Different scaffolds found among the gene hits in cluster {' '.join([h.db_id for h in self.hits])}"
            LOG.error(msg)
            raise ValueError(msg)
        else:
            self.scaff: str = self.hits[0].scaff
        
        # Same for taxon ID
        common_taxon_id = {h.taxon_id for h in self.hits}
        if len(common_taxon_id) > 1:
            msg = f"Different taxon IDs found among the gene hits in cluster {' '.join([h.db_id for h in self.hits])}."
            LOG.error(msg)
            raise ValueError(msg)
        else:
            self.taxon_id: str = self.hits[0].taxon_id
            
        common_taxon_name = {h.taxon_name for h in self.hits}
        if len(common_taxon_name) > 1:
            msg = f"Different taxon names found ammong the gene hits in cluster {' '.join([h.db_id for h in self.hits])}."
            LOG.error(msg)
            raise ValueError(msg)
        else:
            self.taxon_name: str = self.hits[0].taxon_name
        
        return None
    
    
    def __repr__(self) -> str:
        return f"Cluster {self.number}: {len(self.hits)} proteins from {self.scaff} ({self.start} - {self.end}), ({self.strand})\tScore: {self.score}"
    
    
    def as_dict(self) -> dict:
        """
        Convert the Cluster object to a dictionary.
        
        Returns:
            dict: Dictionary with cluster attributes including comma-separated
                  hit IDs and all genomic coordinates.
        """
        return {'hits': ','.join([h.db_id for h in self.hits]),
                'number': self.number,
                'score': self.score,
                'start': self.start,
                'end': self.end,
                'length': self.length,
                'scaff': self.scaff,
                'taxon_id': self.taxon_id,
                'taxon_name': self.taxon_name}
    

class Search(ABC):
    """
    Abstract base class for protein structure searches with cluster identification.
    
    This class manages a FoldSeek-based search workflow (remote or local), including result parsing,
    cluster identification using graph-based algorithms, and output generation.
    Subclasses must implement abstract methods for search execution and result parsing.
    Methods for the cluster identification and output generation are implemented here and shared
    over all subclasses.
    
    Attributes:
        query (list): List of query protein structure file paths.
        params (dict): Search configuration parameters (e.g., max gap, min hits).
        hits (list[Hit]): All identified hits from the search.
        clusters (list[Cluster]): Identified gene clusters passing filters.
        OUTPUT_DIR (Path): Directory for output files.
        TEMP_DIR (Path): Directory for temporary files.
    """
    def __init__(self, query, params = {}, hits = [], clusters = [], 
                output_folder = Path('.'), temp_folder = Path('.')):
        """
        Initialise a Search object.
        
        Args:
            query (list): List of query file paths or query identifiers.
            params (dict, optional): Search parameters dictionary. Defaults to {}.
            hits (list, optional): Pre-loaded Hit objects. Defaults to [].
            clusters (list, optional): Pre-loaded Cluster objects. Defaults to [].
            output_folder (Path, optional): Output directory path. Defaults to '.'.
            temp_folder (Path, optional): Temporary directory path. Defaults to '.'.
        """
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
        """
        Execute the complete search workflow.
        
        Note:
            This method must be implemented by subclasses to orchestrate the entire
            search process including FoldSeek execution and result parsing.
        """
        pass
    
    @abstractmethod
    def run_foldseek(self):
        """
        Run the FoldSeek search tool.
        
        Note:
            This method must be implemented by subclasses to execute FoldSeek either
            remotely or locally with the appropriate parameters, input files and
            target databases.
        """
        pass
    
    @abstractmethod
    def parse_foldseek_results(self):
        """
        Parse FoldSeek output and populate the hits list.
        
        Note:
            This method must be implemented by subclasses to convert raw FoldSeek
            results from the webserver or a local command call into a list of Hit objects.
        """
        pass


    def identify_proximal_groups(self, max_gap: int) -> list[list[Hit]]:
        """
        Identify proximal groups among the hits.
        
        Calculates the distance between all genes on the same scaffold, discards self-hits
        and hit pairs that fail the intergenic distance threshold.
        
        Returns:
            close_groups (list(list[Hit])): Hit pairs of proximal hits that pass the intergenic distance threshold.
            
        Mutates:
            RuntimeError: If there are not hit groups passing the intergenic distance criteria.
        """
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
            
            # Collect if there's a group of proximal hits on this scaffold
            if len(dists) > 0:
                close_groups.append(list(dists.keys()))
                
        ## Abort if there are no proximal hit groups
        if len(close_groups) == 0:
            msg = "No hit groups passed the distance criteria!"
            LOG.error(msg)
            raise RuntimeError(msg)
            
        return close_groups
    
    
    def identify_clusters_from_groups(self, close_groups: list, max_length: int, min_hits: int,
                                      min_covered_queries: int, require: set[str], all_layouts: bool):
        """
        Identify clusters from the proximal hit groups.
        
        Constructs Cluster objects from hit groups that pass all cluster identification thresholds (max cluster length,
        minimum no. hits, minimum no. covered queries, required queries).
        Can also return all cluster layouts that fit the cluster identification thresholds with a less-than-best score.
        
        Returns:
            clusters (list[Cluster]): list of Cluster objects that pass all identification thresholds.
        """
        ## Identify the clusters by finding chains of distance pairs on the same scaffold using a directed network graph
        ## Account for multi-hits and -crossrefs by generating all possible hit chains when encountering pairs on the same genomic location
        LOG.info("Identifying gene clusters from chains of distance pairs passing cluster criteria")
        clusters = []
        for cg in close_groups:
            # Order every hit pair so from up- to downstream
            reordered_cg = [sorted(pair, key = operator.methodcaller('start')) for pair in cg]
            
            # Identify the hit chains
            G = nx.DiGraph()
            G.add_edges_from(reordered_cg)
            chains = list(nx.weakly_connected_components(G))
            
            # Then, identify all possible chains by generating chains for all multi-hit or -crossref combinations
            all_clusters = []
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
                
                # Create Cluster objects from these longest paths
                all_clusters_this_chain = [Cluster(p) for p in all_paths_this_chain]
                
                # Apply intra-cluster filtering criteria
                # Minimum number of hits
                all_clusters_this_chain_filt = [cl for cl in all_clusters_this_chain if len(cl.hits) >= min_hits]
                # Minimum number of covered queries
                covered_queries = {cl: {h.query for h in cl.hits} for cl in all_clusters_this_chain_filt}
                all_clusters_this_chain_filt = [cl for cl,cov_qrs in covered_queries.items()
                                                if len(cov_qrs) >= min_covered_queries and set(require) <= cov_qrs]
                # Maximum cluster length
                all_clusters_this_chain_filt = [cl for cl in all_clusters_this_chain_filt if cl.length <= max_length]
                
                # Keep only the best-scoring cluster layout if not all cluster layouts are requested
                if not(all_layouts) and len(all_clusters_this_chain_filt) > 0:
                    all_clusters_this_chain_filt = [max(all_clusters_this_chain_filt, key = operator.attrgetter('score'))]
                
                # Collect
                all_clusters.append(all_clusters_this_chain_filt)
            
            # Save the clusters for this genomic neighbourhood
            all_clusters = list(it.chain(*all_clusters))
            clusters.append(all_clusters)
        
        # Flatten out all results
        clusters = list(it.chain(*clusters))
        
        return clusters
    
    
    def identify_clusters(self) -> None:
        """
        Identify gene clusters among the hits based on clustering criteria.
        
        This method groups hits by scaffold, calculates intergenic distances,
        filters based on maximum gap and minimum hit thresholds, and uses a
        directed graph to identify chains of unique proximal hits. It then applies
        additional filters for cluster size, query coverage, and length before
        ranking clusters by score.
        
        The method populates self.clusters with identified Cluster objects and
        updates self.hits to contain only hits in identified clusters.
        
        Raises:
            RuntimeError: If no hit groups pass the distance criteria.
            RuntimeError: If no cluster could be identified among the hit groups.
        """
        
        ### Load the requirements from params
        max_gap: int = self.params['max_gap']
        max_length: int = self.params['max_length']
        min_hits: int = self.params['min_hits']
        min_covered_queries: int = self.params['min_cov_qrs']
        require: list = self.params['require']
        all_layouts: bool = self.params['all_layouts']
        
        LOG.debug('Applying the following cluster identification criteria:')
        LOG.debug(f'maximum intergenic gap >= {max_gap}')
        LOG.debug(f'maximum cluster length >= {max_length}')
        LOG.debug(f'minimum number of hits in a cluster >= {min_hits}')
        LOG.debug(f'minimum number of queries covered by a cluster >= {min_covered_queries}')
        LOG.debug(f'queries required to present in a cluster <= {require}')
        
        ### Cluster identification
        ## First find proximal hit groups
        close_groups = self.identify_proximal_groups(max_gap = max_gap)
        ## Then identify clusters
        clusters = self.identify_clusters_from_groups(close_groups, max_length, min_hits,
                                                      min_covered_queries,require, all_layouts)
        
        ## Abort if no clusters have been identified
        if len(clusters) == 0:
            msg = "No cluster could be identified!"
            LOG.error(msg)
            raise RuntimeError(msg)
        
        ## Rank overall by cluster score and add number
        LOG.info('Sorting and renumbering by cluster score')
        clusters.sort(key = operator.attrgetter('score'), reverse = True)
        for idx,cl in enumerate(clusters):
            cl.number = idx+1
        
        ### Save
        self.clusters = clusters
        LOG.info(f"Identified {len(clusters)} gene clusters passing the criteria")
        
        ### Update the hits attribute after filtering at cluster level
        LOG.debug('Discarding hits not present in the identified gene clusters')
        self.hits = [h for cl in self.clusters for h in cl.hits]
        
        return None
   
    
    def generate_tables(self, output_folder: Path) -> None:
        """
        Save hit and cluster lists as tab-separated value (TSV) tables.
        
        Generates two output files:
        - hits.tsv: Table of all hits with their properties.
        - clusters.tsv: Table of all clusters with their properties.
        
        Args:
            output_folder (Path): Directory where output tables will be written.
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
                                                                       'strand', 'taxon_name', 'taxon_id'])
        all_cluster_data_df.write_csv(output_folder / 'clusters.tsv', include_header = True, separator = "\t")
        
        return None
    
    
    def generate_cblaster_session(self) -> Session:
        """
        Generate a cblaster Session object from the identified clusters.
        
        Constructs a cblaster-compatible session containing all search results,
        organised in the same hierarchy as cblaster (by organism, scaffold, and cluster).
        This object can be saved and reloaded for interactive visualisation and analysis
        outside of cfoldseeker.
        
        Returns:
            Session (cblaster.Session): Session holding all information about the identified clusters.
        """
        def get_sequence_length_from_cif(file: Path) -> int:
            """
            Determine the CDS sequence length from a CIF structure file.
            
            Extracts the number of modelled residues from the CIF file header and
            converts to base pair length using codon triplet scaling.
            
            Args:
                file (Path): Path to the PDB CIF format file.
                
            Returns:
                int: Sequence length in base pairs (number of residues * 3).
            """
            structure = MMCIF2Dict(file)
            res_ids = [int(i) for i in structure['_entity_poly_seq.num']]
            return max(res_ids)*3 # codon triplets
        
        def get_clusters_by_id(self, nbs: list) -> list:
            """
            Retrieve cluster objects by their cluster numbers.
            
            Args:
                nbs (list): List of cluster numbers to retrieve.
                
            Returns:
                list: List of Cluster objects matching the requested numbers.
            """
            return [cl for cl in self.clusters if cl.number in nbs]
        
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
        taxon_names = [cl.taxon_name for cl in self.clusters]
        scaffs = [cl.scaff for cl in self.clusters]
        cl_nbs = [cl.number for cl in self.clusters]
        cl_df = pl.DataFrame({'taxon_ids': taxon_ids, 
                              'taxon_names': taxon_names,
                              'scaffolds': scaffs, 
                              'cluster_number': cl_nbs})
        
        # Group by taxon ID and scaffold ID, and cast into a dictionary
        grouped_cl_df = cl_df.group_by(['taxon_ids', 'taxon_names', 'scaffolds']).all()
        grouped_cl_dict = {(row[0], row[1], row[2]): get_clusters_by_id(self, row[3]) for row in grouped_cl_df.iter_rows()}
        
        ## Make the cblaster session fields inside out, i.e. populate organisms first with scaffolds (and other attributes),
        ## then populate the scaffolds with clusters and subjects, then populate the clusters with links to the subjects.
        for (txid, txname, scaff), clusters in grouped_cl_dict.items():
            # Create a new organism instance if there's no one for this taxon ID
            if (txid, txname) in cblaster_organisms.keys():
                this_organism = cblaster_organisms[(txid, txname)]
            else:
                this_organism = {'name': txname,
                                 'strain': "",
                                 'scaffolds': {}}
                cblaster_organisms[(txid, txname)] = this_organism
            
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
    
    