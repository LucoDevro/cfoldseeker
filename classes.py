#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import operator
import itertools as it


class Hit:
    
    def __init__(self, db_id, query, crossref_id = [], crossref_method = '', name = '', 
                 taxon_name = '', taxon_id = 0, db = "", evalue = 1, prob = 1, score = 0, 
                 seqid = 0, qcov = 0, tcov = 0, scaff = '', coords = [], strand = ''):
        
        self.query: str = query #ID of the homologous query protein
        
        # ID attributes
        self.db_id: str = db_id #ID of the hit in its DB
        self.db: str = db #Structure database the hit was found in
        self.crossref_id: list = crossref_id #ID used for crossreffing (either ID from KEGG or GenPept)
        self.crossref_method: str = crossref_method #Method used for crossreffing (either KEGG or GenPept)
        
        # FoldSeek hit properties
        self.name: list = name #name of the hit in the KEGG entry
        self.taxon_name: str = taxon_name #name of the taxon in which this hit was found
        self.taxon_id: int = taxon_id
        self.evalue: float = evalue #evalue of the FoldSeek hit
        self.prob: float = prob #FoldSeek hit probability score
        self.score: int = score #FoldSeek score
        self.seqid: float = seqid #Sequence identity with the query protein
        self.qcov: float = qcov #Query coverage
        self.tcov: float = tcov #Target coverage
        
        # Genomic properties
        self.scaff: str = scaff #RefSeq or GenBank ID of the scaffold encoding the hit
        self.coords: list = coords #list of genomic coordinates of the encoding gene's exons
        self.strand: str = strand #DNA strand the encoding gene is part from
    
    def __repr__(self):
        return f"{self.query} Hit {self.db_id}\t {self.scaff} {self.start()}-{self.end()} ({self.strand})"
    
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
    
    # Returns the intergenic distance between two genes. Negative if they overlap
    def intergenic_distance(self, other_hit):
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
    def same_spot(self, other_hit):
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
    
