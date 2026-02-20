#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from remote import RemoteSearch
from local import LocalSearch


# Queries
query = {
        'fold_sco5086_model_0': '/home/lucas/bin/cfoldseeker_old/AF3_models/fold_sco5086_model_0.cif',
        'fold_sco5087_model_0': '/home/lucas/bin/cfoldseeker_old/AF3_models/fold_sco5087_model_0.cif',
        'fold_sco5088_model_0': '/home/lucas/bin/cfoldseeker_old/AF3_models/fold_sco5088_model_0.cif',
        'fold_sco5089_model_0': '/home/lucas/bin/cfoldseeker_old/AF3_models/fold_sco5089_model_0.cif'
        }
    
# query = {
#           'ery1': '/home/lucas/bin/cfoldseeker/test/prot1.cif',
#           'ery2': '/home/lucas/bin/cfoldseeker/test/prot2.cif',
#           'ery3': '/home/lucas/bin/cfoldseeker/test/prot3.cif'
#           }
    
# Search parameters
params = {'mode': 'remote',
          'db': ['afdb-proteome', 'afdb-swissprot', 'afdb50'],
          # 'db': ['pdb100'],
          'max_eval': 1e-9,
          'min_prob': 0,
          'min_score': 0,
          'min_seqid': 0,
          'min_qcov': 90,
          'min_tcov': 90,
          'max_gap': 5000,
          'max_length': 50000,
          'min_hits': 2,
          'min_cov_qrs': 2,
          'require': []
        }
    
# s = RemoteSearch(query, "uniprot_kegg_genpept.gz", params = params)
# s.run()

s = LocalSearch(query, 'foldseek_db', 'toy_cds.txt', 'tmp', params = params, db_prefix= 'DB/localDB')
s.run()
