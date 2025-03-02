#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:26:13 2024

@author: arnavlal
"""

import sys
import argparse
import time
import pandas as pd
import numpy as np
import collections
from tqdm import tqdm
import pickle

start_time = time.time()
PARSER = argparse.ArgumentParser(prog="ARPA_MCL_formatter.py")
    
PARSER.add_argument("PKL_file", type=str)
PARSER.add_argument("MCL_file", type=str)
PARSER.add_argument("Output", type=str)


if len(sys.argv) == 1:
    sys.exit(0)
    
ARGS = PARSER.parse_args()

pkl = ARGS.PKL_file
mcl = ARGS.MCL_file
pkl_file_name = open(pkl,'rb')
pkl_file = pickle.load(pkl_file_name)
mcl_file = pd.read_csv(mcl, header = None, sep = '\t')

final_clust = collections.defaultdict(list)
final_clust_index = 0
for i in tqdm(np.arange(len(mcl_file))):
    terms = mcl_file.loc[i].dropna()
    for j in terms:
        final_clust[final_clust_index].append(pkl_file[int(j-1)])
    final_clust_index +=1

for i in final_clust.keys():
    final_clust[i] = [x for xs in final_clust[i] for x in xs] #from: https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
    final_clust[i] = [x for xs in final_clust[i] for x in xs]


import json
json.dump( final_clust, open(ARGS.Output+"/Final_MCL_Clusters.json", 'w' ) )


