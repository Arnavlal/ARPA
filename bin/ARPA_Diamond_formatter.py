#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import argparse
import time
import pandas as pd
import numpy as np
import collections
import copy
from tqdm import tqdm

start_time = time.time()
PARSER = argparse.ArgumentParser(prog="ARPA_Diamond_formatter.py")
    
PARSER.add_argument("clust_file", type=str)
PARSER.add_argument("dmnd_file", type=str)
PARSER.add_argument("Output", type=str)


if len(sys.argv) == 1:
    sys.exit(0)
    
ARGS = PARSER.parse_args()

clst = ARGS.clust_file
dmnd = ARGS.dmnd_file
output = ARGS.Output
df_diamond = pd.read_csv(dmnd, header=None) #import diamond
df_diamond = df_diamond[0].str.split('\t', expand=True)
index = copy.deepcopy(np.unique(df_diamond[0]))
df_diamond = df_diamond.set_index(df_diamond[0])

clst_df =  pd.read_csv(clst, header=None)
clst_df = clst_df[0].str.split('\t', expand=True)

clusters = collections.defaultdict(list)
for i in np.arange(len(clst_df)):
    clusters[i].append(list(clst_df.loc[i].dropna()))

import json
final_clust = collections.defaultdict(list)
indexing = 0
for b in tqdm(index):
    a = []
    for i in list(np.unique(df_diamond.loc[b][1])):
        a.append(clusters[int(i)][0])
    a = [item for sublist in a for item in sublist]
    final_clust[indexing] = a
    indexing +=1
json.dump( final_clust, open(output+"/Final_ARPA_align_Clusters.json", 'w' ) )

    
print("Saved Clusters File:--- %s seconds ---" % (time.time() - start_time)) 
print("Created " +str(len(final_clust))+" Clusters")
print(" ")
print(" ")


sys.exit()


