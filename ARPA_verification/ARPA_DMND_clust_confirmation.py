#!/usr/bin/env python3

###################################################################
#ARPA: Alignment-Free Residue Pangenome Analysis
#Arnav Lal, Ahmed M. Moustafa and Paul J. Planet

#University of Pennsylvania, Philadelphia, PA, 19104
#Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA 19104, USA
#Children's Hospital of Philadelphia, Philadelphia, PA, 19104
#Sackler Institute for Comparative Genomics, American Museum of Natural History, New York, NY 10024, USA

#Contact1: arnavlal@sas.upenn.edu
#Contact2: arnavlal@gmail.com
###################################################################

#%%% NEW STUFF

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
import collections
import seaborn as sns
import json

cluster_diamond = pd.DataFrame(columns= ["ARPA","DIAMOND", "ARPA_MCL", "ARPA_DIAMOND"] , index = df_diamond[1], dtype= float)


#%%diamond processing

df_diamond = pd.read_csv('/Diamond Pangenome Output/DIAMOND_FAA2500/s_aur_faa_2500_clust_out.txt', header=None) #import diamond
df_diamond = df_diamond[0].str.split('\t', expand=True)

diamond_index = collections.defaultdict(list)
index = 0
for i in tqdm(np.unique(df_diamond[0])):
    diamond_index[i].append(index)
    index +=1
   
df_diamond['index'] = ''
for i in tqdm(df_diamond.index):
    df_diamond['index'][i] = diamond_index[df_diamond[0][i]][0]

df_diamond = df_diamond.set_index(1)
cluster_diamond["DIAMOND"] = df_diamond["index"]

#%%ARPA Processing
df_ARPA = pd.read_csv('ARPA Pangenome Output/Clusters_ID.tab', header=None) #import arpa
df_ARPA = df_ARPA[0].str.split('\t', expand=True)
ARPA_dict = collections.defaultdict(int)
index = 0
for i in tqdm(df_ARPA.index):
    terms = collections.Counter(df_ARPA.loc[i])
    if terms[None] > 0:
        terms.pop(None)
    for j in terms:
        ARPA_dict[j] = index
    index += 1
    
ARPA_clust_df = pd.DataFrame(list(ARPA_dict.items()), columns=['seq_id','cluster_id'])
ARPA_clust_df = ARPA_clust_df.set_index('seq_id')

cluster_diamond["ARPA"] = ARPA_clust_df["cluster_id"]

#%%ARPA_MCL Processing
MCL = open('/Path to ARPA MCL results/Final_MCL_Clusters.json',) 
df_ARPA_MCL = json.load(MCL)
df_ARPA_MCL = pd.DataFrame.from_dict(df_ARPA_MCL,orient='index')
ARPA_MCL_dict = collections.defaultdict(int)
index = 0
for i in tqdm(df_ARPA_MCL.index):
    terms = collections.Counter(df_ARPA_MCL.loc[i])
    if terms[None] > 0:
        terms.pop(None)
    for j in terms:
        ARPA_MCL_dict[j] = index
    index += 1
    
ARPA_MCL_clust_df = pd.DataFrame(list(ARPA_MCL_dict.items()), columns=['seq_id','cluster_id'])
ARPA_MCL_clust_df = ARPA_MCL_clust_df.set_index('seq_id')

cluster_diamond["ARPA_MCL"] = ARPA_MCL_clust_df["cluster_id"]

#%%% ARPA DIAMOND
ARPA_dmd = open('/Path to ARPA Diamond results/Final_ARPA_align_Clusters.json',)
df_ARPA_dmd = json.load(ARPA_dmd)
df_ARPA_dmd = pd.DataFrame.from_dict(df_ARPA_dmd,orient='index')
ARPA_dmd_dict = collections.defaultdict(int)
index = 0
for i in tqdm(df_ARPA_dmd.index):
    terms = collections.Counter(collections.Counter(df_ARPA_dmd.loc[i].dropna()))
    if terms[None] > 0:
        terms.pop(None)
    for j in terms:
        ARPA_dmd_dict[j] = index
    index += 1
    
ARPA_dmd_clust_df = pd.DataFrame(list(ARPA_dmd_dict.items()), columns=['seq_id','cluster_id'])
ARPA_dmd_clust_df = ARPA_dmd_clust_df.set_index('seq_id')

cluster_diamond["ARPA_DIAMOND"] = ARPA_dmd_clust_df["cluster_id"]

#%%



new_df = cluster_diamond.dropna()
temp = np.array(new_df).astype(int)


from sklearn.metrics.cluster import rand_score
from sklearn.metrics.cluster import adjusted_rand_score

#Rand Score
rand_score(new_df["ARPA_MCL"], new_df["ARPA"])


adjusted_rand_score(new_df["ARPA"], new_df["DIAMOND"])
adjusted_rand_score(new_df["DIAMOND"], new_df["DIAMOND"])
adjusted_rand_score(new_df["ARPA_MCL"], new_df["DIAMOND"])
adjusted_rand_score(new_df["ARPA_DIAMOND"], new_df["DIAMOND"])

adjusted_rand_score(new_df["ARPA"], new_df["ARPA"])
adjusted_rand_score(new_df["DIAMOND"], new_df["ARPA"])
adjusted_rand_score(new_df["ARPA_MCL"], new_df["ARPA"])
adjusted_rand_score(new_df["ARPA_DIAMOND"], new_df["ARPA"])



from sklearn import metrics

#AMI
#metrics.adjusted_mutual_info_score(new_df["ARPA"], new_df["DIAMOND"])

#Homogeneity Score
metrics.homogeneity_score(new_df["DIAMOND"], new_df["ARPA_DIAMOND"])

#Completeness Score
metrics.completeness_score(new_df["DIAMOND"], new_df["ARPA_DIAMOND"])

#V-Measure Score
metrics.v_measure_score(new_df["ARPA"], new_df["DIAMOND"])
metrics.v_measure_score(new_df["DIAMOND"], new_df["DIAMOND"])
metrics.v_measure_score(new_df["ARPA_MCL"], new_df["DIAMOND"])
metrics.v_measure_score(new_df["ARPA_DIAMOND"], new_df["DIAMOND"])

metrics.v_measure_score(new_df["ARPA"], new_df["ARPA"])
metrics.v_measure_score(new_df["DIAMOND"], new_df["ARPA"])
metrics.v_measure_score(new_df["ARPA_MCL"], new_df["ARPA"])
metrics.v_measure_score(new_df["ARPA_DIAMOND"], new_df["ARPA"])


#Fowlkes Mallows Score
metrics.fowlkes_mallows_score(new_df["ARPA"], new_df["DIAMOND"])
metrics.fowlkes_mallows_score(new_df["DIAMOND"], new_df["DIAMOND"])
metrics.fowlkes_mallows_score(new_df["ARPA_MCL"], new_df["DIAMOND"])
metrics.fowlkes_mallows_score(new_df["ARPA_DIAMOND"], new_df["DIAMOND"])

metrics.fowlkes_mallows_score(new_df["ARPA"], new_df["ARPA"])
metrics.fowlkes_mallows_score(new_df["DIAMOND"], new_df["ARPA"])
metrics.fowlkes_mallows_score(new_df["ARPA_MCL"], new_df["ARPA"])
metrics.fowlkes_mallows_score(new_df["ARPA_DIAMOND"], new_df["ARPA"])

#%%
#%%
#%%






