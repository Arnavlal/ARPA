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



######################################################
############  CLUSTER COMPARISON METRICS #############
########  COMPARES PANAROO AND ARPA CLUSTERS #########
######################################################

#Imports
import numpy as np
import pandas as pd
from tqdm import tqdm
import collections


##MODIFY THESE FOLLOWING CODES TO THE DIRECTORY TO YOUR FILES
#files to path
df = pd.read_csv('/Path to Panaroo file/gene_presence_absence.csv', sep=",") #import Panaroo
df_A = pd.read_csv('/Path to ARPA/Clusters_ID.tab', header=None) #import arpa



df_new = df.drop(['Gene', 'Non-unique Gene name', 'Annotation'], axis=1) 
hash_Panaroo = df_new.applymap(hash)
hash_Panaroo1 = hash_Panaroo.to_numpy().astype(float)

#arpa editing
df_A = df_A[0].str.split('\t', expand=True)
ARPA_df = df_A
hash_arpa  = ARPA_df
hash_arpa = hash_arpa.applymap(hash)
hash_arpa1 = hash_arpa.to_numpy().astype(float)


#arpa createa table with unique protein id and cluster number (cluster number is integer)
arpa_hash_clust = np.zeros((len(np.where(hash_arpa1 != float(hash(None)))[0]),2))   
p=0
for i in tqdm(np.arange(0,len(hash_arpa1))):
    a = hash_arpa1[i,hash_arpa1[i,:]!=hash(None)]
    arpa_hash_clust[p:p+len(a),0] = i
    arpa_hash_clust[p:p+len(a),1] = a
    p+=len(a)

#roary createa table with unique protein id and cluster number (cluster number is integer)
panaroo_hash_clust = np.zeros((len(np.where(hash_Panaroo1 != float(hash(None)))[0]),2))   
p=0
for i in tqdm(np.arange(0,len(hash_Panaroo1))):
    a = hash_Panaroo1[i,hash_Panaroo1[i,:]!=hash(None)]
    panaroo_hash_clust[p:p+len(a),0] = i
    panaroo_hash_clust[p:p+len(a),1] = a
    p+=len(a)

del(hash_arpa)
del(hash_arpa1)
del(df)
del(df_A)
del(ARPA_df)
del(hash_Panaroo)
del(hash_Panaroo1)



#use the cluster_ids to two different clusters
dict_arpa_panaroo = collections.defaultdict(list)
for i in tqdm(np.arange(len(panaroo_hash_clust))):
    dict_arpa_panaroo[panaroo_hash_clust[i,1]].append(panaroo_hash_clust[i,0])   
for i in tqdm(np.arange(len(arpa_hash_clust))):
    dict_arpa_panaroo[arpa_hash_clust[i,1]].append(arpa_hash_clust[i,0])
 
keep_proteins_clustered_by_both = []
for i in tqdm(dict_arpa_panaroo.keys()):
    if len(dict_arpa_panaroo[i]) == 2:
        keep_proteins_clustered_by_both.append(dict_arpa_panaroo[i])
keep_proteins_clustered_by_both = np.array(keep_proteins_clustered_by_both)
ARPA_PANAROO_clust = keep_proteins_clustered_by_both
del(keep_proteins_clustered_by_both)


working_group2 = ARPA_PANAROO_clust

from sklearn.metrics.cluster import rand_score
from sklearn.metrics.cluster import adjusted_rand_score

#Rand Score
rand_score(working_group2[:,1], working_group2[:,0])

adjusted_rand_score(working_group2[:,1], working_group2[:,0])

from sklearn import metrics

#AMI
metrics.adjusted_mutual_info_score(working_group2[:,0], working_group2[:,1]) 

#Homogeneity Score
metrics.homogeneity_score(working_group2[:,1], working_group2[:,0]) 

#Completeness Score
metrics.completeness_score(working_group2[:,1], working_group2[:,0]) 

#V-Measure Score
metrics.v_measure_score(working_group2[:,1], working_group2[:,0]) 

#Fowlkes Mallows Score
metrics.fowlkes_mallows_score(working_group2[:,1], working_group2[:,0]) 

#Normalized Mutual Info
metrics.normalized_mutual_info_score(working_group2[:,1], working_group2[:,0]) 


