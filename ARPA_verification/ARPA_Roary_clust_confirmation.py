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




#!!!!
#set your CD to the location of the Clustered_proteins and Clusters_ID  files.
#!!!!

######################################################
############  CLUSTER COMPARISON METRICS #############
######################################################

#COMPARISON OF FILES
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


##MODIFY THESE FOLLOWING CODES TO THE DIRECTORY TO YOUR FILES

df = pd.read_csv('/Path/To/The/File/clustered_proteins', header=None)
df = df[0].str.split(': ', expand=True)
df = df[1].str.split('\t', expand=True)
ROARY_df = df


df = pd.read_csv('/Path/To/The/File/Clusters_ID.tab', header=None)
df = df[0].str.split('\t', expand=True)
ARPA_df = df


hash_arpa  = ARPA_df
hash_arpa = hash_arpa.applymap(hash)

hash_roary  = ROARY_df
hash_roary = hash_roary.applymap(hash)

hash_arpa1 = hash_arpa.to_numpy().astype(float)
hash_roary1 = hash_roary.to_numpy().astype(float)

translation = np.zeros((len(np.where(hash_arpa1 != float(hash(None)))[0]),2))          
translation = np.zeros((len(np.where(hash_arpa1 != float(hash(None)))[0]),2))

arpa_hash_clust = np.zeros((len(np.where(hash_arpa1 != float(hash(None)))[0]),2))   
p=0
for i in np.arange(0,len(hash_arpa1)):
    a = hash_arpa1[i,hash_arpa1[i,:]!=hash(None)]
    arpa_hash_clust[p:p+len(a),0] = i
    arpa_hash_clust[p:p+len(a),1] = a
    p+=len(a)

roary_hash_clust = np.zeros((len(np.where(hash_roary1 != float(hash(None)))[0]),2))   
p=0
for i in np.arange(0,len(hash_roary1)):
    a = hash_roary1[i,hash_roary1[i,:]!=hash(None)]
    roary_hash_clust[p:p+len(a),0] = i
    roary_hash_clust[p:p+len(a),1] = a
    p+=len(a)

b=[]
#This following loop may require some time
for i in np.arange(len(arpa_hash_clust)):
    if len(np.where(roary_hash_clust[:,1] == arpa_hash_clust[i,1])[0])>0.9:
        arpa_hash_clust[i,1] = roary_hash_clust[np.where(roary_hash_clust[:,1] == arpa_hash_clust[i,1])[0],0]
    else:
        b.append(i)

working_group2 = arpa_hash_clust
working_group2 = np.delete(working_group2,b,axis=0)
working_group = working_group2

from sklearn.metrics.cluster import rand_score

#Rand Score
rand_score(working_group2[:,1], working_group2[:,0])

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


#Save file to use in R for adjusted Rand score

np.savetxt('arpa.txt',working_group[:,0].astype(int),fmt='%10.5f')
np.savetxt('roary.txt',working_group[:,1].astype(int),fmt='%10.5f')



######################################################
############ COMBINED CLUSTER ANALYSIS ###############
######################################################


just_arpa = len(np.unique(working_group[:,0]))
just_roary = len(np.unique(working_group[:,1]))
no_nq = len(np.unique(working_group, axis=0))

print(just_arpa)
print(just_roary)
print(no_nq)

#%%
######################################################
################# SUPP FIG 1 #########################
######################################################

#Main Supp. Fig. 1
arpa_count = np.bincount(working_group[:,0].astype(int))
roary_count = np.bincount(working_group[:,1].astype(int))
bins = np.linspace(0,300,301)
from matplotlib import pyplot
pyplot.hist([arpa_count, roary_count], bins, alpha=0.8, label=['ARPA', 'ROARY'])
pyplot.legend(loc='upper right')
