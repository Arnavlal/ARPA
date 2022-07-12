#!/usr/bin/env python3

###################################################################
#ARPA: Alignment-Free Residue Pangenome Analysis
#Creating pangenomes without alignments since 2022

#Arnav Lal, Ahmed Moustafa and Paul J. Planet

#University of Pennsylvania, Philadelphia, PA, 19104
#Children's Hospital of Philadelphia, Philadelphia, PA, 19104

#Contact1: arnavlal@sas.upenn.edu
#Contact2: arnavlal@gmail.com
###################################################################
import os.path
import time
import pandas as pd
import numpy as np
import collections
import matplotlib.pyplot as plt
from Bio import SeqIO
inv_id = collections.defaultdict(list)
import itertools as it



#CHECK ARPA SINGLETONS; Supp. Fig 1 Inset
# 3 Inputs necessary for this code:
    
    
df = pd.read_csv('/Path/To/The/File/clustered_proteins', header=None)
df = df[0].str.split(': ', expand=True)
#df.drop(0, inplace=True, axis=1)
df = df[1].str.split('\t', expand=True)
ROARY_df = df


df = pd.read_csv('/Path/To/The/File/Clusters_ID.tab', header=None)
df = df[0].str.split('\t', expand=True)
ARPA_df = df


Folder ='Path/To/The/File/Special_genomes/'



#%% 
start_time = time.time()
#IMPORTING SEQUENCES
AA = ["G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T"] #list of Amino Acid Residues
os.listdir(Folder)
valid_files =[".faa"] #Look for protein files
seqs = collections.defaultdict(list)
vrcs = collections.defaultdict(list)
genomes = collections.defaultdict(list)
string = collections.defaultdict(list)
identity = collections.defaultdict(list)

num_genome = len(os.listdir(Folder))
def FAAtoNUM(i):
    pathToFile = open(os.path.join(Folder, os.listdir(Folder)[i])) #import file
    for seq_record in SeqIO.parse(pathToFile, """fasta"""):
        a = hash(str(seq_record.seq))
        genomes[a].append(i)
        if len(seqs[a])<1:
            seqs[a].append(str(seq_record.seq))
        if len(vrcs[a])<1:
            aa_count_comparison = np.zeros((20))
            Counts = collections.Counter(seq_record.seq)
            for j in np.arange(0,20):
                aa_count_comparison[j] = Counts.get(AA[j]) or 0 
            vrcs[a].append(aa_count_comparison)
        string[a].append(str(seq_record.description).split(" ",1)[1].split(" [")[0])
        identity[a].append(str(seq_record.description).split(" ",1)[0])  
                                                                                  
print("Finished Initialization: --- %s seconds ---" % (time.time() - start_time)) #TIME: <1 second
print("[=>-----------------]")
print(" ")
print(" ")

for i in np.arange(num_genome):
    FAAtoNUM(i)

print("Finished Sequence Import: --- %s seconds ---" % (time.time() - start_time)) #TIME
print("Imported " +str(sum([len(genomes[x]) for x in genomes if isinstance(genomes[x], list)]))+ " encoded genes from "+str(len(os.listdir(Folder)))+" genomes")

for key, values in identity.items():
    inv_id.update(zip(values, it.repeat(key)))


len_a = []
spcific_singleton =[]
mean = []
median = []


for i in np.arange(len(ARPA_df)):
    if sum(np.array(list(ARPA_df.iloc[i]))!=None)>1:
        print(i)
        continue
    elif sum(np.array(list(ARPA_df.iloc[i]))!=None)==0:
        print(i)
        continue
    else:
        if len(np.where(ARPA_df.iloc[i][0]==ROARY_df)[0])>0:
            spcific_singleton.append(sum(vrcs[inv_id[ROARY_df[np.where(ARPA_df.iloc[i][0]==ROARY_df)[1][0]][np.where(ARPA_df.iloc[i][0]==ROARY_df)[0][0]]]][0]))
            h = np.where(ARPA_df.iloc[i][0]==ROARY_df)[0][0]
            a = np.zeros((sum(np.array(list(ROARY_df.iloc[h]))!=None)))
            for j in np.arange(sum(np.array(list(ROARY_df.iloc[h]))!=None)):
                a[j] = sum(vrcs[inv_id[ROARY_df.iloc[h][j]]][0])
            len_a.append(len(a))
            mean.append(np.mean(a))
            median.append(np.median(a))
            print(i)
                
plt.hist(np.array(spcific_singleton)[np.where(np.array(len_a)!=1)[0]]/np.array(median)[np.where(np.array(len_a)!=1)[0]],bins=np.linspace(0,1.2,200),color='royalblue', rwidth=0.5)

