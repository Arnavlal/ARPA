#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: arnavlal
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import collections
from Bio import SeqIO
from tqdm import tqdm
import time
start_time = time.time()


Folder = '/Path to GCA folder/NCBI-datasets/GCA_folder/'

Folder_list = []
for file in os.listdir(Folder):
    if not file.startswith('.'):   
        Folder_list.append(Folder + file)

AA = ["G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T"]

vrcs = collections.defaultdict(list) 

def FAAtoNUM(i):
    name = Folder_list[i] + str("/protein.faa")    
    pathToFile = open(name)
    index = 0
    for seq_record in SeqIO.parse(pathToFile, """fasta"""):
        if index == 0:
            name_label = seq_record.description.split("[")[1].split("]")[0]
        index +=1
        
        aa_count_comparison = np.zeros((20)) #empty array for the 20-mer
        Counts = collections.Counter(seq_record.seq) #count up each of the letters
        for j in np.arange(0,20): #for each of the 20 amino acids:
            aa_count_comparison[j] = Counts.get(AA[j]) or 0 #add how many of the amino acid residues are in the protein
        vrcs[name_label].append(aa_count_comparison) #add this 20-mer to the has value

print("Importing Genomes:") 
for i in tqdm(np.arange(len(Folder_list))): 
    FAAtoNUM(i) 

#%% Host Species

labels = []
results = np.zeros((len(vrcs),21))
index = 0
for i in vrcs.keys():
    test = np.array(vrcs[i])
    labels.append(i)
    total_counts  = len(test)
    #new_test = test[:,np.random.permutation(np.arange(20))]
    for j in np.arange(20):
        terms = test[:,0:j+1]
        unique_counts = len(np.unique(test[:,0:j+1], axis = 0))
        results[index,j+1] = unique_counts/total_counts
    index +=1


for i in np.arange(len(vrcs)):
    plt.plot(results[i,:], label = labels[i])
plt.grid()
plt.xticks(np.arange(1,21), labels = np.array(AA))
#plt.legend()


#%% E coli 2-d and 3-d

test = np.array(vrcs['Escherichia coli str. K-12 substr. MG1655'])
one_d = test[:,0:1]
two_d = terms = test[:,0:2]
three_d = terms = test[:,0:3]

plt.plot(one_d[:,0], color = "#bcbd22")

plt.plot(two_d[:,0],two_d[:,1], "o", alpha = 0.2, color ="#bcbd22")
plt.xlim(0,100)
plt.ylim(0,100)

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(12, 12))
ax = fig.add_subplot(projection='3d')
ax.scatter(three_d[:,0], three_d[:,1], three_d[:,2], alpha = 0.2, color = "#bcbd22")
plt.show()
ax.set_xlim(0,100)
ax.set_ylim(0,100)
ax.set_zlim(0,100)


#%% Switching up the order

test = np.array(vrcs['Escherichia coli str. K-12 substr. MG1655'])
total_counts = len(test)
# focus on E. coli
n_trial = 1000
results = np.zeros((n_trial,21))
for i in tqdm(np.arange(n_trial)):
    new_test = test[:,np.random.permutation(np.arange(20))]
    for j in np.arange(20):
        terms = new_test[:,0:j+1]
        unique_counts = len(np.unique(new_test[:,0:j+1], axis = 0))
        results[i,j+1] = unique_counts/total_counts
   
#for i in np.arange(n_trial):
#    plt.plot(results[i,:],'dodgerblue', alpha = 0.01, lw = 2)

test = np.array(vrcs['Drosophila melanogaster'])
total_counts = len(test)
# focus on E. coli
n_trial = 1000
results_dm = np.zeros((n_trial,21))
for i in tqdm(np.arange(n_trial)):
    new_test = test[:,np.random.permutation(np.arange(20))]
    for j in np.arange(20):
        terms = new_test[:,0:j+1]
        unique_counts = len(np.unique(new_test[:,0:j+1], axis = 0))
        results_dm[i,j+1] = unique_counts/total_counts

for i in np.arange(n_trial):
    plt.plot(results[i,:],color = "#bcbd22", alpha = 0.008, lw = 2)
for i in np.arange(n_trial):
    plt.plot(results_dm[i,:],'#d62728', alpha = 0.008, lw = 2)
plt.grid()
plt.xticks(np.arange(0,21))

#%% Utility of Amino Acids Counts


labels = []
results = np.zeros((len(vrcs),20))
index = 0
for i in vrcs.keys():
    test = np.array(vrcs[i])
    labels.append(i)
    total_counts  = len(test)
    for j in np.arange(20):
        terms = new_test[:,j]
        unique_counts = len(np.unique(test[:,j], axis = 0))
        results[index,j] = unique_counts/total_counts
    index +=1

#import pandas as pd
#results = pd.DataFrame(results, index = labels, columns=AA)

for i in np.arange(len(results)):
    plt.bar(np.arange(20),results[i,:], label = labels[i], alpha = 0.2)
plt.grid()

plt.xticks(np.arange(0,20), labels = np.array(AA))
#plt.legend()

plt.bar(np.arange(20), np.average(results, axis=0), fill=False) 
plt.errorbar(np.arange(20), np.average(results, axis=0), yerr=results.std(axis=0), fmt="o", color="black")




#%% Comparative distance calculations

test = np.array(vrcs['Escherichia coli str. K-12 substr. MG1655'])

diff = []
jaccard_differences = np.zeros((len(test), len(test)))
for i in tqdm(np.arange(len(test))):
    test_diff = np.sum(abs(test - test[i,:]), axis = 1)
    size = np.sum(test[i,:])
    diff.append(list(test_diff/size))

#https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
diff_cumsum = [x for xs in diff for x in xs]

# INITIAL GRAPH
"""np.histogram(diff_cumsum, bins = np.linspace(0,200,201))
hist = np.histogram(diff_cumsum, bins = np.linspace(0,2,21))
hist[0][0]  = hist[0][0] - len(test)
plt.grid()
plt.bar(np.linspace(0,2,20), hist[0], width = 0.09)
plt.bar(np.linspace(0,2,20), hist[0], width = 0.09, fill=False) 
plt.xlim(-0.05,2.05)
plt.ylim(0,3500000)

x_val = np.linspace(0,2,20)
y_vals = hist[0]


counter = 0
for i in x_val:
    plt.text(i-0.025,y_vals[counter]+60000,y_vals[counter], rotation = 90)
    counter +=1"""

terms_to_plot = []
loc_90 = []
pairwise_comp = np.array(diff)
for i in tqdm(np.arange(len(pairwise_comp))):
    for j in np.arange(len(pairwise_comp)):
        if i < j:
            terms_to_plot.append(min(pairwise_comp[i,j],pairwise_comp[j,i]))
            if min(pairwise_comp[i,j],pairwise_comp[j,i]) < 0.1:
                loc_90.append([i,j])

print(str(len(loc_90)) + " pairwise comparisons with similarity >90%")
print(str(1-np.median(terms_to_plot)) + " median similarity of proteins")



np.histogram(terms_to_plot, bins = np.linspace(0,200,201))
bins = np.linspace(0,1,11)
bins = np.append(bins, np.array([12]))
hist = np.histogram(terms_to_plot, bins = bins)
plt.grid()
plt.bar(np.linspace(0.05,1.05,11), hist[0], width = 0.09, align='center')
plt.bar(np.linspace(0.05,1.05,11), hist[0], width = 0.09, fill=False, align='center') 
plt.xlim(-0.05,1.05)
plt.ylim(0,2200000)

x_val = np.linspace(0.05,0.95,10)
y_vals = hist[0]

counter = 0
for i in x_val:
    plt.text(i-0.0145,y_vals[counter]+60000,y_vals[counter], rotation = 90)
    counter +=1

loc_ninety = np.array(loc_90)
np.savetxt('e_coli_blast.txt', loc_ninety, fmt='%d')

#open blast

with open('/Path to final blast/blast_final_e_ecoli.txt') as f:
    blast_lines = f.read().splitlines()

scores = []
for i in np.arange(len(blast_lines)):
    terms_str = blast_lines[i]
    terms_str_upd = terms_str.split("b")[1].split("\\n")
    terms_str_final = terms_str_upd[0].split("\\t")
    scores.append((float((terms_str_final[4])) * float(terms_str_final[5])/100)/max(float(terms_str_final[1]),float(terms_str_final[3])))
    

"""
    if len(terms_str_upd) > 2:
        terms_str_final = terms_str_upd[0].split("\\t")
        scores.append([float(terms_str_final[2]), float(terms_str_final[3])/(float(terms_str_final[7])- float(terms_str_final[6])+1)])
        #for j in np.arange(len(terms_str_upd)-1):
        #    terms_str_final = terms_str_upd[j].split("\\t")
    else:
        terms_str_final = terms_str_upd[0].split("\\t")
        scores.append([float(terms_str_final[2]), float(terms_str_final[3])/(float(terms_str_final[7])- float(terms_str_final[6])+1)])
"""
scores = np.array(scores)



false_hits = np.unique(loc_ninety[np.invert(np.asarray(scores[:,0] > 95) & np.asarray(scores[:,1] > .99)),:])

lengths = []
for i in np.arange(len(false_hits)):
    lengths.append(sum(test[false_hits[i],:]))


#%% Multi E coli

E_coli_files = ['/Path to E coli files/EHEC EDL933/GCA_000732965.1_ASM73296v1_protein.faa', '/Path to E coli files/K12 MG1655/GCA_000005845.2_ASM584v2_protein.faa', '/Path to E coli files/EHEC Sakai/GCA_000008865.2_ASM886v2_protein.faa', '/Path to E coli files/ETEC H10407/GCA_000210475.1_ASM21047v1_protein.faa']


vrcs_ecoli = collections.defaultdict(list) 

def FAAtoNUM(i):
    name = i 
    pathToFile = open(name)
    index = 0
    for seq_record in SeqIO.parse(pathToFile, """fasta"""):
        if index == 0:
            name_label = seq_record.description.split("[")[1].split("]")[0]
        index +=1
        
        aa_count_comparison = np.zeros((20)) #empty array for the 20-mer
        Counts = collections.Counter(seq_record.seq) #count up each of the letters
        for j in np.arange(0,20): #for each of the 20 amino acids:
            aa_count_comparison[j] = Counts.get(AA[j]) or 0 #add how many of the amino acid residues are in the protein
        vrcs_ecoli[name_label].append(aa_count_comparison) #add this 20-mer to the has value


for i in E_coli_files:
    FAAtoNUM(i)



test = np.concatenate((np.array(vrcs_ecoli['Escherichia coli O157:H7 str. EDL933']), np.array(vrcs_ecoli['Escherichia coli str. K-12 substr. MG1655']),np.array(vrcs_ecoli['Escherichia coli O157:H7 str. Sakai']),np.array(vrcs_ecoli['Escherichia coli ETEC H10407'])))

diff = []
jaccard_differences = np.zeros((len(test), len(test)))
for i in tqdm(np.arange(len(test))):
    test_diff = np.sum(abs(test - test[i,:]), axis = 1)
    size = np.sum(test[i,:])
    diff.append(list(test_diff/size))

diff_cumsum = [x for xs in diff for x in xs]

# INITIAL GRAPH
"""np.histogram(diff_cumsum, bins = np.linspace(0,200,201))
hist = np.histogram(diff_cumsum, bins = np.linspace(0,2,21))
hist[0][0]  = hist[0][0] - len(test)
plt.grid()
plt.bar(np.linspace(0,2,20), hist[0], width = 0.09)
plt.bar(np.linspace(0,2,20), hist[0], width = 0.09, fill=False) 
plt.xlim(-0.05,2.05)
plt.ylim(0,3500000)

x_val = np.linspace(0,2,20)
y_vals = hist[0]


counter = 0
for i in x_val:
    plt.text(i-0.025,y_vals[counter]+60000,y_vals[counter], rotation = 90)
    counter +=1"""

terms_to_plot = []
loc_90 = []
pairwise_comp = np.array(diff)
for i in tqdm(np.arange(len(pairwise_comp))):
    for j in np.arange(len(pairwise_comp)):
        if i < j:
            terms_to_plot.append(min(pairwise_comp[i,j],pairwise_comp[j,i]))
            if min(pairwise_comp[i,j],pairwise_comp[j,i]) < 0.1:
                loc_90.append([i,j])


bins = np.linspace(0,1,11)
bins = np.append(bins, np.array([12]))
hist = np.histogram(terms_to_plot, bins = bins)
plt.grid()
plt.bar(np.linspace(0.05,1.05,11), hist[0], width = 0.09, align='center')
plt.bar(np.linspace(0.05,1.05,11), hist[0], width = 0.09, fill=False, align='center') 
plt.xlim(-0.05,1.05)
plt.ylim(0,55000000)

x_val = np.linspace(0.05,0.95,10)
y_vals = hist[0]

counter = 0
for i in x_val:
    plt.text(i-0.0145,y_vals[counter]+450000,y_vals[counter], rotation = 90)
    counter +=1


#%%

test = np.array(vrcs_ecoli['Escherichia coli str. K-12 substr. MG1655'])

diff = []
jaccard_differences = np.zeros((len(test), len(test)))
for i in tqdm(np.arange(len(test))):
    test_diff = np.sum(abs(test - test[i,:]), axis = 1)
    size = np.sum(test[i,:])
    diff.append(list(test_diff/size))

#https://stackoverflow.com/questions/952914/how-do-i-make-a-flat-list-out-of-a-list-of-lists
diff_cumsum = [x for xs in diff for x in xs]

# INITIAL GRAPH
"""np.histogram(diff_cumsum, bins = np.linspace(0,200,201))
hist = np.histogram(diff_cumsum, bins = np.linspace(0,2,21))
hist[0][0]  = hist[0][0] - len(test)
plt.grid()
plt.bar(np.linspace(0,2,20), hist[0], width = 0.09)
plt.bar(np.linspace(0,2,20), hist[0], width = 0.09, fill=False) 
plt.xlim(-0.05,2.05)
plt.ylim(0,3500000)

x_val = np.linspace(0,2,20)
y_vals = hist[0]


counter = 0
for i in x_val:
    plt.text(i-0.025,y_vals[counter]+60000,y_vals[counter], rotation = 90)
    counter +=1"""

terms_to_plot = []
loc_90 = []
pairwise_comp = np.array(diff)
for i in tqdm(np.arange(len(pairwise_comp))):
    for j in np.arange(len(pairwise_comp)):
        if i < j:
            terms_to_plot.append(min(pairwise_comp[i,j],pairwise_comp[j,i]))
            if min(pairwise_comp[i,j],pairwise_comp[j,i]) < 0.1:
                loc_90.append([i,j])


np.histogram(terms_to_plot, bins = np.linspace(0,200,201))
bins = np.linspace(0,1,11)
bins = np.append(bins, np.array([12]))
hist = np.histogram(terms_to_plot, bins = bins)
plt.grid()
plt.bar(np.linspace(0.05,1.05,11), hist[0], width = 0.09, align='center')
plt.bar(np.linspace(0.05,1.05,11), hist[0], width = 0.09, fill=False, align='center') 
plt.xlim(-0.05,1.05)
plt.ylim(0,2200000)

x_val = np.linspace(0.05,0.95,10)
y_vals = hist[0]

counter = 0
for i in x_val:
    plt.text(i-0.0145,y_vals[counter]+60000,y_vals[counter], rotation = 90)
    counter +=1

#%% Multi E. coli Roary analysis
E_coli_files = ['/Path to E coli files/EHEC EDL933/GCA_000732965.1_ASM73296v1_protein.faa', '/Path to E coli files/K12 MG1655/GCA_000005845.2_ASM584v2_protein.faa', '/Path to E coli files/EHEC Sakai/GCA_000008865.2_ASM886v2_protein.faa', '/Path to E coli files/ETEC H10407/GCA_000210475.1_ASM21047v1_protein.faa']

seqs = collections.defaultdict(list)
e_coli_vrcs = collections.defaultdict(list) 

for i in E_coli_files:
    name = i
    pathToFile = open(name)
    index = 0
    for seq_record in SeqIO.parse(pathToFile, """fasta"""):
        #if index == 0:
        #    name_label = seq_record.description.split("[")[1].split("]")[0]
        #index +=1
        seqs[seq_record.id] = seq_record.seq
        aa_count_comparison = np.zeros((20)) #empty array for the 20-mer
        Counts = collections.Counter(seq_record.seq) #count up each of the letters
        for j in np.arange(0,20): #for each of the 20 amino acids:
            aa_count_comparison[j] = Counts.get(AA[j]) or 0 #add how many of the amino acid residues are in the protein
        e_coli_vrcs[seq_record.id] = aa_count_comparison

multi_roary  = '/Multi E coli Roary Run /clustered_proteins'
multi_roary_pd = pd.read_csv(multi_roary, header=None)

#roary editing
multi_roary_pd = multi_roary_pd[0].str.split(': ', expand=True)
multi_roary_pd = multi_roary_pd[1].str.split('\t', expand=True)


vrcs_order = list(e_coli_vrcs.keys())

seq_clust_id = collections.defaultdict(list) 
index = 0
for i in tqdm(np.arange(len(multi_roary_pd))):
    clust_terms = multi_roary_pd.loc[i][multi_roary_pd.loc[i].notna()].str.split("-",expand = True)[1]
    for j in clust_terms:
        seq_clust_id[j].append(index)
    index +=1

vrcs_order_Clust = np.zeros((1,len(vrcs_order)))
for i in np.arange(len(vrcs_order)):
    if len(seq_clust_id[vrcs_order[i]]) == 0:
        vrcs_order_Clust[0,i] = -1
    else:
        vrcs_order_Clust[0,i] = seq_clust_id[vrcs_order[i]][0]
        
vrcs_order = np.array(vrcs_order)[vrcs_order_Clust[0,:] != -1]
vrcs_order_Clust = vrcs_order_Clust[vrcs_order_Clust != -1]

vrcs_array = np.zeros((20,len(vrcs_order)))
for i in np.arange(len(vrcs_order)):
    vrcs_array[:,i] = e_coli_vrcs[vrcs_order[i]]

vrcs_array = vrcs_array.T

diff = []
jaccard_differences = np.zeros((len(vrcs_array), len(vrcs_array)))
for i in tqdm(np.arange(len(vrcs_array))):
    test_diff = np.sum(abs(vrcs_array - vrcs_array[i,:]), axis = 1)
    size = np.sum(vrcs_array[i,:])
    diff.append(list(test_diff/size))

diff_cumsum = [x for xs in diff for x in xs]

terms_to_plot = []
loc_90 = []
pairwise_comp = np.array(diff)
for i in tqdm(np.arange(len(pairwise_comp))):
    for j in np.arange(len(pairwise_comp)):
        if i < j:
            terms_to_plot.append(min(pairwise_comp[i,j],pairwise_comp[j,i]))
            if min(pairwise_comp[i,j],pairwise_comp[j,i]) < 0.1:
                loc_90.append([i,j])
del(pairwise_comp)
del(diff)

clust_compare = np.zeros((2,len(loc_90)))
for i in np.arange(len(loc_90)):
    term = loc_90[i]
    clust_compare[0,i] = vrcs_order_Clust[term[0]]
    clust_compare[1,i] = vrcs_order_Clust[term[1]]

np.sum(clust_compare[0,:] == clust_compare[1,:])


from itertools import combinations
pairwise_comps = 0
for i in tqdm(np.arange(len(multi_roary_pd))):
    terms = multi_roary_pd.loc[i][multi_roary_pd.loc[i].notna()]
    if len(terms)>1:
        pairwise_comps += len(list(combinations(np.arange(len(terms)), 2)))

"""
TRUE_POS = 25282
FALSE_POS = 3656
TRUE_NEG = 200979616
FALSE_NEG = 2721

ARPA HITS = 28938
ROARY (TRUE) HITS = 28003
TRUE ARPA AND ROARY HITS = INTERSECTION = 25282
ROARY HIT NOT IN ARPA = 28003 - 25282 = 2721 (FN)
NOT ROARY HIT BUT IN ARPA = ARPA HITS NOT IN ROARY = 28938 - 25282 = 3656 (FP)
TOTAL COMPARISONS (((20051*20051)/2)- (20051/2)) = 201011275
NOT ARPA HIT AND NOT ROARY HIT = TOTAL - UNION OF ARPA AND ROARY HIT = 201011275 - (28938 + 28003 - (25282)) = 200979616 (TN)

sensitivity = TP/(TP+FN)
specificity  = TN/(TN+FP)

FN: when arpa says NO HIT but roary says HIT
FP: when ARPA says HIT but roary says NO HIT"""
