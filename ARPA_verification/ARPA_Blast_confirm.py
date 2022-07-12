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


#CODE USE
#Run this CLI program to generate the files needed for BLAST analysis and creation.
#Program structure is similar to ARPA_main.py

###System Imports
import os
import sys
import argparse
import os.path
import time
import pandas as pd
import numpy as np
import collections
import matplotlib.pyplot as plt
from Bio import SeqIO
from scipy.cluster import hierarchy
import subprocess
import pickle

#OS Work
###################################################################
start_time = time.time()
PARSER = argparse.ArgumentParser(
    prog="ARPA.py",
    description="Alignment-Free Residue Pangenomic Analysis HOMOLOG BLAST CONFIRMATION",
)
PARSER.add_argument(
    "-v",
    "--version",
    help="print version and exit",
    action="version",
    version="%(prog)s 1.0",
)
PARSER.add_argument(
    "-ct",
    "--cluster_threshold",
    type=float,
    help="threshold of similarity used in homolog clustering [Default: 0.95]",
)

PARSER.add_argument(
    "-ni",
    "--neighborhood_identity",
    type=int,
    help="threshold of genes to be reclustered in neighborhood analysis [Default: 0.7]",
)

PARSER.add_argument(
    "-nt",
    "--neighborhood_thresh",
    type=int,
    help="threshold of genes conserved within neighborhood (out of 10) [Default: 3]",
)


PARSER.add_argument(
    "-w",
    "--percent_identity",
    type=float,
    nargs="?",
    help="select a blastp percent identity cutoff value [80], range(0,100)",
)
PARSER.add_argument(
    "-c",
    "--percent_coverage",
    type=float,
    nargs="?",
    help="select a blastp percent coverage cutoff value [80], range(0,100)",
)




PARSER.add_argument(
    "sequences_folder", type=str, help="directory with .faa files for pangenome"
)
if len(sys.argv) == 1:
    PARSER.print_help()
    sys.exit(0)
ARGS = PARSER.parse_args()


try:
    GETVERSION = subprocess.Popen("blastp -version", shell=True, stdout=subprocess.PIPE).stdout
    VERSION = GETVERSION.read()
    print("Found blastp (version:{})".format(VERSION.decode().splitlines()[1]))
except:
    PARSER.exit(status=0, message="Error: blastp cannot be found\n")
if ARGS.percent_identity:
    PERCENT_IDENTITY_CUTOFF = ARGS.percent_identity
else:
    PERCENT_IDENTITY_CUTOFF = 80.0
if ARGS.percent_coverage:
    PERCENT_COOVERAGE_CUTOFF = ARGS.percent_coverage
else:
    PERCENT_COOVERAGE_CUTOFF = 80.0

if bool(vars(ARGS)["dont_split_paralogs"]) and bool(vars(ARGS)["neighborhood_identity"]):
    PARSER.exit(status=0, message="Error: You cannot use -ni with -dsp\n")
if bool(vars(ARGS)["dont_split_paralogs"]) and bool(vars(ARGS)["neighborhood_thresh"]):
    PARSER.exit(status=0, message="Error: You cannot use -nt with -dsp\n")


if ARGS.cluster_threshold:
    threshold_percent = ARGS.cluster_threshold
else:
    threshold_percent = 0.95
if threshold_percent < 0.0 or threshold_percent > 1.0:
    PARSER.exit(status=0, message="Error: cluster_threshold values should be in range [0-1]\n")

if ARGS.neighborhood_identity:
    CGN_thresh = ARGS.neighborhood_identity
else:
    CGN_thresh = 0.7
if CGN_thresh < 0.0 or CGN_thresh > 1.0:
    PARSER.exit(status=0, message="Error: neighborhood_identity values should be in range [0-1]\n")

if ARGS.neighborhood_thresh:
    NTHRES = ARGS.neighborhood_thresh-0.1
else:
    NTHRES = 2.9
if NTHRES < 0.0 or NTHRES > 10.0:
    PARSER.exit(status=0, message="Error: neighborhood_thresh values should be in range [0-10]\n")


###TIMING
TIMESTR = time.strftime("%Y%m%d_%H%M%S")
print("ARPA Initialization: --- %s seconds ---" % (time.time() - start_time))
print("[>------------------]")
print(" ")
print(" ")


#####create results folder######
OS_SEPARATOR = os.sep
os.mkdir("ARPA_results_{}".format(TIMESTR))
RESULTS_FOLDER = "ARPA_results_{}{}".format(TIMESTR,OS_SEPARATOR)

print("Created Results Folder: --- %s seconds ---" % (time.time() - start_time))
print("[>------------------]")
print(" ")
print(" ")


Folder = ARGS.sequences_folder
Folder_list = []
try:
    for file in os.listdir(Folder):
        if file.endswith(".faa"):
            Folder_list.append(Folder + file)
    if len(Folder_list) == 0:
        PARSER.exit(
            status=0,
            message="The directory did not have any query faa files\n",
        )
except:
    if Folder.endswith(".faa"):
        PARSER.exit(
            status=0,
            message="Pangenomic Analysis Requires multiple faa files; please direct code to directory with .faa files\n",
        )
    else:
         PARSER.exit(
            status=0,
            message="The directory did not have any query faa files\n",
        )
for file in Folder_list:
    file_obj = open(file, "r")
    line_check = file_obj.readline()
    if not line_check.startswith(">"):
        PARSER.exit(status=0, message="FAA file not in FASTA format\n")
###################################################################
########################  MAIN ANALYSIS   #########################
###################################################################
#%% 

#IMPORTING SEQUENCES
AA = ["G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T"] #list of Amino Acid Residues
os.listdir(Folder)
valid_files =[".faa"] #Look for protein files
num_genome = len(os.listdir(Folder))
def FAAtoNUM(i):
    pathToFile = open(os.path.join(Folder, os.listdir(Folder)[i])) #import file
    allSeqs = []
    for seq_record in SeqIO.parse(pathToFile, """fasta"""):
        allSeqs.append(str(seq_record.seq)) #convert file to list of seq-type objects
    np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
    aa_count_comparison = np.zeros((len(allSeqs),20)) #create array for amino acid counts
    for i in np.arange(0,len(allSeqs)):
        Counts = collections.Counter(allSeqs[i])
        for j in np.arange(0,20):
            aa_count_comparison[i,j] = Counts.get(AA[j]) or 0 
    return aa_count_comparison

def FAAtoSeq(i):
    pathToFile = open(os.path.join(Folder, os.listdir(Folder)[i])) #import file
    allSeqs = []
    for seq_record in SeqIO.parse(pathToFile, """fasta"""):
        allSeqs.append(str(seq_record.seq)) #convert file to list of seq-type objects
    np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)
    return allSeqs


def FAAtoSTR(i):
    pathToFile = open(os.path.join(Folder, os.listdir(Folder)[i])) #import file
    allSeqs = []
    for seq_record in SeqIO.parse(pathToFile, """fasta"""):
        allSeqs.append(str(seq_record.description).split(" ",1)[1].split(" [")[0])
    return allSeqs
print("Finished Initialization: --- %s seconds ---" % (time.time() - start_time)) #TIME: <1 second
print("[=>-----------------]")
print(" ")
print(" ")

placeholder = [0]*1 #file to contain all genes (listed with the 20 Amino Acid counts)
for i in np.arange(0,len(os.listdir(Folder))):
    placeholder.append(FAAtoNUM(i)) #append each new genome to the list with all gene counts
placeholder.pop(0) #delete the first of list

placeholder2 = [0]*1 #file to contain all genes (listed with the 20 Amino Acid counts)
for i in np.arange(0,len(os.listdir(Folder))):
    placeholder2.append(FAAtoSeq(i)) #append each new genome to the list with all gene counts
placeholder2.pop(0) #delete the first of list
sequences = np.hstack(placeholder2)

value_of_genome_breakpoints = np.zeros((1,len(os.listdir(Folder))))
for i in np.arange(0,len(os.listdir(Folder))):
    value_of_genome_breakpoints[0,i] = len(placeholder[i])#find the index where the proteins are from different genomes.
cumulative_array = np.vstack(placeholder) #combine the amino acid counts for all of the genomes.
cumulative_breakpoints = np.cumsum(value_of_genome_breakpoints) #find genome breakpoints at the ends of the genomes
index = np.zeros((len(cumulative_array))) 
index[0:cumulative_breakpoints[1].astype(int)]=0
for i in np.arange(1,len(cumulative_breakpoints)):
    index[cumulative_breakpoints[i-1].astype(int):cumulative_breakpoints[i].astype(int)]=i #tell the array which genes belong to which genome.
placeholder_names = [0]*1#file to contain all genes (listed with the 20 Amino Acid counts)
for i in np.arange(0,len(os.listdir(Folder))):
    placeholder_names.append(FAAtoSTR(i)) #append each new genome to the list with all gene counts
placeholder_names.pop(0) #delete the first of list
cumulative_array_str = np.hstack(placeholder_names)

print("Finished Sequence Import: --- %s seconds ---" % (time.time() - start_time)) #TIME
print("Imported " +str(len(cumulative_array))+ " encoded genes from "+str(len(os.listdir(Folder)))+" genomes")
print("[====>--------------]")
print(" ")
print(" ")


#ALLELE COMPRESSION
##remove duplicates for analysis.
#%%
genomes = pd.DataFrame(cumulative_array)
pangenome = genomes.drop_duplicates()
hash_pangenome=pangenome.apply(lambda x: hash(tuple(x)), axis = 1)
hash_pangenome = pd.DataFrame(hash_pangenome)
hash_pangenome = pd.DataFrame.set_index(hash_pangenome, [pd.Index(np.arange(0,len(hash_pangenome)))])
hash_pangenome=hash_pangenome.squeeze()
len_ca = len(cumulative_array)
del(cumulative_array)
hash_genome=pd.DataFrame(genomes.apply(lambda x: hash(tuple(x)), axis = 1))
count = collections.Counter(hash_genome)
hash_genome = hash_genome.squeeze()
hash_genome_translation = hash_genome.to_numpy().astype('float')
translation = np.zeros((4,len(hash_genome_translation)))
translation = np.zeros((4,len(hash_genome_translation)))
translation[0,:] = hash_genome_translation
translation[1,:] = np.arange(len(translation[0,:]))
hash_genome_set = pd.DataFrame.set_index(hash_genome, [pd.Index(index)])
hash_genome_set = hash_genome_set.squeeze()
AApangenome = pangenome.to_numpy().astype('float')
AApangenome_hash = np.zeros((len(AApangenome), 1))
AApangenome_hash[:,0] = hash_pangenome.to_numpy().astype('float')
AApangenome_count = np.zeros((len(AApangenome), 1))
for i in np.arange(len(AApangenome)):
    AApangenome_count[i,0] = count[hash_pangenome[i]]
AApangenome  = np.concatenate((AApangenome,AApangenome_hash), axis=1)
AApangenome  = np.concatenate((AApangenome,np.zeros((len(AApangenome_hash),1))), axis=1)
AApangenome  = np.concatenate((AApangenome,AApangenome_count), axis=1)

del(AApangenome_count)
#%%

print("Finished Allele Compression:--- %s seconds ---" % (time.time() - start_time))
print("Compression to  " +str(len(AApangenome_hash))+ " encoded genes")
print("Compression ratio is:  " +str(len_ca/len(AApangenome_hash)))
print("[======>------------]")
print(" ")
print(" ")

#CLUSTERING OF HOMOLOGOUS GROUPS
breakpoints = np.zeros((1,len(AApangenome)))
l=1
val=0
new_arr = AApangenome*0
points_prior = 0
for i in np.arange(0,len(AApangenome)):
    if len(AApangenome)>0.5:
        test = np.sum(abs(AApangenome[:,0:20]-AApangenome[0,0:20]),axis=1)
        thresh = (1-threshold_percent)*np.sum(AApangenome[0,0:20])
        competitive_thresh = (1-((threshold_percent)*(threshold_percent)))*((2-threshold_percent)*np.sum(AApangenome[0,0:20]))
        points = np.where(test<thresh)[0]
        points_2 = np.where((test<competitive_thresh) & (test>thresh))[0]
        points_dup = points
        if len(points_2)>0:
            first_point = 0
            for j in np.arange(len(points_2)):
                vals = AApangenome[points]
                ctest = np.sum(abs(vals[:,0:20]-AApangenome[points_2[j],0:20]),axis=1)
                cthresh = ((1-threshold_percent))*np.sum(AApangenome[points_2[j],0:20])
                if len(np.where(ctest<cthresh)[0])>0: #if this outer ring genome is 95% related to some inner core one,
                    if first_point ==0:
                        retest = np.zeros((1,23))
                        retest[0,:] = AApangenome[points_2[j],:]
                    else:
                        retest = np.vstack((retest, AApangenome[points_2[j],:]))
                    first_point +=1
                    points_dup = np.concatenate((points_dup, points_2[j]), axis=None) #add to OG list
            if first_point>0:
                AApangenome[points_dup,21] = l
                new_arr[points_prior:points_prior+len(points_dup),:] = AApangenome[points_dup]
                points_prior += len(points_dup)
                breakpoints[0,i] =  points_prior
                AApangenome = np.delete(AApangenome, points_dup, axis=0)
                if len(retest)>0:
                    c2test = np.sum(abs(retest[0,0:20]-AApangenome[:,0:20]),axis=1)
                    c2thresh = ((1-threshold_percent))*np.sum(retest[0,0:20])
                    c2competitive_thresh = (1-((threshold_percent)*(threshold_percent)))*((2-threshold_percent)*np.sum(retest[0,0:20]))
                    c2points = np.where(c2test<c2thresh)[0]
                    c2points_2 = np.where((c2test<c2competitive_thresh) & (c2test>c2thresh))[0]
                    for k in np.arange(len(c2points_2)):
                        val = np.sum(abs(AApangenome[c2points_2[k],0:20]-AApangenome[c2points,0:20]),axis=1)
                        val_thresh = (1-threshold_percent)*np.sum(AApangenome[c2points_2[k],0:20])
                        if len(np.where(val<val_thresh)[0])>0:
                            retest = np.vstack((retest, AApangenome[c2points_2[k],:]))
                    AApangenome[c2points,21] = l
                    new_arr[points_prior:points_prior+len(c2points),:] = AApangenome[c2points]
                    points_prior += len(c2points)
                    breakpoints[0,i] =  points_prior
                    AApangenome = np.delete(AApangenome, c2points, axis=0)
                    retest = np.delete(retest, 0, axis=0)
                else:
                    break
                l+=1
            else:
                AApangenome[points,21] = l
                new_arr[points_prior:points_prior+len(points),:] = AApangenome[points]
                points_prior += len(points)
                breakpoints[0,i] =  points_prior
                AApangenome = np.delete(AApangenome, points, axis=0)
                l+=1
        else:
            AApangenome[points,21] = l
            new_arr[points_prior:points_prior+len(points),:] = AApangenome[points]
            points_prior += len(points)
            breakpoints[0,i] =  points_prior
            AApangenome = np.delete(AApangenome, points, axis=0)
            l+=1
    else:
        break

print("Proteins Clustered:--- %s seconds ---" % (time.time() - start_time)) 
print("Created " +str(l)+" Clusters")
print("[============>------]")
print(" ")
print(" ")




#consider using representative of OG and checking the cluster identities.
most_represented = np.zeros((int(new_arr[-1,21]), 24))
for i in np.arange(new_arr[-1,21]):
    genome = new_arr[np.where(new_arr[:,21]==i+1)[0],:]
    m_r = np.where(genome[:,22]==max(genome[:,22]))[0]
    if len(m_r)>1:
        most_represented[int(i),0:23] = genome[int(m_r[0]),:]
    else:
        most_represented[int(i),0:23] = genome[int(m_r),:]
    most_represented[int(i),23] = i+1
    

#Most Represented [Default]
for i in np.arange(new_arr[-1,21]):
    new_arr[np.where(new_arr[:,21]==most_represented[int(i),21])[0],22] = np.sum(abs(most_represented[int(i),0:20] - new_arr[np.where(new_arr[:,21]==most_represented[int(i),21])[0],0:20]), axis=1)

print("Generated Allele Scores: --- %s seconds ---" % (time.time() - start_time)) #TIME

test_keys = list(new_arr[:,20])
test_values = list(new_arr[:,21])
test_values2 = list(new_arr[:,22])
res = {test_keys[i]: test_values[i] for i in range(len(test_keys))}
res2 = {test_keys[i]: test_values2[i] for i in range(len(test_keys))}

for i in np.arange(len(translation[0,:])):
    translation[2,i] = res[translation[0,i]]
    translation[3,i] = res2[translation[0,i]]

seq_info = translation[1:4,:].T 
seq_info = seq_info[seq_info[:, 0].argsort()]
placeholder = np.zeros((len(seq_info),1))
placeholder[:,0] = index
seq_w_genome = np.concatenate((seq_info, placeholder), axis=1)
              

print("Starting Addition of alle scores to File: --- %s seconds ---" % (time.time() - start_time)) #TIME
Homolog_Split = np.zeros((num_genome,int(max(translation[2,:]))))*np.nan
Homolog_metadata = np.zeros((3,int(max(translation[2,:])))).astype(str)
print("Starting Metadata: --- %s seconds ---" % (time.time() - start_time))
list_names = [[] for i in range(len(most_represented)+1)]
for i in np.arange(len(seq_w_genome[:,1])):
    list_names[int(seq_w_genome[i,1]-1)].append(cumulative_array_str[int(seq_w_genome[i,0])])
for i in np.arange(len(most_represented)):
    ue,count_terms = np.unique(list_names[i], return_counts=True)
    sorted_b = ue[np.argsort(count_terms)[::-1]]
    if len(sorted_b)>3:
        Homolog_metadata[:,i-1] = sorted_b[0:3]
    else:
        if len(sorted_b)==2:
            Homolog_metadata[0:2,i-1] = sorted_b[0:2]
        else:
            Homolog_metadata[0:1,i-1] = sorted_b[0:1]
print("Ending Metadata: --- %s seconds ---" % (time.time() - start_time))
for i in np.arange(len(seq_w_genome)):
    Homolog_Split[int(seq_w_genome[i,3]),int(seq_w_genome[i,1])-1] = seq_w_genome[i,2]
hcount_nan = np.sum(np.isnan(Homolog_Split),axis=0)
HG_split_ordered = Homolog_Split[:,np.argsort(hcount_nan)]
HG_metadata = Homolog_metadata[:,np.argsort(hcount_nan)] 
dn = hierarchy.dendrogram(hierarchy.linkage(np.isnan(HG_split_ordered), 'single'))
print("Dendogram analysis done: --- %s seconds ---" % (time.time() - start_time)) #TIME
plt.close()
HG_split_ordered_DN = HG_split_ordered*0
for i in np.arange(len(dn['leaves'])):
    HG_split_ordered_DN[i,:] = HG_split_ordered[dn['leaves'][i],:]
print("Done with clustering.")

def countList(lst1, lst2):
    return [sub[item] for item in range(len(lst2))
                      for sub in [lst1, lst2]]

#cov  = [[] for i in range(int(max(seq_w_genome[:,1]))-1)]
#per_id = [[] for i in range(int(max(seq_w_genome[:,1]))-1)]

#outer_cov  = [[] for i in range(int(max(seq_w_genome[:,1]))-1)]
#outer_per_id = [[] for i in range(int(max(seq_w_genome[:,1]))-1)]

cov = []
per_id = []

outer_cov = []
outer_per_id = []

for i in np.arange(1,len(seq_w_genome[:,1])):
    group = sequences[seq_w_genome[np.where(seq_w_genome[:,1]==i)[0],0].astype(int)]
    if len(group)>1:
        group_names = list(np.arange(len(group)).astype(str))
        group_names = [">" + s for s in group_names]
        total_File = countList(group_names, group)
        with open("group.faa", "w") as f:
            for s in total_File:
                f.write(str(s) +"\n")
        
        ingroup_names = list(np.arange(1).astype(str))
        ingroup_names = [">" + s for s in ingroup_names]
        total_File = ingroup_names.append(group[0])
        with open("ingroup.faa", "w") as f:
            for s in ingroup_names:
                f.write(str(s) +"\n")
        
        rest = sequences[seq_w_genome[np.where(seq_w_genome[:,1]!=i)[0],0].astype(int)]
        rest_names = list(np.arange(len(rest)).astype(str))
        rest_names = [">" + s for s in rest_names]
        total_File = countList(rest_names, rest)
        with open("rest.faa", "w") as f:
            for s in total_File:
                f.write(str(s) +"\n")
        intra_blast = os.system("blastp -query {} -subject {} -evalue 1e-06 -max_target_seqs {} -max_hsps 1 -outfmt '6 qseqid sacc evalue qcovs pident' -out {}".format("ingroup.faa", "group.faa",num_genome, "ig_g.txt"))
        inter_blast = os.system("blastp -query {} -subject {} -max_target_seqs {} -max_hsps 1 -outfmt '6 qseqid sacc evalue qcovs pident' -out {}".format("ingroup.faa", "rest.faa",len_ca, "ig_r.txt"))
        a = np.loadtxt("ig_g.txt")
        if len(a)>1:
            cov.append(list(a[:,3]))
            per_id.append(list(a[:,4]))

        a2 = np.loadtxt("ig_r.txt")
        if len(a2)>1:
            if len(np.shape(a2))>1:
                outer_cov.append(list(a2[:,3]))
                outer_per_id.append(list(a2[:,4]))
            else:
                outer_cov.append(a2[3])
                outer_per_id.append(a2[4])
        print("Cluster: " + str(i))
    
print("saving")
with open("cov.pickle", "wb") as fp:   #Pickling
   pickle.dump(cov, fp)
print("saving coverage (cov)")
with open("per_id.pickle", "wb") as fp:   #Pickling
   pickle.dump(per_id, fp)
print("saving percent id (per id)")
with open("outer_cov.pickle", "wb") as fp:   #Pickling
   pickle.dump(outer_cov, fp)
print("saving coverage for proteins outside cluster (outer_cov)")
with open("outer_per_id.pickle", "wb") as fp:   #Pickling
   pickle.dump(outer_per_id, fp)
print("saving percentage id for proteins outside cluster (outer_per_id)")

sys.exit()
