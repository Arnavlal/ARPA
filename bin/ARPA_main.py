#!/usr/bin/env python3


###################################################################
##################      INITIALIZATION     ########################
###################################################################

# ALIGNMENT - FREE CLUSTERING WITH ARPA
#ARPA: Alignment-Free Residue Pangenome Analysis


#Arnav Lal, Ahmed Moustafa Andries Feder, and Paul J. Planet
#University of Pennsylvania, Philadelphia, PA, 19104
#Children's Hospital of Philadelphia, Philadelphia, PA, 19104

#Contact1: arnavlal@alumni.upenn.edu 
#Contact2: arnavlal@sas.upenn.edu
#Contact3: planetp@chop.edu
###################################################################

###System Imports
import os
import sys
import argparse
import os.path
import time
import pandas as pd
import numpy as np
import collections
import csv
import copy
import matplotlib.pyplot as plt
from Bio import SeqIO
from matplotlib.pyplot import figure
from tqdm import tqdm
from itertools import compress
from scipy.cluster.hierarchy import ward, leaves_list
from scipy.spatial.distance import pdist

###################################################################

start_time = time.time()
PARSER = argparse.ArgumentParser(
    prog="ARPA_main.py",
    description="Alignment-Free Residue Pangenomic Analysis")
PARSER.add_argument(
    "-v",
    "--version",
    help="print version and exit",
    action="version",
    version="%(prog)s 1.0")
PARSER.add_argument(
    "-ct",
    "--cluster_threshold",
    type=float,
    help="threshold of similarity used in homolog clustering [Default: 0.95]")
PARSER.add_argument(
    "-b",
    "--binary",
    action='store_true',
    help="give pangenomic output as presence/absence (alternative is numerical est. of protein difference)")
PARSER.add_argument(
    "-dsp",
    "--dont_split_paralogs",
    action='store_true',
    help="Do not split paralogous groups (no neighborhood analysis)")
PARSER.add_argument(
    "-ps",
    "--paralog_syntany_analysis_size",
    type=int,
    help="Int size of gene neighborhood analysis on each side of gene for paralog separation [Default: 5]")
PARSER.add_argument(
    "-gr",
    "--graphical_mode",
    action='store_true',
    help="Output graphical pangenome plots at the cost of significant additional running time")

PARSER.add_argument(
    "-ld",
    "--large_dataset",
    action='store_true',
    help="Pangenome file saved as a Feather file, not CSV")

################################################################

PARSER.add_argument("-o","--output_file", type=str, help="location of output folder")

PARSER.add_argument("sequences_folder", type=str, help="directory name with list of '.faa' files for pangenome analysis")

   
if len(sys.argv) == 1:
    PARSER.print_help()
    print(" --------------- ")
    print("")
    print("")
    print("Analysis runs may take the form of:")
    print("")
    print("ARPA_main.py -ct 0.95 -gr /user/path/to/faa/files/")
    print("")
    print("--- OR ---")
    print("")
    print("ARPA_main.py -o /loc/of/output/folder/ /user/path/to/faa/files/")
    print("")
    print("--- OR ---")
    print("")
    print("ARPA_main.py -dsp -gr /user/path/to/faa/files/")
    print("")
    print("")
    print(" --------------- ")
    print("")
    sys.exit(0)
    
ARGS = PARSER.parse_args()

if ARGS.paralog_syntany_analysis_size:
    neighborhood_size = ARGS.paralog_syntany_analysis_size
    if ARGS.paralog_syntany_analysis_size < 0.5:
        PARSER.exit(status=0, message="Invalid size of neighborhood for paralog syntany analysis [must be integer greater than 0]")
else:
    neighborhood_size = 5

if ARGS.cluster_threshold:
    threshold_percent = ARGS.cluster_threshold
else:
    threshold_percent = 0.95
if threshold_percent < 0.0 or threshold_percent > 1.0:
    PARSER.exit(status=0, message="Error: cluster_threshold values should be in range [0-1]\n")

###  create a list of the .faa files within the folder  ####

if ARGS.sequences_folder.endswith("/"): #if the folder name doesn't have "/" at the end, then:
    Folder = ARGS.sequences_folder #name of the folder from human input above
else:
    Folder = ARGS.sequences_folder
    Folder += str("/") #add the "/" to the end. 

Folder_list = [] #list that will have all the paths + names ("/path/to/file/file.faa") for all sequences for analysis
try:
    for file in os.listdir(Folder):
        if file.endswith(".faa"): #within the folder, just take .faa files
            if not file.startswith('._'):   
                Folder_list.append(Folder + file) #add the path to the file names.
    if len(Folder_list) == 0: #if there are no faa files:
        PARSER.exit(status=0, message="The directory did not have any query faa files\n") # then exit the code
except:
    if Folder.endswith(".faa"): #if someone called ARPA by specifying a FAA file, then exit and tell them:
        PARSER.exit(
            status=0,
            message="Pangenomic Analysis Requires multiple faa files; please direct code to directory with .faa files\n",
        )
    else: #if someone did something even wilder, then give the generic exit:
         PARSER.exit(
            status=0,
            message="The directory did not have any query faa files\n",
        )

dsp = ARGS.dont_split_paralogs
        
#####  check each faa file to make sure they are FAA (primairly look for astarting ">"  #####

print("Checking FAA files: --- %s seconds ---" % "{:.2f}".format(time.time() - start_time)) #Variables are defined and everything is ready to run

### Error in reading this file
for file in tqdm(Folder_list): #for each of the files in the list of paths + names,
    file_obj = open(file, "r") #open each file
    line_check = file_obj.readline()
    if not str(line_check).startswith(">"): #if some of them does not start with ">" then complain.   
        PARSER.exit(status=0, message="Not every FAA file is in FASTA format\n")



### timing
print(" ")
print(" ")
TIMESTR = time.strftime("%Y%m%d_%H%M%S") #Current Year, Month, Day, Hour, Minutes, and Seconds
print("ARPA Initialization: --- %s seconds ---" % (time.time() - start_time)) #End of the initialization
print(" ")
print(" ")

OS_SEPARATOR = os.sep
#####  create results folder  ######
if ARGS.output_file:
    if ARGS.output_file.endswith("/"):
        os.mkdir(str(ARGS.output_file) + "ARPA_results_{}".format(TIMESTR))
        RESULTS_FOLDER = str(ARGS.output_file) + "ARPA_results_{}{}".format(TIMESTR,OS_SEPARATOR)
    else:
        os.mkdir(str(ARGS.output_file) + str(OS_SEPARATOR) + "ARPA_results_{}".format(TIMESTR))
        RESULTS_FOLDER = str(ARGS.output_file) + str(OS_SEPARATOR) + "ARPA_results_{}{}".format(TIMESTR,OS_SEPARATOR)
else:
    
    os.mkdir("ARPA_results_{}".format(TIMESTR)) #make a new directory in the working directory to store analysis files
    RESULTS_FOLDER = "ARPA_results_{}{}".format(TIMESTR,OS_SEPARATOR) #results folder uses he current timestring to create a new file that is time-stamped
    
    print("Created Results Folder: --- %s seconds ---" % (time.time() - start_time)) #created the results folder
    print(" ")
    print(" ")


###TO LOWER MEMORY FOOTPRINT
del(file)

###################################################################
########################     IMPORT       #########################
###################################################################
#%% 
threshold = 1-threshold_percent
#IMPORTING SEQUENCES
AA = ["G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T"] #list of Amino Acid Residues
#os.listdir(Folder)
valid_files =[".faa"] #Look for protein files

seqs = collections.defaultdict(list)
vrcs = collections.defaultdict(list) #Create a dictionary with each of the 20-number of vrcs for each unique protein
genomes = collections.defaultdict(list) #list their genome information
string = collections.defaultdict(list) #list their string characteristics
identity = collections.defaultdict(list) #keep their technical identifying number

if dsp == False:
    location = collections.defaultdict(dict) #dictionary for arranging paralogs, with indices for the location of each gene within each genome 
    hash_paralog = collections.defaultdict(dict)
    identity_per_genome = collections.defaultdict(dict)


num_genome = len(Folder_list) #number of genomes for analysis

def FAAtoNUM(i): #define this function that adds one genome (FAA) to the collections default dicts defined above.
    pathToFile = open(Folder_list[i]) #import file with the open command
    if dsp == False:
        index_paralog = 0
    check_dict = collections.defaultdict(list)
    for seq_record in SeqIO.parse(pathToFile, """fasta"""): #use the SeqIO tool to parse the opened file, one protein at a time
        #print(seq_record.description)
        #print(str(seq_record.description).split(" ",1)[1].split(" [")[0])
        a = hash(str(seq_record.seq)) # create a hash value of the protein (identification)
        genomes[a].append(i) #using the has value identifier, add to the genomes dictionary under the protein's hash that one genome contains it
        string[a].append(str(seq_record.description))   ## FIX THIS #.split(" ",1)[1].split(" [")[0]) #identify the name of the protein and save to hash value
        
        
        ###FIX THE ABOVE^^^
        
        identity[a].append(str(seq_record.description).split(" ",1)[0])# identify the technical identifying number and save to the hash value
        check_dict[a].append(np.random.rand())
        if dsp == False:
            identity_per_genome[i][a] = str(seq_record.description).split(" ",1)[0]
            if len(check_dict[a])>1:
                temp = hash_paralog[i][a]
                b = list(np.append(temp,index_paralog))
                hash_paralog[i][a] = b
            else:
                hash_paralog[i][a] = [index_paralog]
            location[i][index_paralog] = a
            index_paralog+=1
        if len(seqs[a])<1: #Sequence name is not kept, but this and the following line could retain it.
            seqs[a].append(str(seq_record.seq))
        if len(vrcs[a])<1: #if there is not already a 20-mer of numbers of that hash value, then save that 20-mer
            aa_count_comparison = np.zeros((20)) #empty array for the 20-mer
            Counts = collections.Counter(seq_record.seq) #count up each of the letters
            for j in np.arange(0,20): #for each of the 20 amino acids:
                aa_count_comparison[j] = Counts.get(AA[j]) or 0 #add how many of the amino acid residues are in the protein
            vrcs[a].append(aa_count_comparison) #add this 20-mer to the has value

print("Finished Initialization: --- %s seconds ---" % "{:.2f}".format(time.time() - start_time)) #Variables are defined and everything is ready to run
print(" ")
print(" ")

if dsp == False:
    genome_lengths = np.zeros((num_genome))
    for i in np.arange(num_genome):
        genome_lengths[i] = len([1 for line in open(Folder_list[i]) if line.startswith(">")])

print("Importing Genomes:") #starting to actually import the genomes into the DefaultDict one at a time
for i in tqdm(np.arange(num_genome)): #for each of the genomes in the list of genome paths
    FAAtoNUM(i) #run the defined function; this will result in a number of dictionaries with information from all the pangenomes
print("Finished Sequence Import: --- %s seconds ---" % "{:.2f}".format(time.time() - start_time)) #Time at which import is done
print("Imported " +str(sum([len(genomes[x]) for x in genomes if isinstance(genomes[x], list)]))+ " encoded genes from "+str(len(os.listdir(Folder)))+" genomes")


###  Vrc processing   ####
amino_acid_group = pd.DataFrame.from_dict(vrcs,orient='index') #turn all of the 20-mer vrcs from dictionary into a pandas dataframe; each row contains a 20-mer array and is indexed with the hash values

AAgroup = np.stack(amino_acid_group.to_numpy()[:,0]) #keep just the vrcs without the hash values into a numpy array
AA_hash = np.zeros((len(AAgroup),1)).astype(int) #make an empty array to put in the associated hash values
AA_hash[:,0] = np.array(amino_acid_group.index.values.tolist()) #fill the empty hash array with all of the hash values
AApangenome  = np.concatenate((AAgroup.astype(int),AA_hash.astype(int)), axis=1).astype(int) # Cols 1-20 are the vrcs and Col 21 is the hash values
AApangenome  = np.concatenate((AApangenome.astype(int),np.zeros((len(AAgroup),1)).astype(int)), axis=1).astype(int) #add empty column
AApangenome  = np.concatenate((AApangenome.astype(int),np.zeros((len(AAgroup),1)).astype(int)), axis=1).astype(int) #add empty column
AApangenome  = np.concatenate((AApangenome.astype(int),np.zeros((len(AAgroup),1)).astype(int)), axis=1).astype(int) #add empty column
for i in np.arange(len(AApangenome)):
    AApangenome[i,23] = len(genomes[AApangenome[i,20]]) #column 22: number of times the genes show up (not the number of genomes)
AApangenome[:,21] = np.sum(AApangenome[:,0:20],axis=1) #column 20: size of protein

print("Finished Allele Metadata:--- %s seconds ---" % "{:.2f}".format(time.time() - start_time))
print("Compression to  " +str(len(vrcs))+ " encoded genes")
print("Compression ratio is:  " +str("{:.2f}".format(sum([len(genomes[x]) for x in genomes if isinstance(genomes[x], list)])/len(vrcs))))
print(" ")
print(" ")

###TO LOWER MEMORY FOOTPRINT

del(AA_hash)
del(AAgroup)
del(amino_acid_group)

###################################################################
##############     QUICK 'N DIRTY CLUSTERING       ################
###################################################################
#%%
#cut the proteins into ones of similar enough sizes

print("Creating Rapid and Inflated Pre-Clusters:--- %s seconds ---" % "{:.2f}".format(time.time() - start_time))
print(" ")
print(" ")

import collections
def flatten(x): #flatten a list of lists into one list
    if isinstance(x, collections.abc.Iterable):
        return [a for i in x for a in flatten(i)]
    else:
        return [x]

def NE_clust(Near_Exact,divisor,addition): #cluster proteins that are similar to one another
    NE = Near_Exact 
    NE_df = pd.DataFrame(np.array((NE[:,0:20]+addition)/divisor).astype(int)) #create a dataframe of the data, keeping only vRCs, modify by grouping them by approximate size, creating a new identifier
    NE_hash = NE_df.apply(lambda x: hash(tuple(x)), axis=1).to_numpy().astype(float) #create a hash of these new identifiers
    NE_dict = collections.defaultdict(list) #create a dictionary for these idientifiers
    for i in np.arange(len(NE_hash)): #for each of the hash values, 
        NE_dict[NE_hash[i]].append(Near_Exact[i,20]) #link the hash to the vRC in the dictionary
    indices = np.unique(NE_hash) #keep only unique hash values
    return(NE_dict, indices)

def dict_inv(NE_dict,indices): #how to invert a dictionary
    NE_dict_inv = collections.defaultdict(list)
    for i in np.arange(len(indices)):
        a = flatten(NE_dict[indices[i]])
        for j in np.arange(len(a)):
            NE_dict_inv[a[j]].append(indices[i])
    return(NE_dict_inv)
        

def NE_clust_w_comb(Near_Exact, NE_dict, OG_NE_dict_inv, divisor,addition):
    NE = Near_Exact
    NE_df = pd.DataFrame(np.array((NE[:,0:20]+addition)/divisor).astype(int)) #create a dataframe,  keeping only vRCs, modify by grouping them by approximate size, creating a new identifier
    NE_hash = NE_df.apply(lambda x: hash(tuple(x)), axis=1).to_numpy().astype(float) #create a hash of these new identifiers
    NE_dict_itv = collections.defaultdict(list) #create a new dictionary for these identifiers
    for i in np.arange(len(NE_hash)): #for each of the hash values, 
        NE_dict_itv[NE_hash[i]].append(Near_Exact[i,20]) #link the hash to the vRC in the dictionary
    indices = np.unique(NE_hash) #keep only unique hash values
    comparison = np.zeros((len(NE),3)).astype(int) 
    comparison[:,0] = NE[:,20]
    NEW_NE_dict_inv = dict_inv(NE_dict_itv,indices)
    for i in np.arange(len(NE)):
        comparison[i,1] = OG_NE_dict_inv[comparison[i,0]][0]
        comparison[i,2] = (-1)*NEW_NE_dict_inv[comparison[i,0]][0]
    #unq = np.unique(comparison[:,1:3], axis=0)
    clust = 1
    clustering = collections.defaultdict(list)
    cluster_id = collections.defaultdict(list)
    cluster_fixing = []
    for i in np.arange(len(comparison)):
        if len(cluster_id[comparison[i,1]])==0 and len(cluster_id[comparison[i,2]])==0:
            clustering[clust].append(comparison[i,0])
            cluster_id[comparison[i,1]].append(clust)
            cluster_id[comparison[i,2]].append(clust)                
            clust+=1
        elif len(cluster_id[comparison[i,1]])==1 and len(cluster_id[comparison[i,2]])==0:
            clustering[cluster_id[comparison[i,1]][0]].append(comparison[i,0])
            cluster_id[comparison[i,2]].append(cluster_id[comparison[i,1]][0])
        elif len(cluster_id[comparison[i,1]])==0 and len(cluster_id[comparison[i,2]])==1:
            clustering[cluster_id[comparison[i,2]][0]].append(comparison[i,0])
            cluster_id[comparison[i,1]].append(cluster_id[comparison[i,2]][0])
        elif len(cluster_id[comparison[i,1]])==1 and len(cluster_id[comparison[i,2]])==1:
            clustering[cluster_id[comparison[i,2]][0]].append(comparison[i,0])
            if cluster_id[comparison[i,1]][0] != cluster_id[comparison[i,2]][0]:
                cluster_fixing.append(np.array([cluster_id[comparison[i,1]][0],cluster_id[comparison[i,2]][0]]))
        else:
            print('major error~!')
            print('System aborting; error in preclustering')
            sys.exit()
            
    #cluster_combinations = collections.defaultdict(list)
    #for i in np.arange(len(cluster_fixing)):

    cf = np.vstack(cluster_fixing)
    e=0
    #deleted_cfs = collections.defaultdict(list)
    for i in np.arange(len(cf)):
        if cf[i,0]==cf[i,0]:
            continue
        else:
            clustering[cf[i,0]].append(clustering[cf[i,1]])
            e+=len(clustering[cf[i,1]])
            print(e)
            #deleted_cfs[cf[i,1]].append(cf[i,0])
            cf[cf==cf[i,1]] = cf[i,0]
            del clustering[cf[i,1]]
    cluster = collections.defaultdict(list)
    clust = 1
    for i in clustering.keys():
        cluster[clust].append(clustering[i])
        clust+=1
    ind = np.arange(clust-1)+1
    return(cluster, ind)

iterations = 5

itv_clust, itv_ind = NE_clust(AApangenome,iterations,0)
print("Iteration " +str(1)+" of "+str(iterations))
print("Formed " +str(len(itv_clust))+" Hypothesized Clusters")
print(" ")
print(" ")
itv_dict_inv = dict_inv(itv_clust,itv_ind)
for i in np.arange(1,iterations):
    itv_clust,itv_ind = NE_clust_w_comb(AApangenome, itv_clust, itv_dict_inv, iterations,i)
    itv_dict_inv = dict_inv(itv_clust,itv_ind)
    print("Iteration " +str(i+1)+" of "+str(iterations))
    print("Formed " +str(len(itv_clust))+" Hypothesized Clusters")
    print(" ")
    print(" ")

###TO LOWER MEMORY FOOTPRINT
del(AApangenome)
del(itv_dict_inv)
del(itv_ind)

print("Created Hypothesized Inflated Pre-Clusters:--- %s seconds ---" % "{:.2f}".format(time.time() - start_time))
print("Formed " +str(len(itv_clust))+" Clusters")
print(" ")
print(" ")


###################################################################
####################     TRUE CLUSTERING       ####################
###################################################################
#%%
#%%

NE_dict = itv_clust #take cluster list (list of hashes) and move it to different name
del(itv_clust) #delete the old one
clusters_checked = collections.defaultdict(list) #confirmed accurate clusters
reps = collections.defaultdict(list) #precluster representatives
current_clust  = len(NE_dict)+1 #current cluster number

print("Conducting True Clustering, Step 1 of 2:")
for i in tqdm(np.arange(len(NE_dict))):
    temp = flatten(NE_dict[i+1]) #get a list of the hashes corresponding to the proteins within the cluster
    if len(temp)>1: #if the number of genes within the cluster is >1, then:
        testing_metrics_1 = np.sum(vrcs[temp[0]][0]) #using the first protein within the cluster, testing metric 1 is the length of this first protein
        a = np.zeros((len(temp),20)) #a will become a stacked list of the proteins in the cluster
        for j in np.arange(len(temp)):#along the first index of a, 
            a[j,0:20] = vrcs[temp[j]][0] #add in the vrc for each of the proteins
        a_sum = np.sum(abs(a-vrcs[temp[0]][0]),axis=1) #a_sum is the difference between the proteins in the cluster to the 1st protein in the cluster
        testing_metrics_3 = np.max(a_sum)/testing_metrics_1 #testing metric 3 is the maximum percentage length away a protein is from the 1st; ideally <5%
        mtc=a_sum/testing_metrics_1 #mtc holds all % difference values
        if testing_metrics_3<=threshold: #if all is well and the entire cluster is within the [95%] threshold from the first protein, then:
            clusters_checked[i+1].append(temp) #add the cluster to the clusters_checked array
            if testing_metrics_3<(threshold/2): #if the  max % length is lower than 0.5*threshold, then it is the sole representative.
                reps[i+1].append(temp[0])
            else:
                reps[i+1].append(temp[0]) #if not, add the first sequence to the list of representatives, and then:
                for k in np.arange(len(temp)):
                    if testing_metrics_3>(threshold/2):
                        if len(np.shape(a_sum)) == 1:
                            val = np.argmax(a_sum)
                        else:
                            val = np.argmax((a_sum.min(axis=0)))
                        a_sum = np.vstack((a_sum, np.sum(abs(a-a[val]),axis=1)))
                        reps[i+1].append(temp[val])
                        testing_metrics_3 = np.max(a_sum.min(axis=0))/testing_metrics_1
                    else:
                        continue
        else:
            clusters_checked[i+1].append(list(compress(temp,a_sum/testing_metrics_1<=threshold)))
            temp2 = list(compress(temp,a_sum/testing_metrics_1<=threshold))
            testing_metrics_1 = np.sum(vrcs[temp2[0]][0]) #length
            a = np.zeros((len(temp2),20))
            for j in np.arange(len(temp2)):
                a[j,0:20] = vrcs[temp2[j]][0]
            a_sum = np.sum(abs(a-vrcs[temp2[0]][0]),axis=1)
            testing_metrics_3 = np.max(a_sum)/testing_metrics_1
            if testing_metrics_3<(threshold/2):
                reps[i+1].append(temp2[0])
            else:
                reps[i+1].append(temp2[0])
                for k in np.arange(len(temp2)):
                    if testing_metrics_3>(threshold/2):
                        if len(np.shape(a_sum)) == 1:
                            val = np.argmax(a_sum)
                        else:
                            val = np.argmax((a_sum.min(axis=0)))
                        a_sum = np.vstack((a_sum, np.sum(abs(a-a[val]),axis=1)))
                        reps[i+1].append(temp2[val])
                        testing_metrics_3 = np.max(a_sum.min(axis=0))/testing_metrics_1
                    else:
                        continue
            new_group = list(compress(temp, mtc>threshold))
            ng2 = new_group
            for l in np.arange(len(ng2)):
                if len(new_group)<1:
                    continue
                else:
                    temp = new_group
                    testing_metrics_1 = np.sum(vrcs[temp[0]][0]) #length
                    a = np.zeros((len(temp),20))
                    for j in np.arange(len(temp)):
                        a[j,0:20] = vrcs[temp[j]][0]
                    a_sum = np.sum(abs(a-vrcs[temp[0]][0]),axis=1)
                    testing_metrics_3 = np.max(a_sum)/testing_metrics_1 #max length away (%)
                    mtc=(a_sum/testing_metrics_1)>threshold
                    if testing_metrics_3<=threshold:
                        clusters_checked[current_clust].append(temp)
                        if testing_metrics_3<(threshold/2):
                            reps[current_clust].append(temp[0])
                            current_clust+=1
                        else:
                            reps[current_clust].append(temp[0])
                            for k in np.arange(len(temp)):
                                if testing_metrics_3>(threshold/2):
                                    if len(np.shape(a_sum)) == 1:
                                        val = np.argmax(a_sum)
                                    else:
                                        val = np.argmax((a_sum.min(axis=0)))
                                    a_sum = np.vstack((a_sum, np.sum(abs(a-a[val]),axis=1)))
                                    reps[current_clust].append(temp[val])
                                    testing_metrics_3 = np.max(a_sum.min(axis=0))/testing_metrics_1
                                else:
                                    continue
                            current_clust+=1
                        a=0
                        new_group = list(compress(temp,mtc))
                    else:
                        clusters_checked[current_clust].append(list(compress(temp,a_sum/testing_metrics_1<=threshold)))
                        mtc2 = (a_sum/testing_metrics_1)>threshold
                        temp2 = list(compress(temp,a_sum/testing_metrics_1<=threshold))
                        testing_metrics_1 = np.sum(vrcs[temp2[0]][0]) #length
                        a = np.zeros((len(temp2),20))
                        for j in np.arange(len(temp2)):
                            a[j,0:20] = vrcs[temp2[j]][0]
                        a_sum = np.sum(abs(a-vrcs[temp2[0]][0]),axis=1)
                        testing_metrics_3 = np.max(a_sum)/testing_metrics_1
                        mtc = a_sum/testing_metrics_1
                        if testing_metrics_3<(threshold/2):
                            reps[current_clust].append(temp2[0])
                            current_clust+=1
                        else:
                            reps[current_clust].append(temp2[0])
                            for k in np.arange(len(temp2)):
                                if testing_metrics_3>(threshold/2):
                                    if len(np.shape(a_sum)) == 1:
                                        val = np.argmax(a_sum)
                                    else:
                                        val = np.argmax((a_sum.min(axis=0)))
                                    a_sum = np.vstack((a_sum, np.sum(abs(a-a[val]),axis=1)))
                                    reps[current_clust].append(temp2[val])
                                    testing_metrics_3 = np.max(a_sum.min(axis=0))/testing_metrics_1
                                else:
                                    continue
                            current_clust+=1
                        new_group = list(compress(temp,mtc2))
                        a = len(new_group)
                        #mtc = np.vstack((mtc,a_sum/testing_metrics_1)).min(axis=0)
    else:
        clusters_checked[i+1].append(temp)
        reps[i+1].append(temp[0])

del(NE_dict)
del(a_sum)

tuple_dict = collections.defaultdict(list)
c=[]
for i in reps.keys():
    temp = reps[i]
    for j in np.arange(len(temp)):
        tuple_dict[np.sum(vrcs[temp[j]][0])].append([i,temp[j]])
        c.append(np.sum(vrcs[temp[j]][0]))
del(reps)
index = pd.DataFrame(c).drop_duplicates().astype('int64')
del(c)
index = index.sort_values(by=0, ascending = False)
index = index.reset_index()[0]
combinations = collections.defaultdict(list)
identifier2 = 0
ini = index[0]+10
first = 0
print("Conducting True Clustering, Step 2 of 2:")
for i in tqdm(index): #index is the list of lengths of rep sequences
    if i < ini:
        rep_len_sze = threshold*i #representative length deviation allowed
        sze = index[(index>=(i-rep_len_sze)) & (index<=(ini))].to_list() #find all rep lens that fit this
        length = 0 #number of reps to be considered
        for l in sze:
            length+=len(tuple_dict[l]) #tuple_dict actually contains the reps, so this is the addition of reps for analysis
        if length == 0:
            pass #finished if no reps
        else:
            first+=1
            analysis = np.zeros((length, 22))#array for comparison of vrcs
            count=0
            ini = i-rep_len_sze #this makes sure that any prior analyses are not repeated.
            for m in sze:
                temp = tuple_dict[m]
                for n in np.arange(len(temp)):
                    analysis[count,0:20] = vrcs[temp[n][1]][0] #add vrcs and other values to analysis framework
                    analysis[count,20] = temp[n][0]
                    analysis[count,21] = sum(vrcs[temp[n][1]][0])
                    count+=1
            #for future delete redundancy,and stop re-adding to dict. 
            for o in np.arange(len(analysis)):
                analysis2 = np.sum(abs(analysis[:,0:20]-analysis[o,0:20]),axis=1)
                indices_to_add = analysis[:,20][analysis2<(threshold*sum(analysis[o,0:20]))] #terms similar to the "o"th array term
                for p in np.arange(len(indices_to_add)):
                    if analysis[o,20] != indices_to_add[p]:
                        combinations[min(indices_to_add[p],analysis[o,20])].append(max(indices_to_add[p],analysis[o,20]))
            if first == 1:
                analysis3 = analysis
            else:
                for q in np.arange(len(analysis3)):
                    analysis4 = np.sum(abs(analysis[:,0:20]-analysis3[q,0:20]),axis=1)
                    indices_to_add = analysis[:,20][analysis4<(threshold*sum(analysis3[q,0:20]))]
                    for r in np.arange(len(indices_to_add)):
                        if analysis3[q,20] != indices_to_add[r]:
                            combinations[min(indices_to_add[r],analysis3[q,20])].append(max(indices_to_add[r],analysis3[q,20]))
                analysis3 = analysis

###TO LOWER MEMORY FOOTPRINT
del(analysis)
del(analysis2)
del(analysis3)
del(analysis4)
del(tuple_dict)

new_dict = collections.defaultdict(list)
for i in combinations.keys():
    new_dict[i] = list(set(combinations[i]))
del(combinations)

d=0
e = []
for i in new_dict.keys():
    a = new_dict[i]
    #b = flatten(clusters_checked[i+1])
    d +=len(a)

comparison = np.zeros((d,2)).astype(int)
ind = 0
for i in new_dict.keys():
    a = new_dict[i]
    for j in np.arange(len(a)):
        comparison[ind,0] = i
        comparison[ind,1] = a[j]
        ind+=1
count = 0
delete_array = []
for i in np.arange(len(comparison)):
    if comparison[i,0] == comparison[i,1]:
        count+=1
    else:
        clusters_checked[comparison[i,0]].append(clusters_checked[comparison[i,1]])
        delete_array.append(comparison[i,1])
        comparison[comparison==comparison[i,1]] = comparison[i,0]
del(comparison)
del(new_dict)

keep_array= []
for i in np.arange(len(clusters_checked)):
    if (i+1 in delete_array) == False:
        keep_array.append(i+1)
    else:
        continue
del(delete_array)
print(" ")
print(" ") 
print("Finished with Clustering:--- %s seconds ---" % "{:.2f}".format(time.time() - start_time))
print(str(len(keep_array))+" clusters formed")
print(" ")
print(" ") 

if dsp == True:
    with open(os.path.join(RESULTS_FOLDER,"Clusters_ID.tab"), 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        c = 0
        for b in np.arange(len(keep_array)):
            a = []
            for i in flatten(clusters_checked[keep_array[b]]):
                a.append((identity[i]))
            a = [item for sublist in a for item in sublist]
            c+=1
            tsv_writer.writerow(a)
    print("Saved Clusters File:--- %s seconds ---" % (time.time() - start_time)) 
    print("Created " +str(c)+" Clusters")
    print(" ")
    print(" ")
    

###################################################################
###############      PARALOG SEPARATION       #####################
###################################################################
#%% 
#Count number of paralogs


ini = []
split_clusters = collections.defaultdict(list)
identity_clusters = collections.defaultdict(list)
#sequence_clusters = collections.defaultdict(list)
index3 = 0
test_ind_time = 0


if dsp == False:
    print("Finding Relevant Paralogs:--- %s seconds ---" % "{:.2f}".format(time.time() - start_time))
    print(" ")
    print(" ") 
    inverted_clusters_checked = dict_inv(clusters_checked,keep_array)
    
    """
    def Paralog_hunter(test):
        location_dictionary = collections.defaultdict(list)
        for k in np.arange(len(test)):
            for m in np.arange(2*neighborhood_size+1):
                if test[k,1] < neighborhood_size:
                   try:
                       location_dictionary[k].append(inverted_clusters_checked[location[test[k,0]][test[k,1]+m-neighborhood_size]][0])
                   except:
                       new_val = m + test[k,1] - neighborhood_size + int(genome_lengths[int(test[k,0])])
                       location_dictionary[k].append(inverted_clusters_checked[location[test[k,0]][new_val]][0])
                else:  
                    try:
                        location_dictionary[k].append(inverted_clusters_checked[location[test[k,0]][test[k,1]+m-neighborhood_size]][0])
                    except:
                        new_val = m + test[k,1] - neighborhood_size - int(genome_lengths[int(test[k,0])])
                        location_dictionary[k].append(inverted_clusters_checked[location[test[k,0]][new_val]][0])
        return location_dictionary 
    
    """
    
    def Paralog_hunter(test):
        location_dictionary = collections.defaultdict(list)
        for k in range(len(test)):
            genome_number = test[k, 0]
            protein_location = test[k, 1]
            genome_length = int(genome_lengths[int(genome_number)])
            
            if genome_length < 2*neighborhood_size-0.1:
                for m in location[genome_number].keys():
                    protein = location[genome_number][m]
                    location_dictionary[k].append(inverted_clusters_checked[protein][0])
            else:
                for m in range(-neighborhood_size, neighborhood_size + 1):
                    index = protein_location + m
        
                    # Wrap around if index is out of bounds
                    if index < 0:
                        index += genome_length
                    elif index >= genome_length:
                        index -= genome_length
        
                    protein = location[genome_number][index]
                    location_dictionary[k].append(inverted_clusters_checked[protein][0])
        return location_dictionary


    keep_array2 = []
    index = 0
    for b in tqdm(keep_array):
        a = []
        for i in flatten(clusters_checked[b]):
            a.append((genomes[i]))
        a2 = [item for sublist in a for item in sublist]
        values = list(collections.Counter(a2).values())
        test = np.array(values)>1
        if sum(test)>0 and len(values)>1: #no paralog splitting for homolog group with genes just from the same genome  
            test_ind = [test][0]
            keys = np.array(list(collections.Counter(a2).keys()))[test_ind]
            del(values)
            del(test)  
            ini.append(b)
            #terms = Paralog_list_create(a2, b) #array of genome number and location of homolog within that genome (# homolog by 2)
            #term_loc = Paralog_list_create(b)
            term_loc = collections.defaultdict(list)
            for j in flatten(clusters_checked[b]):
                for k in genomes[j]:
                    if len(hash_paralog[k][j])>1:
                        if len(term_loc[k])>0:
                            if hash_paralog[k][j] in term_loc[k]:
                                pass
                            else:
                                term_loc[k].append(hash_paralog[k][j])
                        else:
                            term_loc[k].append(hash_paralog[k][j])
                    else:
                        term_loc[k].append(hash_paralog[k][j])
            save = []
            for j in keys:
                terms = flatten(term_loc[j])
                test = np.zeros((len(terms),2))
                for k in np.arange(len(terms)):
                    test[k,0] = j
                test[:,1] = terms
                location_dict = Paralog_hunter(test) #get their neighborhood
                s_interpretation = set(location_dict[0]) #take first term
                for n in np.arange(1,len(test)): # for the rest of the homolog neighborhoods
                    s_interpretation = s_interpretation & set(location_dict[n]) # if there is at least some difference in neighborhood, then this will tell us
                if len(s_interpretation) < 2*neighborhood_size+1: #threshold for similarity of paralog neighborhoods, currently flexibility of 1     
                    save.append(j)
                else:
                    continue
            

            if len(save)>0:
                index2 = 0
                paralog_array = np.zeros((3,len(a2))).astype(int)
                index = 0
                for j in term_loc.keys():
                    add = flatten(term_loc[j])
                    paralog_array[0,index:index+len(add)] = j
                    paralog_array[1,index:index+len(add)] = add
                    index += len(add)
                for j in np.arange(len(paralog_array[0,:])):
                    paralog_array[2,j] = location[paralog_array[0,j]][paralog_array[1,j]]
                
                
                pa_dict = Paralog_hunter(paralog_array[0:2,:].T)
                #neighbor_array = np.zeros((11,len(a_hash_2))).astype(int)
                #for i in pa_dict.keys():
                #    neighbor_array[:,i] = pa_dict[i] 
                
                repeat = collections.defaultdict(list)
                i = save[0]
                temp_paralog_dict = collections.defaultdict(list)
                special_paralog_array = np.zeros((2,len(flatten(term_loc[i]))))
                special_paralog_array[0,:] = i
                special_paralog_array[1,:] = flatten(term_loc[i])
                spa_dict = Paralog_hunter(special_paralog_array.T)
                term_exclude = set(spa_dict[0])
                #term_all = set(spa_dict[0])
                #TEST_time2 = time.time() - TEST_start
                for j in spa_dict.keys():
                    term_exclude = term_exclude & set(spa_dict[j])
                    #term_all = term_all.union(set(spa_dict[j]))
                for j in spa_dict.keys():
                    spa_dict[j] = list(set(spa_dict[j]).difference(term_exclude))
                for j in pa_dict.keys():
                    temp_score = np.zeros((len(spa_dict.keys())))
                    for k in np.arange(len(spa_dict.keys())):
                        temp_score[k] = len(set(pa_dict[j]) & set(spa_dict[list(spa_dict.keys())[k]]))
                    temp_paralog_dict[index2 + temp_score.argmax()].append(j)
                index2 += len(spa_dict.keys())
                what_to_analyze = collections.defaultdict(list)
                #wta_analysis_info = collections.defaultdict(dict)
                wta = []
                final = []
                for j in temp_paralog_dict.keys():
                    multi_counts = [k for k, v in collections.Counter(paralog_array[0,temp_paralog_dict[j]]).items() if v > 1]
                    terms = list(set(save) & set(multi_counts))
                    if len(list(collections.Counter(paralog_array[0,temp_paralog_dict[j]]).keys()))==1:
                        final.append(j)
                    elif terms == list([save[0]]):
                        final.append(j)
                    elif len(terms)>0:
                        wta.append(j)
                        what_to_analyze[j].append(terms)
                        #for k in terms:
                        #    wta_analysis_info[j][k] = 
                    else:
                        final.append(j)
                if len(wta)>0:
                    while len(wta)>0: 
                    #for i in save:
                        for i in wta:
                            #for j in temp_paralog_dict[i]:
                                
                            save_term = what_to_analyze[i][0]
                            analysis_term = save_term[0]
                            
                            special_paralog_array = paralog_array[0:2,temp_paralog_dict[i]][:,paralog_array[0,temp_paralog_dict[i]]==analysis_term]
                            spa_dict = Paralog_hunter(special_paralog_array.T)
                            term_exclude = set(spa_dict[0])
                            #term_all = set(spa_dict[0])
                            #TEST_time2 = time.time() - TEST_start
                            for j in spa_dict.keys():
                                term_exclude = term_exclude & set(spa_dict[j])
                                #term_all = term_all.union(set(spa_dict[j]))
                            for j in spa_dict.keys():
                                spa_dict[j] = list(set(spa_dict[j]).difference(term_exclude))
                            
                            analysis_terms = collections.Counter()
                            for j in temp_paralog_dict[i]:
                                temp_score = np.zeros((len(spa_dict.keys())))
                                for k in np.arange(len(spa_dict.keys())):
                                    temp_score[k] = len(set(pa_dict[j]) & set(spa_dict[list(spa_dict.keys())[k]]))
                                temp_paralog_dict[index2 + temp_score.argmax()].append(j)
                                analysis_terms[index2 + temp_score.argmax()]+=1
                            
                            index2 += len(spa_dict.keys())
                            
                            analysis_keys = list(analysis_terms.keys())
                            for j in analysis_keys:
                                multi_counts = [k for k, v in collections.Counter(paralog_array[0,temp_paralog_dict[j]]).items() if v > 1]
                                terms = list(set(save) & set(multi_counts))
                                if temp_paralog_dict[i] == temp_paralog_dict[j]:
                                    repeat[j].append(analysis_term)
                                    terms = list(set(terms).difference({analysis_term}))
                                if len(repeat[i])>0:
                                    terms = list(set(terms).difference(set(flatten(repeat[i]))))
                                    repeat[j].append(flatten(repeat[i]))
                                if len(terms)>0:
                                    wta.append(j)
                                    what_to_analyze[j].append(terms)
                                else:
                                    final.append(j)      
                            del what_to_analyze[i]
                            del temp_paralog_dict[i]
                            wta = wta[1:len(wta)]
                else:
                    pass
                for i in final:
                    for j in temp_paralog_dict[i]:
                        split_clusters[index3].append(np.array([paralog_array[0,j],paralog_array[2,j]]))
                        identity_clusters[index3].append(identity_per_genome[paralog_array[0,j]][paralog_array[2,j]])
                    index3 +=1
            else:
                keep_array2.append(b)
                test_ind_time +=1
        else:
            keep_array2.append(b)
    
    print("Paralogs Separated:--- %s seconds ---" % "{:.2f}".format(time.time() - start_time))
    print(str(len(keep_array) - len(keep_array2)) +" of "+str(len(keep_array))+" original clusters split to form " + str(len(identity_clusters))+  " additional clusters.")
    print(" ")
    print(" ")
    
    
    
    with open(os.path.join(RESULTS_FOLDER,"Clusters_ID.tab"), 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter='\t')
        c = 0
        for b in np.arange(len(keep_array2)):
            a = []
            for i in flatten(clusters_checked[keep_array2[b]]):
                a.append((identity[i]))
            a = [item for sublist in a for item in sublist]
            c+=1
            tsv_writer.writerow(a)
        for c in identity_clusters.keys():
            a = identity_clusters[c]
            tsv_writer.writerow(a)
    length_ic = len(identity_clusters) + len(keep_array2)
    del(identity_clusters)
    print("Saved Clusters File:--- %s seconds ---" % (time.time() - start_time)) 
    print("Created " +str(length_ic)+" Clusters")
    print(" ")
    print(" ") 
            
###TO LOWER MEMORY FOOTPRINT
if dsp == False:
    hash_paralog.clear()
    identity_per_genome.clear()
    location.clear()
    inverted_clusters_checked.clear()
    genome_lengths = []

identity.clear() 

###################################################################
####################     DATA OUTPUT       ########################
###################################################################
#%%

if dsp == True:
    if ARGS.binary == True:
        index = 0
        Homolog_Split = np.zeros((num_genome,c))
        Homolog_metadata = np.zeros((3,c)).astype(object)
        Homolog_seq = np.zeros((1,c)).astype(object)
        for b in np.arange(len(keep_array)):
            a = []
            d = []
            g = []
            for i in flatten(clusters_checked[keep_array[b]]):
                a.append((genomes[i]))
                d.append((string[i]))
                if len(g) == 0:
                    g.append(seqs[i])
            a = [item for sublist in a for item in sublist]
            d = [item for sublist in d for item in sublist]
            for i in a:
                Homolog_Split[i,index]=1
            ue,count_terms = np.unique(d, return_counts=True)
            sorted_b = ue[np.argsort(count_terms)]
            if len(sorted_b)>3:
                Homolog_metadata[:,index] = sorted_b[0:3]
            else:
                if len(sorted_b)==2:
                    Homolog_metadata[0:2,index] = sorted_b[0:2]
                else:
                    Homolog_metadata[0:1,index] = sorted_b[0:1]
            Homolog_seq[0,index] = g[0][0]
            index+=1
        hcount_zero = -1*np.sum(Homolog_Split,axis=0)
        HG_split_ordered = Homolog_Split[:,np.argsort(hcount_zero)]
        HG_metadata = Homolog_metadata[:,np.argsort(hcount_zero)] 
        Homolog_seq = Homolog_seq[:,np.argsort(hcount_zero)]

    if ARGS.binary == False:
        index = 0
        Homolog_Split = np.zeros((num_genome,c))
        Homolog_metadata = np.zeros((3,c)).astype(object)
        Homolog_seq = np.zeros((1,c)).astype(object)
        Homolog_Split = Homolog_Split*np.nan
        for b in np.arange(len(keep_array)):
            a = []
            d = []
            e = []
            f = []
            g = []
            for i in flatten(clusters_checked[keep_array[b]]):
                a.append((genomes[i]))
                f.append(len(genomes[i]))
                d.append((string[i]))
                e.append((vrcs[i]))
                if len(g) == 0:
                    g.append(seqs[i])
            d = [item for sublist in d for item in sublist]
            e = np.array(e)[:,0,:]
            f = f.index(max(f))
            e = np.sum(abs(e[:,0:20]-e[f,0:20]),axis=1)
            for i in np.arange(len(a)):
                for j in a[i]:
                    Homolog_Split[j,index]=e[i]
            ue,count_terms = np.unique(d, return_counts=True)
            sorted_b = ue[np.argsort(count_terms)]
            if len(sorted_b)>3:
                Homolog_metadata[:,index] = sorted_b[0:3]
            else:
                if len(sorted_b)==2:
                    Homolog_metadata[0:2,index] = sorted_b[0:2]
                else:
                    Homolog_metadata[0:1,index] = sorted_b[0:1]
            Homolog_seq[0,index] = g[0][0]
            index+=1
        hcount_zero = np.count_nonzero(np.isnan(Homolog_Split),axis=0)
        HG_split_ordered = Homolog_Split[:,np.argsort(hcount_zero)]
        HG_metadata = Homolog_metadata[:,np.argsort(hcount_zero)] 
        Homolog_seq = Homolog_seq[:,np.argsort(hcount_zero)]

else:
    if ARGS.binary == True:
        index = 0
        Homolog_Split = np.zeros((num_genome,len(keep_array2) + len(split_clusters)))
        Homolog_metadata = np.zeros((3,len(keep_array2) + len(split_clusters))).astype(object)
        Homolog_seq = np.zeros((1,len(keep_array2) + len(split_clusters))).astype(object)
        for b in np.arange(len(keep_array2)):
            a = []
            d = []
            g = []
            for i in flatten(clusters_checked[keep_array[b]]):
                a.append((genomes[i]))
                d.append((string[i]))
                if len(g) == 0:
                    g.append(seqs[i])
            a = [item for sublist in a for item in sublist]
            d = [item for sublist in d for item in sublist]
            for i in a:
                Homolog_Split[i,index]=1
            ue,count_terms = np.unique(d, return_counts=True)
            sorted_b = ue[np.argsort(count_terms)]
            if len(sorted_b)>3:
                Homolog_metadata[:,index] = sorted_b[0:3]
            else:
                if len(sorted_b)==2:
                    Homolog_metadata[0:2,index] = sorted_b[0:2]
                else:
                    Homolog_metadata[0:1,index] = sorted_b[0:1]
            Homolog_seq[0,index] = g[0][0]
            index+=1
        for c in split_clusters.keys():
            a = []
            b = []
            d = []
            g = []
            for i in np.arange(len(split_clusters[c])):
                a.append(split_clusters[c][i][0])
                b.append(split_clusters[c][i][1])
            pre_d = np.unique(b, return_counts=True)[0]
            for i in np.arange(len(pre_d)):
                d.append(string[pre_d[i]])
                if len(g) == 0:
                    g.append(seqs[pre_d[i]])
            d = [item for sublist in d for item in sublist]
            for i in a:
                Homolog_Split[i,index]=1
            ue,count_terms = np.unique(d, return_counts=True)
            sorted_b = ue[np.argsort(count_terms)]
            
            if len(sorted_b)>3:
                Homolog_metadata[:,index] = sorted_b[0:3]
            else:
                if len(sorted_b)==2:
                    Homolog_metadata[0:2,index] = sorted_b[0:2]
                else:
                    Homolog_metadata[0:1,index] = sorted_b[0:1]
            Homolog_seq[0,index] = g[0][0]
            index+=1 
        hcount_zero = -1*np.sum(Homolog_Split,axis=0)
        HG_split_ordered = Homolog_Split[:,np.argsort(hcount_zero)]
        HG_metadata = Homolog_metadata[:,np.argsort(hcount_zero)]
        Homolog_seq = Homolog_seq[:,np.argsort(hcount_zero)]

    if ARGS.binary == False:
        index = 0
        Homolog_Split = np.zeros((num_genome,len(keep_array2) + len(split_clusters)))
        Homolog_metadata = np.zeros((3,len(keep_array2) + len(split_clusters))).astype(object)
        Homolog_seq = np.zeros((1,len(keep_array2) + len(split_clusters))).astype(object)
        Homolog_Split = Homolog_Split*np.nan
        for b in tqdm(np.arange(len(keep_array2))):
            a = []
            d = []
            e = []
            f = []
            g = []
            for i in flatten(clusters_checked[keep_array[b]]):
                a.append((genomes[i]))
                f.append(len(genomes[i]))
                d.append((string[i]))
                e.append((vrcs[i]))
                if len(g) == 0:
                    g.append(seqs[i])
            d = [item for sublist in d for item in sublist]
            e = np.array(e)[:,0,:]
            f = f.index(max(f))
            e = np.sum(abs(e[:,0:20]-e[f,0:20]),axis=1)#/sum(e[f,0:20]) #difference versus proportion
            for i in np.arange(len(a)):
                for j in a[i]:
                    Homolog_Split[j,index]=e[i]
            ue,count_terms = np.unique(d, return_counts=True)          
            sorted_b = ue[np.argsort(count_terms)]
            if len(sorted_b)>3:
                Homolog_metadata[:,index] = sorted_b[0:3]
            else:
                if len(sorted_b)==2:
                    Homolog_metadata[0:2,index] = sorted_b[0:2]
                else:
                    Homolog_metadata[0:1,index] = sorted_b[0:1]
            Homolog_seq[0,index] = g[0][0]
            index+=1
        for c in tqdm(split_clusters.keys()):
            a = []
            b = []
            d = []
            e = []
            g = []
            for i in np.arange(len(split_clusters[c])):
                a.append(split_clusters[c][i][0])
                b.append(split_clusters[c][i][1])
            pre_e,f = np.unique(b, return_counts=True)
            for i in np.arange(len(pre_e)):
                e.append((vrcs[pre_e[i]]))
                d.append(string[pre_e[i]])
                if len(g) == 0:
                    g.append(seqs[pre_e[i]])
            e = np.array(e)[:,0,:]
            d = [item for sublist in d for item in sublist]
            f = list(f).index(max(f))
            e = np.sum(abs(e[:,0:20]-e[f,0:20]),axis=1)
            temp_dict = collections.defaultdict(list)
            for i in np.arange(len(pre_e)):
                temp_dict[pre_e[i]].append(e[i])
            for i in np.arange(len(split_clusters[c])):
                Homolog_Split[split_clusters[c][i][0],index]=temp_dict[split_clusters[c][i][1]][0]
            ue,count_terms = np.unique(d, return_counts=True)
            sorted_b = ue[np.argsort(count_terms)]
            if len(sorted_b)>3:
                Homolog_metadata[:,index] = sorted_b[0:3]
            else:
                if len(sorted_b)==2:
                    Homolog_metadata[0:2,index] = sorted_b[0:2]
                else:
                    Homolog_metadata[0:1,index] = sorted_b[0:1]
            Homolog_seq[0,index] = g[0][0]
            index+=1
        hcount_zero = np.count_nonzero(np.isnan(Homolog_Split),axis=0)
        HG_split_ordered = Homolog_Split[:,np.argsort(hcount_zero)]
        HG_metadata = Homolog_metadata[:,np.argsort(hcount_zero)] 
        Homolog_seq = Homolog_seq[:,np.argsort(hcount_zero)]

###TO LOWER MEMORY FOOTPRINT
string.clear()
seqs.clear()
vrcs.clear()
del(keep_array)
        
Folder_species_name_list = []
for i in np.arange(len(Folder_list)):
    Folder_species_name_list.append(Folder_list[i].split("/")[-1])

if ARGS.large_dataset == False:
    print("Starting to save File: --- %s seconds ---" % (time.time() - start_time)) #TIME
    df_met = pd.DataFrame(HG_metadata.T)
    df_met = df_met.rename(columns={0: "gene name", 1: "alternative gene name 1",2:"alternative gene name 2"})
    df_seq = pd.DataFrame(Homolog_seq.T)
    df_seq = df_seq.rename(columns={0: "Sample sequence"})
    a2 = Folder_species_name_list
    df = pd.DataFrame(HG_split_ordered.T, columns = a2)
    df = df.fillna('')
    df_concat = pd.concat([df_met, df_seq], axis=1)
    df_concat = pd.concat([df_concat, df], axis=1)
    df_concat.to_csv(os.path.join(RESULTS_FOLDER,"Pangenome_Numerical.csv"), index=False)
    print("Results Saved:--- %s seconds ---" % (time.time() - start_time))

          
if ARGS.large_dataset == True:
    print("Starting to save File: --- %s seconds ---" % (time.time() - start_time)) #TIME
    df_met = pd.DataFrame(HG_metadata.T)
    df_met = df_met.rename(columns={0: "gene name", 1: "alternative gene name 1",2:"alternative gene name 2"})
    df_seq = pd.DataFrame(Homolog_seq.T)
    df_seq = df_seq.rename(columns={0: "Sample sequence"})
    a2 = Folder_species_name_list
    df = pd.DataFrame(HG_split_ordered.T, columns = a2)
    #df = df.fillna("")
    df_concat = pd.concat([df_met, df_seq], axis=1)
    df_concat = pd.concat([df_concat, df], axis=1)
    df_concat.astype(str).to_feather(os.path.join(RESULTS_FOLDER,"Pangenome_Numerical.feather"))
    print("Results Saved:--- %s seconds ---" % (time.time() - start_time)) 


if ARGS.graphical_mode == True:
    pass
else:
    HG_split_ordered[~np.isnan(HG_split_ordered)] = 1
    HG_split_ordered[np.isnan(HG_split_ordered)] = 0
    hist_vals = np.sum(HG_split_ordered,axis=0)
    core_acc, bins = np.histogram(hist_vals, bins = np.array([0,0.15*num_genome,0.95*num_genome,0.99*num_genome,max(hist_vals)+1]))
    print("")
    print("")
    print('ARPA found:')
    print(str(core_acc[-1]) + ' CORE genes      (99% -- 100%)')
    print(str(core_acc[-2]) + ' SOFT CORE genes (95% -- 100%)')
    print(str(core_acc[-3]) + ' SHELL genes     (15% -- 95%)')
    print(str(core_acc[-4]) + ' CLOUD genes     (0% -- 15%)')
    print('-----------------------')
    print(str(sum(core_acc)) + ' TOTAL genes     (0% -- 100%)')
    print(" ")
    print(" ")
    print(" ")
    print("Done! Thank you for using ARPA.")
    
    sys.exit()

###################################################################
#####################     GRAPHICS       ##########################
###################################################################
#%%


### TO FIX
### ADD GRAPHICS HERE TO THE BINARY ONE TOO
###


if ARGS.binary == True:
    print("Starting to save figure: --- %s seconds ---" % (time.time() - start_time)) #TIME
    plt.imsave(os.path.join(RESULTS_FOLDER,'Pangenome_Visual.png'), HG_split_ordered, cmap='turbo', vmin=0, vmax=1)
    figure(figsize=(120, 60), dpi=80)
    print("Starting to save figure 2: --- %s seconds ---" % (time.time() - start_time)) #TIME
    plt.imshow(HG_split_ordered, cmap='plasma', vmin=0, vmax=30,interpolation='nearest', aspect='auto')
    plt.xlabel('Gene Clusters')
    plt.ylabel('Genomes')
    plt.gca().axes.yaxis.set_ticklabels([])
    plt.savefig(os.path.join(RESULTS_FOLDER,'Pangenome_F1.png'),dpi=400)
    print("Finished Saving Figures: --- %s seconds ---" % (time.time() - start_time)) #TIME


          
if ARGS.binary == False:
    indexes = np.sum(~np.isnan(HG_split_ordered),axis=0)
    for i in np.arange(len(indexes)):
        if indexes[i] < 0.95*num_genome:
            start_acc = i
            break
    for i in np.arange(len(indexes)):
        if indexes[i] < 0.15*num_genome:
            end_acc = i
            break
        else:
            end_acc = indexes[-1]
  
    Z = ward(pdist(np.isnan(HG_split_ordered)))
    list_out = leaves_list(Z)
    a = os.listdir(Folder)
    a2 = [0] * num_genome
    HG_split_ordered_DN = HG_split_ordered*0
    for i in np.arange(len(list_out)):
        HG_split_ordered_DN[i,:] = HG_split_ordered[list_out[i],:]
        a2[i] = a[list_out[i]]
   
    ## OPTION: IF YOU WANT TO CLUSTER WITH FASTTREE INSTEAD OF DENDOGRAM
   
    #Accessory_HG_split_phylogenetics = pd.DataFrame(HG_split_ordered[:,start_acc:end_acc].T)
    #Accessory_HG_split_phylogenetics[~Accessory_HG_split_phylogenetics.isnull()] = 1
    #Accessory_HG_split_phylogenetics[Accessory_HG_split_phylogenetics.isnull()] = 0
    #file = open(os.path.join(RESULTS_FOLDER,"accessory_presence_absence.phy"), "w")
    #for i in np.arange(num_genome):
    #    file.write(">"+str(Folder_species_name_list[i])+ "\n")
    #    file.write(str("".join(map(str,Accessory_HG_split_phylogenetics[i].to_list())))+ "\n")
    #file.close()
    #print("Running Fast Tree on Accessory genes to sort Pangenome: --- %s seconds ---" % (time.time() - start_time)) #TIME
    #input_file = str(os.path.join(RESULTS_FOLDER,"accessory_presence_absence.phy"))
    #output_file = str(os.path.join(RESULTS_FOLDER,"accessory_alignment_output"))
    #os.system('fasttree -nt '+ str(input_file) +' > '+ str(output_file))
    #from Bio import Phylo
    #tree = Phylo.read(str(output_file),"newick")
    #location_name = collections.defaultdict(list)
    #for i in np.arange(len(Folder_species_name_list)):
    #    location_name[Folder_species_name_list[i]].append(i)
    #HG_split_ordered_DN = np.zeros((len(HG_split_ordered),len(HG_split_ordered[0,:])))
    #index = 0
    #for i in tree.get_terminals(): 
    #    HG_split_ordered_DN[index,:] = HG_split_ordered[location_name[i.name],:]
    #    index +=1
    
    print("Starting to save Pangenome figure: --- %s seconds ---" % (time.time() - start_time)) #TIME
    #plt.imsave(os.path.join(RESULTS_FOLDER,'Pangenome_Visual.png'), HG_split_ordered_DN, cmap='turbo', vmin=0, vmax=1)
    figure(figsize=(120, 60), dpi=80)
    #print("Starting to save figure 2: --- %s seconds ---" % (time.time() - start_time)) #TIME
    plt.imshow(HG_split_ordered_DN, cmap='turbo', vmin=0, vmax=20,interpolation='nearest', aspect='auto')
    plt.xlabel('Gene Clusters')
    plt.ylabel('Genomes')
    plt.gca().axes.yaxis.set_ticklabels([])
    plt.savefig(os.path.join(RESULTS_FOLDER,'Pangenome_F1.png'),dpi=80)
    print("Finished Saving Pangenome Figures: --- %s seconds ---" % (time.time() - start_time)) #TIME

import seaborn as sns

#Pangenome Collector's Curve

if ARGS.binary == False:
    #plt.figure()
    HG_split_ordered2 = copy.deepcopy(HG_split_ordered)
    HG_split_ordered2[~np.isnan(HG_split_ordered)] = 1
    HG_split_ordered2[np.isnan(HG_split_ordered)] = 0

else:
    HG_split_ordered2 = HG_split_ordered

pangenome_size = np.zeros((4,num_genome))
pangenome_size_conserv = np.zeros((4,num_genome))

upper_lim = np.array([1,0.99,0.98,0.9])
lower_lim = np.array([1,2,3,4])
for i in tqdm(np.arange(1, num_genome+1)):
    counts = np.sum(HG_split_ordered2[0:i,:],axis=0)
    c = collections.Counter(counts)
    for j in np.arange(len(upper_lim)):
        pangenome_size_conserv[j,i-1] = sum(counts > (upper_lim[j]*i)-0.01)
    for j in np.arange(len(lower_lim)):
        pangenome_size[j,i-1] = sum(counts > (lower_lim[j])-0.01)
figure(figsize=(100, 50), dpi=80)        
sns.set_theme(style="whitegrid")
sns.despine()
colors_top = ['maroon','firebrick','indianred','lightcoral']
colors_bottom = ['darkblue','blue','royalblue','dodgerblue']

for i in np.arange(len(lower_lim)):
    plt.plot(np.arange(i,num_genome), pangenome_size.T[i:num_genome+1,i],lw=9-i,c=colors_top[i], label='No. genes present in > ' +str(lower_lim[i])+' genomes')

for i in np.arange(len(upper_lim)):   
    plt.plot(np.arange(num_genome),pangenome_size_conserv.T[:,i], lw = 9-i, c=colors_bottom[i], label='No. genes conserved in ' +str(int(100*upper_lim[i]))+'% of genomes')
plt.xlim(-(num_genome/20),num_genome+(num_genome/60))
sze = min(max(4,num_genome/10),70)

plt.legend(fontsize=sze)
plt.xticks(fontsize = sze)
plt.yticks(fontsize = sze)
plt.savefig(os.path.join(RESULTS_FOLDER,'Pangenome_Collectors_Curve.png'),dpi=100)


    



#Histogram of clusters
plt.figure()
hist_vals = np.sum(HG_split_ordered2,axis=0)
hist_term, bins = np.histogram(hist_vals, bins = np.arange(num_genome+1)-0.01)
hist_term[-1]
sns.set_theme(style="whitegrid", palette="pastel")

sns.histplot(x=hist_vals, color='Darkblue',bins = np.arange(max(hist_vals)+1)-0.01)
plt.ylim(0,hist_term[-1]+50)
plt.xlabel('Number of Genomes in gene cluster')
plt.savefig(os.path.join(RESULTS_FOLDER,'Group_Histogram.png'),dpi=150)

hist_vals = np.sum(HG_split_ordered2,axis=0)
core_acc, bins = np.histogram(hist_vals, bins = np.array([0,0.15*num_genome,0.95*num_genome,0.99*num_genome,max(hist_vals)+1]))


print("")
print("")
print('ARPA found:')
print(str(core_acc[-1]) + ' CORE genes      (99% -- 100%)')
print(str(core_acc[-2]) + ' SOFT CORE genes (95% -- 100%)')
print(str(core_acc[-3]) + ' SHELL genes     (15% -- 95%)')
print(str(core_acc[-4]) + ' CLOUD genes     (0% -- 15%)')
print('-----------------------')
print(str(sum(core_acc)) + ' TOTAL genes     (0% -- 100%)')




sys.exit()
