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
from Bio import SeqIO
from tqdm import tqdm
from itertools import compress

###################################################################

start_time = time.time()
PARSER = argparse.ArgumentParser(
    prog="ARPA_main_clust.py",
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
    print("ARPA_main_clust.py -ct 0.95 /user/path/to/faa/files/")
    print("")
    print("")
    print(" --------------- ")
    print("")
    sys.exit(0)
    
ARGS = PARSER.parse_args()

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
        os.mkdir(str(ARGS.output_file) + "ARPA_Align_results")
        RESULTS_FOLDER = str(ARGS.output_file) + "ARPA_Align_results{}".format(OS_SEPARATOR)
    else:
        os.mkdir(str(ARGS.output_file) + str(OS_SEPARATOR) + "ARPA_Align_results")
        RESULTS_FOLDER = str(ARGS.output_file) + str(OS_SEPARATOR) + "ARPA_Align_results{}".format(OS_SEPARATOR)
else:
    os.mkdir("ARPA_Align_results") #make a new directory in the working directory to store analysis files
    RESULTS_FOLDER = "ARPA_Align_results{}".format(OS_SEPARATOR) #results folder uses he current timestring to create a new file that is time-stamped
    
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

num_genome = len(Folder_list) #number of genomes for analysis

def FAAtoNUM(i): #define this function that adds one genome (FAA) to the collections default dicts defined above.
    pathToFile = open(Folder_list[i]) #import file with the open command
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

sys.exit()

