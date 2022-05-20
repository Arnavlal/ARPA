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
from matplotlib.pyplot import figure
from scipy.cluster.hierarchy import ward, leaves_list
from scipy.spatial.distance import pdist
from sklearn.cluster import DBSCAN


#OS Work
###################################################################
start_time = time.time()
PARSER = argparse.ArgumentParser(
    prog="ARPA.py",
    description="Alignment-Free Residue Pangenomic Analysis",
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
    "-m",
    "--metric",
    type=str,
    help="metric for genes: 'mr' diff. From most represented, 'aa' diff. from all alleles, 'aan' diff. from all alleles normalized by counts [Default: 'mr']",
)

PARSER.add_argument(
    "-dsp",
    "--dont_split_paralogs",
    action='store_true',
    help="Do not split paralogous groups (no neighborhood analysis)",
)

PARSER.add_argument(
    "-ld",
    "--large_dataset",
    action='store_true',
    help="Dendogram sorting step avoided, full image not printed, and pangenome file saved as Pickle file, not CSV",
)

PARSER.add_argument(
    "-ni",
    "--neighborhood_identity",
    type=float,
    help="threshold of genes to be reclustered in neighborhood analysis [Default: 0.7]",
)

PARSER.add_argument(
    "-nt",
    "--neighborhood_thresh",
    type=float,
    help="threshold of genes conserved within neighborhood (out of 10) [Default: 3]",
)

PARSER.add_argument(
    "sequences_folder", type=str, help="directory with .faa files for pangenome"
)
if len(sys.argv) == 1:
    PARSER.print_help()
    sys.exit(0)
ARGS = PARSER.parse_args()

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
if ARGS.metric:
    if ARGS.metric != "mr" and ARGS.metric != "aa" and ARGS.metric != "aan":
        PARSER.exit(status=0, message="Error: unrecognized string inputted into metric analysis\n")


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
#FIX THIS
#ONLY LOOK AT FAA FILES NOT .DS_STORE   
#!!!
#!!!
#!!!   
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

del(value_of_genome_breakpoints)
del(placeholder)
del(AA)

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
count = collections.Counter(hash_genome)
for i in np.arange(len(AApangenome)):
    AApangenome_count[i,0] = count[hash_pangenome[i]]
AApangenome  = np.concatenate((AApangenome,AApangenome_hash), axis=1)
AApangenome  = np.concatenate((AApangenome,np.zeros((len(AApangenome_hash),1))), axis=1)
AApangenome  = np.concatenate((AApangenome,AApangenome_count), axis=1)


del(AApangenome_count)
del(hash_genome)
del(hash_genome_translation)
del(hash_genome_set)
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
                    if first_point == 0:
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


diff  = np.diff(new_arr[:,21])
breakpoints = np.where(diff==1)[0]
most_represented = np.zeros((int(new_arr[-1,21]), 24))
mval = int(new_arr[-1,21])
for i in np.arange(new_arr[-1,21]):
    if i ==0:
        genome = new_arr[0:breakpoints[0]+1,:]
        maximum = np.where(genome[:,22] == np.max(genome[:,22]))[0][0]
        new_arr[0:breakpoints[0]+1,22] = np.sum(abs(genome[:,0:20] - genome[int(maximum),0:20]), axis=1)
    elif i ==new_arr[-1,21]-1:
        genome = new_arr[breakpoints[int(i-1)]+1:int(len(new_arr)),:]
        maximum = np.where(genome[:,22] == np.max(genome[:,22]))[0][0]
        new_arr[breakpoints[int(i-1)]+1:int(len(new_arr)),22] = np.sum(abs(genome[:,0:20] - genome[int(maximum),0:20]), axis=1)
    else:
        genome = new_arr[breakpoints[int(i-1)]+1:breakpoints[int(i)]+1,:]
        maximum = np.where(genome[:,22] == np.max(genome[:,22]))[0][0]
        new_arr[breakpoints[int(i-1)]+1:breakpoints[int(i)]+1,22] = np.sum(abs(genome[:,0:20] - genome[int(maximum),0:20]), axis=1)

print("Generated Allele Scores: --- %s seconds ---" % (time.time() - start_time)) #TIME

#FIX THIS
#MAKE AS FAST AS MR    
#!!!
#!!!
#!!!    
if ARGS.metric == "aa":
    for i in np.arange(new_arr[-1,21]):
        for j in np.arange(len(new_arr[np.where(new_arr[:,21]==most_represented[int(i),21])[0],22])):
            new_arr[min(np.where(new_arr[:,21]==most_represented[int(i),21])[0])+j,22] = np.sum(abs(new_arr[min(np.where(new_arr[:,21]==most_represented[int(i),21])[0])+j,0:20]-new_arr[np.where(new_arr[:,21]==most_represented[int(i),21])[0],0:20]))/len(new_arr[np.where(new_arr[:,21]==most_represented[int(i),21])[0],0:20])

if ARGS.metric == "aan":
    for i in np.arange(new_arr[-1,21]):
        for j in np.arange(len(new_arr[np.where(new_arr[:,21]==most_represented[int(i),21])[0],22])):
            new_arr[min(np.where(new_arr[:,21]==most_represented[int(i),21])[0])+j,22] = np.sum(abs(new_arr[min(np.where(new_arr[:,21]==most_represented[int(i),21])[0])+j,0:20]-new_arr[np.where(new_arr[:,21]==most_represented[int(i),21])[0],0:20]))/len(new_arr[np.where(new_arr[:,21]==most_represented[int(i),21])[0],0:20])/most_represented[int(i),0:20]



test_keys = list(new_arr[:,20])
test_values = list(new_arr[:,21])
test_values2 = list(new_arr[:,22])
res = {test_keys[i]: test_values[i] for i in range(len(test_keys))}
res2 = {test_keys[i]: test_values2[i] for i in range(len(test_keys))}

del(new_arr)

for i in np.arange(len(translation[0,:])):
    translation[2,i] = res[translation[0,i]]
    translation[3,i] = res2[translation[0,i]]

seq_info = translation[1:4,:].T 
seq_info = seq_info[seq_info[:, 0].argsort()]
placeholder = np.zeros((len(seq_info),1))
placeholder[:,0] = index
seq_w_genome = np.concatenate((seq_info, placeholder), axis=1)

del(res)
del(res2)
del(placeholder)
del(test_keys)
del(test_values)
del(test_values2)
del(seq_info)

if ARGS.dont_split_paralogs==True:
    if ARGS.large_dataset==True:
        print("Starting Addition of allele scores to File: --- %s seconds ---" % (time.time() - start_time)) #TIME
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
        a = os.listdir(Folder)
        print("Starting to save figure: --- %s seconds ---" % (time.time() - start_time)) #TIME
        plt.figure(figsize=(60, 50), dpi=80)
        plt.imshow(HG_split_ordered, cmap='turbo', vmin=0, vmax=30,interpolation='nearest', aspect='auto')
        plt.colorbar()
        plt.savefig(os.path.join(RESULTS_FOLDER,'Wholistic_Pangenome.png'))
        print("Starting to save File: --- %s seconds ---" % (time.time() - start_time)) #TIME
        df_met = pd.DataFrame(HG_metadata.T)
        df_met = df_met.rename(columns={0: "gene name", 1: "alternative gene name 1",2:"alternative gene name 2"})
        df = pd.DataFrame(HG_split_ordered.T, columns = a)
        df = df.fillna('')
        df.to_pickle(os.path.join(RESULTS_FOLDER,"Pangenome_Numerical.pickle"))
        print("Results:--- %s seconds ---" % (time.time() - start_time)) 
        print("[==================>]")
        print(" ")
        print(" ")
        print("Done! Thank you for using ARPA.")
        sys.exit()
    else:
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
        Z = ward(pdist(np.isnan(HG_split_ordered)))
        list_out = leaves_list(Z)
        #dn = hierarchy.dendrogram(hierarchy.linkage(np.isnan(HG_split_ordered), 'single'), labels=os.listdir(Folder))
        a = os.listdir(Folder)
        a2 = [0] * num_genome
        print("Dendogram analysis done: --- %s seconds ---" % (time.time() - start_time)) #TIME
        plt.close()
        HG_split_ordered_DN = HG_split_ordered*0
        for i in np.arange(len(list_out)):
            HG_split_ordered_DN[i,:] = HG_split_ordered[list_out[i],:]
            a2[i] = a[list_out[i]]
        print("Starting to save figures: --- %s seconds ---" % (time.time() - start_time)) #TIME
        plt.imsave(os.path.join(RESULTS_FOLDER,'Pangenome_Visual.png'), HG_split_ordered_DN, cmap='turbo', vmin=0, vmax=30)
        figure(figsize=(120, 60), dpi=80)
        print("Starting to save figure 2: --- %s seconds ---" % (time.time() - start_time)) #TIME
        plt.imshow(HG_split_ordered_DN, cmap='turbo', vmin=0, vmax=30,interpolation='nearest', aspect='auto')
        plt.colorbar()
        plt.savefig(os.path.join(RESULTS_FOLDER,'Wholistic_Pangenome.png'))
        print("Starting to save File: --- %s seconds ---" % (time.time() - start_time)) #TIME
        df_met = pd.DataFrame(HG_metadata.T)
        df_met = df_met.rename(columns={0: "gene name", 1: "alternative gene name 1",2:"alternative gene name 2"})
        df = pd.DataFrame(HG_split_ordered_DN.T, columns = a2)
        df = df.fillna('')
        df.to_csv(os.path.join(RESULTS_FOLDER,"Pangenome_Numerical.csv"), index=False)
        df_concat = pd.concat([df_met, df], axis=1)
        df_concat.to_csv(os.path.join(RESULTS_FOLDER,"Pangenome_Numerical.csv"), index=False)
        print("Results:--- %s seconds ---" % (time.time() - start_time)) 
        print("[==================>]")
        print(" ")
        print(" ")
        print("Done! Thank you for using ARPA.")
        sys.exit()


#%%
#%%            


#LOOK AT THE GENE NEIGHBORHOOD
#essentially take the concatenated version of the list of genes within the genome
#and then look at the specific genes (their orthologous group) within the genome.
#For one OG:
def PG_finder(OG):
    a  = np.zeros((len(np.where(seq_w_genome[:,1]==OG)[0]),16))
    for i in np.arange(len(a)):
        a[i,0:4] = seq_w_genome[np.where(seq_w_genome[:,1]==OG)[0][i],:]
        if a[i,0] < 6:
            num_ovr = 5-a[i,0]
            new_bp = min(np.where(seq_w_genome[:,3]==1)[0])-1
            a[i,4:4+int(num_ovr)] = seq_w_genome[int(new_bp)-int(num_ovr):int(new_bp),1]
            a[i,4+int(num_ovr):9] = seq_w_genome[0:int(a[i,0]),1]
            a[i,11:16] = seq_w_genome[int(a[i,0])+1:int(a[i,0])+6,1]
        elif a[i,0] > len(seq_w_genome)-6:
            num_at = len(seq_w_genome)-a[i,0]-1
            new_bp = min(np.where(seq_w_genome[:,3]==max(seq_w_genome[:,3]))[0])
            a[i,11+int(num_at):16] = seq_w_genome[int(new_bp):int(new_bp)+5-int(num_at),1]
            a[i,11:11+int(num_at)] = seq_w_genome[int(a[i,0])+1:int(a[i,0])+int(num_at)+1,1]
            a[i,4:9] = seq_w_genome[int(a[i,0])-5:int(a[i,0]),1]
        elif min(abs(a[i,0]-(cumulative_breakpoints)))<6:
            if a[i,0]<cumulative_breakpoints[0]:
                num_at = abs(cumulative_breakpoints[int(a[i,3])]-a[i,0])-1
                num_ovr = 5-num_at
                a[i,11+int(num_at):16] = seq_w_genome[0:0+int(num_ovr),1]
                a[i,11:11+int(num_at)] = seq_w_genome[int(a[i,0])+1:int(a[i,0])+1+int(num_at),1]
                a[i,4:9] = seq_w_genome[int(a[i,0])-5:int(a[i,0]),1]

            elif abs(a[i,0]-cumulative_breakpoints[int(a[i,3])])<6:#back of genome
                num_at = abs(cumulative_breakpoints[int(a[i,3])]-a[i,0])-1
                num_ovr = 5-num_at
                a[i,11+int(num_at):16] = seq_w_genome[int(cumulative_breakpoints[int(a[i,3]-1)]):int(cumulative_breakpoints[int(a[i,3]-1)])+int(num_ovr),1]
                a[i,11:11+int(num_at)] = seq_w_genome[int(a[i,0])+1:int(a[i,0])+1+int(num_at),1]
                a[i,4:9] = seq_w_genome[int(a[i,0])-5:int(a[i,0]),1]
            else:#front of genome
                if a[i,0]-cumulative_breakpoints[int(a[i,3]-1)]<6:
                    num_at = a[i,0]-cumulative_breakpoints[int(a[i,3])-1]
                    num_ovr = 5-num_at
                    a[i,4:9-int(num_at)] = seq_w_genome[int(cumulative_breakpoints[int(a[i,3])])-int(num_ovr):int(cumulative_breakpoints[int(a[i,3])]),1]
                    a[i,9-int(num_at):9] = seq_w_genome[int(cumulative_breakpoints[int(a[i,3])-1]):int(a[i,0]),1]
                    a[i,11:16] = seq_w_genome[int(a[i,0])+1:int(a[i,0])+6,1]
                else:
                   a[i,4:9] = seq_w_genome[int(a[i,0])-5:int(a[i,0]),1]
                   a[i,11:16] = seq_w_genome[int(a[i,0])+1:int(a[i,0])+6,1] 
        else:
            a[i,4:9] = seq_w_genome[int(a[i,0])-5:int(a[i,0]),1]
            a[i,11:16] = seq_w_genome[int(a[i,0])+1:int(a[i,0])+6,1]
    return a

#NEIGHBORHOOD SPLITTING OF PARALOGS

OGcount= 0
Orthologous_Split = np.zeros((num_genome,int(len(AApangenome_hash))-1))*np.nan
Ortholog_metadata = np.zeros((3,int(len(AApangenome_hash))-1)).astype(str)

for n in np.arange(1,int(max(translation[2,:]))):
    a = PG_finder(n)
    #condenser
    test_val,count = np.unique(np.vstack((a[:,4:9],a[:,11:16])), return_counts=True)
    count_sort_ind = np.argsort(-count)
    test_val = test_val[count_sort_ind]
    protein_test = np.zeros((len(test_val),24))
    protein_test = most_represented[(test_val -1).astype(int),:]
    q = 0.1
    for j in np.arange(len(protein_test)):
        if len(protein_test)>0.5:
            test = np.sum(abs(protein_test[:,0:20]-protein_test[0,0:20]),axis=1)
            thresh = (1-CGN_thresh)*np.sum(protein_test[0,0:20])
            point = np.where(test<thresh)[0]
            points = protein_test[point,23].astype(int)
            for k in np.arange(len(points)):
                a[:,4:9] = np.where(a[:,4:9]==points[k], q, a[:,4:9])
                a[:,11:16] = np.where(a[:,11:16]==points[k], q, a[:,11:16])
            protein_test = np.delete(protein_test, point, axis=0)
            q +=1
        else:
            continue
    a[:,4:9] = a[:,4:9]-0.1
    a[:,11:16] = a[:,11:16]-0.1
    #condenser end
    if len(a)>1:
        a_scores = np.zeros((len(a), len(a)))
        a_mod = np.delete(np.delete(a[:,4:16],5,1),5,1)
        #pairwise comp
        for i in np.arange(len(a)):
            test_group = a_mod[i,:]
            duplicate_test, counts = np.unique(test_group, return_counts=True)
            for j in np.arange(len(duplicate_test)):
                group_id = a[:,4:16]==duplicate_test[j]
                group_id = np.sum(group_id, axis=1)
                group_id[group_id>counts[j]] = counts[j]
                a_scores[:,i]+=group_id
        clust = DBSCAN(eps=0.3, min_samples=1, metric='euclidean', algorithm='auto').fit(a_scores)
        placeholder = np.zeros((len(a),1))
        placeholder[:,0] = clust.labels_
        together = np.concatenate((a[:,4:16],placeholder), axis=1)
        together = together[together[:, -1].argsort()]
        a_scores = np.zeros((len(a), len(a)))
        a_mod = together[:,0:12]
        a_mod = np.delete(np.delete(a_mod,5,1),5,1)
        #pairwise comp
        for i in np.arange(len(a_mod)):
            test_group = a_mod[i,:]
            duplicate_test, counts = np.unique(test_group, return_counts=True)
            for j in np.arange(len(duplicate_test)):
                group_id = a_mod==duplicate_test[j]
                group_id = np.sum(group_id, axis=1)
                group_id[group_id>counts[j]] = counts[j]
                a_scores[:,i]+=group_id
        index = np.arange(np.max(clust.labels_)+1)
        for i in np.arange(np.max(clust.labels_)):
            index[i+1] = np.max(np.where(together[:,12]==i)[0])
        index +=1
        index[0]=0
        index = np.append(index,np.array(len(a_scores)))
        a_dup = a_mod*0
        ind = np.arange(np.max(clust.labels_)+1)
        k=0
        cutoffs = []
        for i in np.arange(len(index)-1):
            if np.sum(ind==i)>0:
                rm = []
                for j in ind: #add the "-1"?
                    if len(ind)>0:
                        if np.average(a_scores[index[i]:index[i+1],index[j]:index[j+1]])>NTHRES:
                            a_dup[k:k+index[j+1]-index[j],:] = a_mod[index[j]:index[j+1],:]
                            k+=index[j+1]-index[j]
                            rm.append(np.where(ind==j)[0][0])
                            #print(rm)
                sel = []
                addnl_terms = len(rm)
                while addnl_terms>0:
                    mr = rm
                    for l in mr:
                        for m in ind:
                            if len(sel)>0:
                                if np.sum(sel==np.where(ind==m)[0][0])>0:
                                    continue
                                else:
                                    if np.sum(mr==np.where(ind==m)[0][0])>0:
                                        continue
                                    else:
                                        if np.average(a_scores[index[ind[l]]:index[ind[l]+1],index[m]:index[m+1]])>NTHRES:
                                            a_dup[k:k+index[m+1]-index[m],:] = a_mod[index[m]:index[m+1],:]
                                            k+=index[m+1]-index[m]
                                            rm.append(np.where(ind==m)[0][0])
                                            sel.append(np.where(ind==m)[0][0])
                            else:
                                if np.sum(mr==np.where(ind==m)[0][0])>0:
                                    continue
                                else:
                                    if np.average(a_scores[index[ind[l]]:index[ind[l]+1],index[m]:index[m+1]])>NTHRES:
                                        a_dup[k:k+index[m+1]-index[m],:] = a_mod[index[m]:index[m+1],:]
                                        k+=index[m+1]-index[m]
                                        rm.append(np.where(ind==m)[0][0])
                                        sel.append(np.where(ind==m)[0][0])
                    addnl_terms = len(rm)-len(mr)
                    rm = mr
                cutoffs.append(k)
                ind = np.delete(ind,rm, axis=0)
            else:
                continue
        
        a_ind  = np.delete(np.delete(a,9,1),9,1)             
        for i in np.arange(len(a)):
            a[i,1] = np.where(np.sum(a_ind[i,4:14] == a_dup, axis=1)==10)[0][0]
        
        cutoffs = np.array(cutoffs)-.1
        for i in np.arange(len(cutoffs)):
            if i==0:
                a[:,1] = np.where((a[:,1]>-0.1) & (a[:,1]<cutoffs[0]), 0, a[:,1])
            else:
                a[:,1] = np.where((a[:,1]>cutoffs[i-1]) & (a[:,1]<cutoffs[i]), i, a[:,1])
        size = np.append(cutoffs[0]+0.1,np.diff(cutoffs))
        for i in np.arange(len(size)):
            Orthologous_Split[a[np.where(a[:,1]==i)[0],3].astype(int),int(OGcount)] = a[np.where(a[:,1]==i)[0],2]
            #Orthologous_Split[a[:,3].astype(int),int(OGcount)] = a[np.where(a[:,1]==i)[0],2]
            b = cumulative_array_str[seq_w_genome[a[np.where(a[:,1]==i)[0],0].astype(int),0].astype(int)]
            ue,count_terms = np.unique(b, return_counts=True)
            sorted_b = ue[np.argsort(count_terms)[::-1]]
            if len(sorted_b)>3:
                Ortholog_metadata[:,int(OGcount)] = sorted_b[0:3]
            else:
                if len(sorted_b)==2:
                    Ortholog_metadata[0:2,int(OGcount)] = sorted_b[0:2]
                else:
                    Ortholog_metadata[0:1,int(OGcount)] = sorted_b[0:1]
            OGcount+=1



Orthologous_Split = Orthologous_Split[:,0:OGcount+1]
Ortholog_metadata = Ortholog_metadata[:,0:OGcount+1]
print("Neighborhood-Based Paralog Separation:--- %s seconds ---" % (time.time() - start_time)) 
print("Paralog Separation Created " +str(OGcount)+" Clusters")
print("[=================>-]")
print(" ")
print(" ")

#This still removes the groups of size 1.
#Also need to create a metric for the neighborhood variation.


count_nan = np.sum(np.isnan(Orthologous_Split),axis=0)
OG_split_ordered = Orthologous_Split[:,np.argsort(count_nan)]
OG_metadata = Ortholog_metadata[:,np.argsort(count_nan)] 
Z = ward(pdist(np.isnan(OG_split_ordered)))
list_out = leaves_list(Z)
dn = hierarchy.dendrogram(hierarchy.linkage(np.isnan(OG_split_ordered), 'single'))
a = os.listdir(Folder)
a2 = [0] * num_genome
plt.close()
OG_split_ordered_DN = OG_split_ordered*0
for i in np.arange(len(dn['leaves'])):
    OG_split_ordered_DN[i,:] = OG_split_ordered[list_out[i],:]
    a2[i] = a[list_out[i]]
plt.imsave(os.path.join(RESULTS_FOLDER,'Pangenome_Visual.png'), OG_split_ordered_DN, cmap='turbo', vmin=0, vmax=30)
figure(figsize=(120, 60), dpi=80)
plt.imshow(OG_split_ordered_DN, cmap='turbo', vmin=0, vmax=30,interpolation='nearest', aspect='auto')
plt.colorbar()
plt.savefig(os.path.join(RESULTS_FOLDER,'Wholistic_Pangenome.png'))
df_met = pd.DataFrame(OG_metadata.T)
df_met = df_met.rename(columns={0: "gene name", 1: "alternative gene name 1",2:"alternative gene name 2"})
df = pd.DataFrame(OG_split_ordered_DN.T, columns = os.listdir(Folder))
df = df.fillna('')
df_concat = pd.concat([df_met, df], axis=1)
df_concat.to_csv(os.path.join(RESULTS_FOLDER,"Pangenome_Numerical.csv"), index=False)
print("Results:--- %s seconds ---" % (time.time() - start_time)) 
print("[==================>]")
print(" ")
print(" ")
print("Done! Thank you for using ARPA.")
sys.exit()

