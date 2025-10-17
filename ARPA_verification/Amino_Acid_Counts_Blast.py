#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: arnavlal
"""



import numpy as np
import matplotlib.pyplot as plt
import os
import collections
from Bio import SeqIO
from tqdm import tqdm
import time
start_time = time.time()


loc_ninety_input = np.loadtxt('/Users/arnavlal/e_coli_blast.txt', dtype=int)
e_coli_path = '/Path to E coli files/GCA_000005845.2/protein.faa'

seqs = collections.defaultdict(list) 


name = e_coli_path
pathToFile = open(name)
index = 0
for seq_record in SeqIO.parse(pathToFile, """fasta"""):
    if index == 0:
        name_label = seq_record.description.split("[")[1].split("]")[0]
    index +=1
    seqs[name_label].append(seq_record.seq)

import subprocess
from Bio.SeqRecord import SeqRecord
results = []
for i in np.arange(len(loc_ninety_input)):
    my_records = SeqRecord(seqs['Escherichia coli str. K-12 substr. MG1655'][loc_ninety_input[i,0]], id="query1")
    SeqIO.write(my_records, "/Users/arnavlal/query.faa", "fasta")
    my_records = SeqRecord(seqs['Escherichia coli str. K-12 substr. MG1655'][loc_ninety_input[i,1]], id="query2")
    SeqIO.write(my_records, "/Users/arnavlal/query2.faa", "fasta")
    result = subprocess.run(["blastp -query /Users/arnavlal/query.faa -subject /Users/arnavlal/query2.faa -outfmt '6 qseqid qlen sseqid slen length pident qstart qend sstart send evalue bitscore qcovhsp qcovs'"], shell=True, capture_output=True)
    results.append(result.stdout)
    
with open('blast_final_e_ecoli.txt', 'w') as f:
    for line in results:
        f.write(f"{line}\n")
