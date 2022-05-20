[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# ARPA: Alignment-Free Residue Pangenome Analysis
By: Arnav Lal, Ahmed Moustafa, and Paul J. Planet

<img width="1128" alt="Screen Shot 2022-05-07 at 2 26 53 678PM" src="https://user-images.githubusercontent.com/66033960/167267136-a3e3d491-bcb0-4f95-8c36-2b286d2f2ddd.png">


## Introduction
ARPA utilizes the amio acid residue counts of proteins for pangenomic protein clustering. <br/>

In current version, clustering ("-dsp") has been the primary focus.


### Installation and Dependencies
* [Python3.x](https://www.python.org/downloads/)
* NumPy
* SciPy (Scipy.cluster)
* Matplotlib
* Pandas
* Bio (SeqIO)
* Pickle (if creating pangenomes for large datasets)
* General Dependencies: time, os, sys, argparse

```
$git clone https://github.com/Arnavlal/ARPA
$cd ARPA
$chmod +x *.py
$export PATH=$PATH:/path/to/folder/having/ARPA
```

###Command-line Options

```
usage: ARPA.py [-h] [-v] [-ct CLUSTER_THRESHOLD] [-m METRIC] [-dsp] [-ld] [-ni NEIGHBORHOOD_IDENTITY] [-nt NEIGHBORHOOD_THRESH]
               sequences_folder

Alignment-Free Residue Pangenomic Analysis

positional arguments:
  sequences_folder      directory with .faa files for pangenome

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         print version and exit
  -ct CLUSTER_THRESHOLD, --cluster_threshold CLUSTER_THRESHOLD
                        threshold of similarity used in homolog clustering [Default: 0.95]
  -m METRIC, --metric METRIC
                        metric for genes: 'mr' diff. From most represented, 'aa' diff. from all alleles, 'aan' diff. from all
                        alleles normalized by counts [Default: 'mr']
  -dsp, --dont_split_paralogs
                        Do not split paralogous groups (no neighborhood analysis)
  -ld, --large_dataset  Dendogram sorting step avoided, full image not printed, and pangenome file saved as Pickle file, not CSV
  -ni NEIGHBORHOOD_IDENTITY, --neighborhood_identity NEIGHBORHOOD_IDENTITY
                        threshold of genes to be reclustered in neighborhood analysis [Default: 0.7]
  -nt NEIGHBORHOOD_THRESH, --neighborhood_thresh NEIGHBORHOOD_THRESH
                        threshold of genes conserved within neighborhood (out of 10) [Default: 3]
  ```
