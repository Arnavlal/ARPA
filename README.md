[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# ARPA: Alignment-Free Residue Pangenome Analysis
**Arnav Lal, Ahmed Moustafa, and Paul J. Planet**

- University of Pennsylvania, Philadelphia, PA 19104, USA
- Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA 19104, USA
- Children's Hospital of Philadelphia, Philadelphia, PA 19104, USA
- Sackler Institute for Comparative Genomics, American Museum of Natural History, New York, NY 10024, USA

<img width="1128" alt="Screen Shot 2022-05-07 at 2 26 53 678PM" src="https://user-images.githubusercontent.com/66033960/167267136-a3e3d491-bcb0-4f95-8c36-2b286d2f2ddd.png">


## Introduction
ARPA utilizes the amio acid residue counts of proteins for pangenomic protein clustering. <br/>

<img width="1724" alt="Graphical_Summary" src="https://user-images.githubusercontent.com/66033960/178545390-c41fb071-a489-49c6-af43-f1364035b508.png">

**Workflow of Pangenome Construction.**

Amino acid sequences from imported genomes (.faa files) are compressed to yield unique alleles. Subsequent encoded genes are iteratively clustered through bulk single linkage clustering by comparing residue counts (vRC), with cluster boundaries relying upon biological differences in vRCs of different proteins. Resulting homologous proteins are separated by comparing gene neighborhoods, through a process that utilizes results from the prior clustering step, to yield orthologous groups.

In current version, clustering ("-dsp") has been the primary focus over paralog separation.

## Literature

A preprint of ARPA covering motivation, benchmarking and applications is available on bioRxiv at DOI: [10.1101/2022.06.03.494761](https://doi.org/10.1101/2022.06.03.494761)


## Installation and Dependencies
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

### Command-line Options

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
  
  
### Common Errors
The most common observed arror is the addition of non-".faa" files to the genome directory. Doing so will result in error, and this may be especially true for Mac users, wherein a ".DS_Store" hidden file may be added. Proceed to the directory (use ls -a to see contents of directory) and remove the ".DA_store" file:
  ```
  cd /path/to/folder/having/ARPA_genomes
  rm .DS_Store
  ```
