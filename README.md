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

Amino acid sequences from imported genomes (.faa files) are compressed to yield unique alleles. Subsequent encoded genes are iteratively clustered through bulk single linkage clustering by comparing residue counts (vRC; a 20-value vectors in which each term is the total number of each of the 20 amino acid), with cluster boundaries relying upon biological differences in vRCs of different proteins. Resulting homologous proteins are separated by comparing gene neighborhoods, through a process that utilizes results from the prior clustering step, to yield orthologous groups.

**In current first version of ARPA, clustering ("-dsp") has been the primary focus over paralog separation.**

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
$cd ARPA/bin
$chmod +x *.py
$export PATH=$PATH:/path/to/folder/having/ARPA/bin
```
## Usage
### Command-line Options

```
usage: ARPA.py [-h] [-v] [-ct CLUSTER_THRESHOLD] [-dsp] [-cc] [-ld] [-ni NEIGHBORHOOD_IDENTITY]
               [-nt NEIGHBORHOOD_THRESH]
               sequences_folder

Alignment-Free Residue Pangenomic Analysis

positional arguments:
  sequences_folder      directory with .faa files for pangenome

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         print version and exit
  -ct CLUSTER_THRESHOLD, --cluster_threshold CLUSTER_THRESHOLD
                        threshold of similarity used in homolog clustering [Default: 0.95]
  -dsp, --dont_split_paralogs
                        Do not split paralogous groups (no neighborhood analysis)
  -cc, --check_clusters
                        Report the clustered proteins through their identifiers
  -ld, --large_dataset  Dendogram sorting step avoided, full image not printed, and pangenome file saved as
                        Pickle file, not CSV
  -ni NEIGHBORHOOD_IDENTITY, --neighborhood_identity NEIGHBORHOOD_IDENTITY
                        threshold of genes to be reclustered in neighborhood analysis [Default: 0.7]
  -nt NEIGHBORHOOD_THRESH, --neighborhood_thresh NEIGHBORHOOD_THRESH
                        threshold of genes conserved within neighborhood (out of 10) [Default: 3]
  ```
### Input 
ARPA requires an input folder of “.faa” processed amino acid files corresponding to sequenced genomes. Refer to the folder when directing the ARPA code to the set of genomes to construct the pangenome (i.e., "/path/to/folder/with/genomes/"  but NOT "/path/to/folder/with/genomes/*.faa") Encoded proteins are converted into a 1-by-20 row vector (vRC) with each term corresponding to the number of specific amino acid residues present within the sequence. These 20 numbers will be used identify encoded genes during pangenome clustering. ARPA then  compresses the database by removing proteins with identical VRCs. 

### Clustering
ARPA clusters genomes using a bulk single linkage clustering (see preprint for details: [bioRxiv](https://doi.org/10.1101/2022.06.03.494761)). Clustered proteins can be listed if desired using the "-check_clusters" option, especially if using [WhatsGNU](https://doi.org/10.1186/s13059-020-01965-w) [modified files](https://github.com/ahmedmagds/WhatsGNU#command-line-options-for-whatsgnu_database_customizerpy) to compare against [Roary](https://doi.org/10.1093/bioinformatics/btv421).  The default use of a 95% threshold of protein similarity is used unless another value is specificied.

### Paralog Separation [Preliminary]
Paralog separation can also be performed, in a preliminary manner. This involves the separation of paralogous groups through analysis of conserved gene neighborhoods. ARPA uses pre-built clusters from homolog clustering as a starting point for neighborhood analysis. ARPA identifies neighborhoods for all homologs (as opposed to only splitting paralogous groups only when a single genome contains multiple homologs) and splits them into paralogous groups when their neighborhoods differ. This  strategy avoids missed paralogs from differential loss of paralogs within one genome. In comparing neighborhhoods, a similar method of clustering is utilized, with a user-defined value of neighborhood identity (default 70% gene vRCS identity and  >3 of 10 neighborhood genes are similar to be orthologs). 

### Outputs
ARPA typically generates three output files, two of which are visualizations of the pangenome. The most data-rich output is a “.csv” (or “.pickle” for large pangenomics) file, which lists proteins, labelled with observed protein identifiers, and a pangenome output table. The pangenomes table reports the presence or absence of proteins but also lists the “allele deviation scores,” which is calculated as the number of residues different from the protein sequence that is most represented within the homolog group. This output can be utilized within a program like [Scoary](https://doi.org/10.1186/s13059-016-1108-8) or can be modified for phylogenetics (see [bioRxiv preprint](https://doi.org/10.1101/2022.06.03.494761), Fig 3). Non-large database analyses will cluster the visual pangenome into a most similar genomes to reveal patterns within similar genomes through the use of dendogram based solely upon the presence and absence of genes. 

  
## Verification Code
Benchmarking and analytic code is provided within the zipped folder "ARPA_Verification"
Use these codes to recreate analytic and benchmarking results. One of the four codes is a CLI script ("ARPA_Blast_Analysis.py"). For the remaining code, substitute pathways asscociated with your computer into the required subsections and run the code within an IDE (like Spyder). 



## Common Errors
Section in change, with errors that will be worked on to be resolved or avoided in the future. 

The most common observed arror is the addition of non-".faa" files to the genome directory. Doing so will result in error, and this may be especially true for Mac users, wherein a ".DS_Store" hidden file may be added to the folder. Proceed to the directory (use ls -a to see contents of directory) and remove the ".DS_store" file:
  ```
  cd /path/to/folder/having/ARPA_genomes
  rm .DS_Store
  ```
