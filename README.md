[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# ARPA: Alignment-Free Residue Pangenome Analysis
By: Arnav Lal, Ahmed Moustafa, and Paul J. Planet

<img width="1128" alt="Screen Shot 2022-05-07 at 2 26 53 678PM" src="https://user-images.githubusercontent.com/66033960/167267136-a3e3d491-bcb0-4f95-8c36-2b286d2f2ddd.png">


## Introduction
ARPA utilizes the amio acid residue counts of proteins for comparative gene analysis. Pangenomic creation is <br/>
WhatsGNU compresses proteins database based on exact match to much fewer number of proteins that differ by at least one amino acid. WhatsGNU will save a copy of the compressed database in two formats; database.txt and database.pickle for faster subsequent uses.<br/>

### Installation and Dependencies
* [Python3.x](https://www.python.org/downloads/)
* NumPy
* SciPy (Scipy.cluster)
* Matplotlib
* Pandas
* Bio (SeqIO)
* General Dependencies: time, os, sys, argparse

```
$git clone https://github.com/Arnavlal/ARPA
$cd ARPA
$chmod +x *.py
$export PATH=$PATH:/path/to/folder/having/ARPA
```
