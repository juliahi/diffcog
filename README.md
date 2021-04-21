# DIFFCOG
-  Heuristic algorithm to build *DIFFerential Contigs from Overlap (string) Graph* build for two sets of samples from RNA-seq experiments.

This repository contains code to build differential contigs with String graph
using two sets of samples containing NGS paired-end RNA-seq reads.


## Requirements
* *conda* -- Miniconda, a minimal version of conda can be installed from https://docs.conda.io/en/latest/miniconda.html

## Installation:
* Clone this repository:

```
git clone https://github.com/juliahi/diffcog.git
cd diffcog  
```

* Create conda environment -- this will install all required packages:  snakemake, Python packages: networkx, numpy, and other software: seqkit and sga (String Graph Assembler: https://github.com/jts/sga/blob/master/README.md)

```
conda env create -f diffcog.yaml
```

* Activate environment:

```
conda activate diffcog
```

* Run on provided test data specified in *config.py* to check if installation was succesfull:

```
./run_all.sh
```

### Alternative installation 
If you don't want to use conda for managing dependencies, you can put path to SGA into *config.py* file and provide snakemake, seqkit and Python with installed dependencies available from path. 

## Run DIFFCOG
Specify sets of reads in **config.py** file by specifying:
* two sets of sample names and input files' suffixes for paired-end reads
* input and output directory
* optional: other parameters such as number of simplification iterations or minimal contig length.

Run whole analysis with **./run_all.sh**: build graph, simplify it and run heuristics with Snakemake. For parts of workflow see next sections.

### Build graph using SGA -- directory run_sga

Use *SGA* to build String graph in specific way to be able to retrieve
read counts from each condition. It assumes reads are the same length (i.e. untrimmed).

This also prepares counts for each condition for each node in graph.

``` cd run_sga
snakemake all
```

Although it is possible to run heuristics on graph in ASQG format provided from other source we do not recommend that.

### Simplify graph -- directory graph

Simplify graph for RNA-seq (assuming strand-specific) and save simplified graph for further analysis. This needs graph in ASQG format and table with counts in tsv, which is prepared in run_sga/counts_from_sga.py, but can be provided from other source.

To run simplification step run:
``` cd overlap_graph
snakemake simplify
```

### Run heuristics -- directory graph

Available heuristic methods to construct contigs from simplified string graph with specified fold-change threshold FC:
* longest - greedily extend path by adding longest vertex
* longestfc - greedily add longest vertex while foldchange is larger than FC (or smaller than 1/FC)
* bestfc - greedily extend path by vertex keeping the furthest fold-change from FC

To run heuristics specified in *config.py* file run *snakemake heuristics* in graph directory:

``` cd overlap_graph
snakemake heuristics
```

Reads are normalized for each condition by sum(condition)/10^6.
