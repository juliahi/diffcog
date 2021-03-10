# Heuristic algorithm to build differential contigs from Overlap/String graph build for two sets of samples from RNA-seq experiments.

This repository contains code to build differential contigs with String graph
using two sets of samples containing NGS reads.
They need to be specified in *config.py* file
by specifying input and output directory, two sets of sample names and input files' suffixes for paired-end reads. It also 

## Build graph using SGA -- run_sga

Use *SGA* to build String graph in specific way to be able to retrieve
read counts from each condition. It assumes reads are the same length (i.e. untrimmed).

This also prepares counts for each condition for each node in graph.

## Simplify graph -- graph

Simplify graph for RNA-seq (assuming strand-specific) and save simplified graph for further analysis. This needs graph in ASQG format and table with counts in tsv, which is prepared in run_sga/counts_from_sga.py, but can be provided from other source.

## Run heuristics

Available heuristic methods to construct contigs from simplified string graph with specified fold-change threshold FC:
* longest - greedily extend path by adding longest vertex
* longestfc - greedily add longest vertex while foldchange is larger than FC (or smaller than 1/FC)
* bestfc - greedily extend path by vertex keeping the furthest fold-change from FC


