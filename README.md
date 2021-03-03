# Heuristic algorithm to build differential contigs from Overlap/String graph build for two sets of samples from RNA-seq experiments.

This repository contains code to build differential contigs with String graph
using two sets of samples containing NGS reads.
They need to be specified in config.py file
by specifying input and output directory, two sets of sample names and input files' suffixes for paired-end reads.

## Build graph using SGA -- run_sga

Use *SGA* to build String graph in specific way to be able to retrieve
read counts from each condition. It assumes reads are the same length (i.e. untrimmed).

## Simplify graph -- graph

Simplify graph for RNA-seq (assuming strand-specific) and save simplified graph for further analysis.

## Run heuristics

Choose heuristic method to construct contigs from simplified string graph with specified fold-change threshold FC:
* longest - greedily extend path by adding longest vertex
* longestfc - greedily add longest vertex while foldchange is larger than FC (or smaller than 1/FC)
* bestfc - greedily extend path by vertex keeping the furthest fold-change from FC
