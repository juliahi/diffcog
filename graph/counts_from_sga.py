import SGAFilter
import cPickle
import time
import graph_stats
import sys
import os

from common import *

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("graph", type = file, #type = argparse.FileType('r'), 
                    help="input graph asqg in asqg format")
parser.add_argument("-c", "--control_counts", type = file, 
                    help="counts from control samples")
parser.add_argument("-t", "--treated_counts", type = file, 
                    help="counts from treated samples")
parser.add_argument("-m", "--merged_counts", type = file, 
                    help="counts from merging conditions")
parser.add_argument("-g", "--merged_graph", type = file, 
                    help="input graph of merged conditions in asqg format")
parser.add_argument("-o", "--output", type = file, 
                    help="output file to save counts") 

parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")

args = parser.parse_args()


if __name__ == "__main__":

    for f in [args.graph, args.c, args.t, args.m, args.g]:
        if not os.path.exists(f):
            raise FileNotFoundError("File %s not found"%f)
    try:
        with open(f, 'w') as x:
            pass
    except e:
        raise(e)

        
    if args.v:
        print "Loading graph:", filename, give_time()
        sys.stdout.flush()

    sg = SGAFilter.SgaGraph(conditions)
    sg.init_graph(args.graph)
    
    if args.v:
        print "Finished loading graph", give_time()
        sys.stdout.flush()

    sg.add_duplicates_asqg(args.g)
    sg.add_duplicates_fasta(args.m, reverse=True)
    sg.add_duplicates_fasta(args.c)
    sg.add_duplicates_fasta(args.t)

    #sg.finish_loading_counts()
    if args.v:
        print("Finished loading duplicates")

    sg.save_counts(args.o)

    if args.v:
        print("Finished saving graph")


