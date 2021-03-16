#from __future__ import absolute_import

from graph import SGAFilter
import sys
import os
import argparse

from config import *


parser = argparse.ArgumentParser()
parser.add_argument("graph", #type = file, #type = argparse.FileType('r'), 
                    help="input graph asqg in asqg format")
parser.add_argument("control_counts", #type = file, 
                    help="counts from control samples")
parser.add_argument("treated_counts", #type = file, 
                    help="counts from treated samples")
parser.add_argument("merged_counts", #type = file, 
                    help="counts from merging conditions")
parser.add_argument("merged_graph", #type = file, 
                    help="input graph of merged conditions in asqg format")
parser.add_argument("output", type = str, 
                    help="output file to save counts") 

parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")

args = parser.parse_args()


conditions = {k: 1 for k in SAMPLES_TREATED}
conditions.update({k:0 for k in SAMPLES_CONTROL})


if __name__ == "__main__":
    for f in [args.graph, args.control_counts, args.treated_counts, args.merged_graph, args.graph]:
        if not os.path.exists(f):
            raise FileNotFoundError("File %s not found"%f)
    try:
        with open(args.output, 'w') as x:
            pass
    except Exception:
        exit("Cannot create file %s" % args.output)

    
    if args.verbose:
        print("Loading graph:", args.graph)
        sys.stdout.flush()

    sg = SGAFilter.SgaGraph(conditions)
    sg.init_graph(args.graph)
    
    if args.verbose:
        print("Finished loading graph")
        sys.stdout.flush()

    sg.add_duplicates_asqg(args.graph)
    sg.add_duplicates_fasta(args.merged_graph, reverse=True)
    sg.add_duplicates_fasta(args.control_counts)
    sg.add_duplicates_fasta(args.treated_counts)

    #sg.finish_loading_counts()
    if args.verbose:
        print("Finished loading duplicates")

    sg.save_counts(args.output)

    if args.verbose:
        print("Finished saving graph")


