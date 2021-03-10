import SGAFilter
import cPickle
import time
import graph_stats
import sys


from common import *

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("graph", type = file, #type = argparse.FileType('r'), 
                    help="input graph asqg in asqg format")
parser.add_argument("-c", "--counts", type = file, 
                    help="counts for each vertex")
parser.add_argument("-og", "--output_graph", type = file, 
                    help="output file to save simplified graph") 
parser.add_argument("-oc", "--output_counts", type = file, 
                    help="output file to save counts") 

parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")

args = parser.parse_args()

if __name__ == "__main__":

    for f in [args.graph, args.c]:
        if not os.path.exists(f):
            raise FileNotFoundError("file %s not found" % f)

    if args.v:
        print("Loading graph:", filename)
        sys.stdout.flush()

    try: with open(f, 'w') as x: pass
    except e: raise(e)

    sg = SGAFilter.SgaGraph(conditions)
    sg.init_graph(args.graph)
    sg.load_graph()
    sg.load_counts(args.c)
    
    if args.v:
        print("Finished loading graph")
        graph_stats.short_summary(sg)
        sys.stdout.flush()


    ### simplify graph and compute graph statistics
    sg.compress_simple_paths()
    sg.remove_short_islands(MIN_LENGTH)
    if args.v:
        print("Finished simplifying paths")
        graph_stats.short_summary(sg)
        sys.stdout.flush()
    
    # REMOVE DEAD-ENDS and simplify
    for i in xrange(DEADENDS_REMOVE_ROUNDS):
        sg.remove_deadends_by_length(DEADENDS_MIN_LENGTH)
        if args.v:
            print("Finished dead-ends round %d" % i)
            graph_stats.short_summary(sg)
            sys.stdout.flush()

    sg.compress_simple_paths()
    sg.remove_short_islands(MIN_LENGTH)

    if args.v:
        print "Finished all simplification", give_time()
        graph_stats.short_summary(sg)
        sys.stdout.flush()
    
    ### save graph in ASQG
    sg.write_to_asqg(args.og, contigs=None, rename=True)
    sg.save_counts(args.oc)

    if args.v:
        print("Finished saving graph")

