import SGAFilter
import graph_stats
import sys
import argparse
import os
from datetime import datetime

from common import *


parser = argparse.ArgumentParser()
parser.add_argument("graph", #type = file, #type = argparse.FileType('r'), 
                    help="input graph asqg in asqg format")
parser.add_argument("counts", #type = file, 
                    help="counts for each vertex")
parser.add_argument("output_graph", #type = file, 
                    help="output file to save simplified graph") 
parser.add_argument("output_counts", #type = file, 
                    help="output file to save counts") 

def check_positive(value):
    try:
        ivalue = int(value)
    except:
        argparse.ArgumentTypeError("%s must be positive integer" % value)
    if ivalue <= 0:
        raise argparse.ArgumentTypeError("%s must be positive integer" % value)
    return ivalue

parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
parser.add_argument("-r", "--DEADENDS_REMOVE_ROUNDS", type=int, default=10, choices=range(1, 11),
                    help="number of iterations of deadends removal")
parser.add_argument("-d", "--DEADENDS_MIN_LENGTH", type=check_positive, default="100",
                    help="length of dead-end to remove")
parser.add_argument("-m", "--MIN_LENGTH", type=check_positive, default="200",
                    help="minimal contig length")

args = parser.parse_args()



def print_log(info):
    if args.verbose:
        print(str(datetime.now()) + ": " + info)
        sys.stdout.flush()

if __name__ == "__main__":

    for f in [args.graph, args.counts]:
        if not os.path.exists(f):
            raise FileNotFoundError("file %s not found" % f)

    print_log("Loading graph:%s \n with parameters: DEADENDS_REMOVE_ROUNDS = %d, DEADENDS_MIN_LENGTH = %d, \
         MINIMAL CONTIG LENGHT = %d" % (args.graph, args.DEADENDS_REMOVE_ROUNDS, args.DEADENDS_MIN_LENGTH, args.MIN_LENGTH))

    try: 
        with open(args.output_graph, 'w') as x: pass
        with open(args.output_counts, 'w') as x: pass
    except e: raise(e)

    sg = SGAFilter.SgaGraph()
    sg.init_graph(args.graph, verbose=args.verbose)
    sg.load_graph(args.counts, verbose=args.verbose)
    
    print_log("Finished loading graph")

    if args.verbose:
        graph_stats.short_summary(sg)
        sys.stdout.flush()

    sg.remove_short_islands(args.MIN_LENGTH, verbose=args.verbose)

    ### simplify graph and compute graph statistics
    sg.compress_simple_paths(verbose=args.verbose)
    sg.remove_short_islands(args.MIN_LENGTH, verbose=args.verbose)


    print_log("Finished simplifying paths")
    if args.verbose:
        graph_stats.short_summary(sg)
        sys.stdout.flush()
    
    # REMOVE DEAD-ENDS and simplify
    for i in range(1, args.DEADENDS_REMOVE_ROUNDS+1):
        sg.remove_deadends_by_length(args.DEADENDS_MIN_LENGTH, verbose=args.verbose)

        print_log("Finished dead-ends round %d" % i)

        if args.verbose:
            graph_stats.short_summary(sg)
            sys.stdout.flush()

    sg.compress_simple_paths(verbose=args.verbose)
    sg.remove_short_islands(args.MIN_LENGTH, verbose=args.verbose)

    print_log("Finished all simplification")
    if args.verbose:
        graph_stats.short_summary(sg)
        sys.stdout.flush()
    
    ### save graph in ASQG
    sg.write_to_asqg(args.output_graph, contigs=None, rename=True)
    sg.save_counts(args.output_counts)

    print_log("Finished saving graph")

