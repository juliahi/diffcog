import cPickle
import heuristics
from common import *
from common_parameters import *
from compare_heuristics import *
from compare_assemblies import load_sga_scaffold
import os

from graph_stats import get_path_sequences

print "Loading graph", my_renamed_graph_pickle, give_time()
sg = cPickle.load(open(my_renamed_graph_pickle, 'rb'))
print "Finished loading graph", give_time()
try:
    os.mkdir(os.path.dirname(longest_fasta))
except OSError:
    pass


# Heuristic 1: take longest
def init_longest(sg, min_len, min_fc):
    print "Start heuristic 1: longest", give_time()
    longest = heuristics.take_longest(sg)
    longest_seqs = get_path_sequences(sg, longest)
    longest_filter, longest_filter_seqs = filter_by_fc(sg, longest, longest_seqs, min_fc=min_fc)
    longest_filter200, longest_filter200_seqs = filter_by_seqlen(longest_filter, longest_filter_seqs, min_len)
    save_paths_to_fasta(sg, longest_filter200, longest_filter200_seqs, longest_fasta)

    print "paths before FC filter\t %d" % len(longest)
    print "after fc filter"
    stats(sg, longest_filter, longest_filter_seqs)
    print "after lengths filter"
    stats(sg, longest_filter200, longest_filter200_seqs)

    return longest_filter200, longest_filter200_seqs


# Heuristic 2: take longest until foldchange
def init_longest_fc(sg, min_len, min_fc):
    # Heuristic 2: take longest until foldchange
    print "Start heuristic 2: longestfc", give_time()
    longest_fc = heuristics.take_longest_minfc(sg, min_fc)
    longest_fc_seqs = get_path_sequences(sg, longest_fc)
    longest_fc200, longest_fc200_seqs = filter_by_seqlen(longest_fc, longest_fc_seqs, min_len)
    save_paths_to_fasta(sg, longest_fc200, longest_fc200_seqs, longestfc_fasta)

    print "max_lenfc before length filter"
    stats(sg, longest_fc, longest_fc_seqs)
    print "max_lenfc after length filter"
    stats(sg, longest_fc200, longest_fc200_seqs)

    return longest_fc200, longest_fc200_seqs


# Heuristic 3: take best foldchange until > FC
def init_best_fc(sg, min_len, min_fc):
    print "Start heuristic 3: bestfc", give_time()

    best_fc = heuristics.take_best_fc(sg, min_fc)
    best_fc_seqs = get_path_sequences(sg, best_fc)
    best_fc200, best_fc200_seqs = filter_by_seqlen(best_fc, best_fc_seqs, min_len)
    save_paths_to_fasta(sg, best_fc200, best_fc200_seqs, bestfc_fasta)
    
    print "best_fc before length filter"
    stats(sg, best_fc, best_fc_seqs)
    print "best_fc after lengths filter"
    stats(sg, best_fc200, best_fc200_seqs)

    return best_fc200, best_fc200_seqs


def init_sga(sg, filename, min_len, min_fc):
    sga_paths, sga_seqs = load_sga_scaffold(filename, sg)
    sga_paths200, sga_seqs200 = filter_by_seqlen(sga_paths, sga_seqs , min_len)

    sga_paths_filter, sga_seqs_filter = filter_by_fc(sg, sga_paths200, sga_seqs200, min_fc=min_fc)
    save_paths_to_fasta(sg, sga_paths_filter, sga_seqs_filter, sga_fasta)

    print "sga before length filter"
    stats(sg, sga_paths, sga_seqs)
    print "sga after length filter/before fc"
    stats(sg, sga_paths200, sga_seqs200)
    print "sga after fc filter"
    stats(sg, sga_paths_filter, sga_seqs_filter)
    return sga_paths_filter, sga_seqs_filter


init_longest(sg, MIN_LENGTH, MIN_FC)
init_longest_fc(sg, MIN_LENGTH, MIN_FC)
init_best_fc(sg, MIN_LENGTH, MIN_FC)
init_sga(sg, sgafile, MIN_LENGTH, MIN_FC)








