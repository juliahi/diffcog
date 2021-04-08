
import heuristics
import SGAFilter
import argparse

from graph_stats import get_path_sequences, get_path_count, get_path_coverage
from common import *


parser = argparse.ArgumentParser()
parser.add_argument("graph", #type = file, #type = argparse.FileType('r'), 
                    help="input graph asqg in asqg format")
parser.add_argument("counts", #type = file, 
                    help="counts for each vertex")
parser.add_argument("-o", "--output_dir", help="directory to write output FASTA files") 

parser.add_argument("-r", "--heuristics", nargs='+', required=True, choices=("longest", "longestfc", "bestfc"), #type = file, 
                    help="types of heuristics") 
parser.add_argument("-fc", "--fold_change", type=float,  default=2,
                    help="minimal fold change of differential contigs") 
parser.add_argument("-l", "--minimal_length", type=int, default=200,
                    help="minimal contig length") 
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")

args = parser.parse_args()

def stats(graph, paths, seqs):
    counts = list(map(lambda x: get_path_count(graph, x), paths))

    unnormed = (sum(map(lambda x: x[0], counts)) * graph.sums[0] + \
           sum(map(lambda x: x[1], counts)) * graph.sums[1]) / 1000000
    print("No. paths: \tNo. vertices: \tLength (bp): \tCounts0: \tCounts1: \tCounts:")
    print(f"{len(paths)}\t\t{sum([len(path) for path in paths])} " + 
        f"\t{sum(map(len, seqs))}\t{sum(map(lambda x: x[0], counts))} " +
        f"\t{sum(map(lambda x: x[1], counts))}\t{sum(map(lambda x: sum(x), counts))}")
    print("\t counts from graph [fraction]")
    print(f"{sum(map(lambda x: sum(x), counts)) / sum(graph.node_nreads_sums())}")


def save_paths_to_fasta(sg, paths, seqs, out_file, tsv_file):        # writing to FASTA
        counts = 0.
        length = 0
        total = 0
        with open(out_file, 'w+') as f:
            with open(tsv_file, 'w+') as tsv:
                for path, seq in zip(paths, seqs):
                    counts_tmp = get_path_count(sg, path)
                    counts += counts_tmp[0] + counts_tmp[1]
                    length += len(seq)
                    total += 1
                    f.write('>%s counts1=%f counts2=%f foldchange=%f\n%s\n' % (
                            PATH_SEP.join(path), counts_tmp[0], counts_tmp[1],
                            foldchange(counts_tmp[0], counts_tmp[1]), seq))
                    tsv.write(f"{PATH_SEP.join(path)}\t{counts_tmp[0]}\t{counts_tmp[1]}\n")
        if args.verbose:
            print("Saved %d sequences with %d counts of total length %d" % (total, counts, length))


def load_info_from_fasta(filename):
    paths = []
    seqs = []
    counts = []

    for name, seq in iter_fasta(filename):
        name, c1, c2, fc = name.split()
        paths.append(name.split(PATH_SEP))
        counts.append((float(c1.split('=')[1]), float(c2.split('=')[1])))
        seqs.append(seq)
    return paths, seqs, counts

def load_paths_from_fasta(filename):
    paths = []
    seqs = []

    for name, seq in iter_fasta(filename):
        name = name.strip().split()
        paths.append(name.split(PATH_SEP))
        seqs.append(seq)
    return paths, seqs

### Filtering selected paths ###

def filter_by_fc(graph, paths, sequences, min_fc, new_names=False):
        if len(paths) == 0:
            return []

        counts = 0
        filtered_paths = []
        filtered_seqs = []
        for path, s in zip(paths, sequences):
            counts_tmp = get_path_count(graph, path, new_names)
            if not foldchange_compare(counts_tmp[0], counts_tmp[1], min_fc):
                continue
            counts += counts_tmp[0] + counts_tmp[1]
            filtered_paths.append(path)
            filtered_seqs.append(s)

        print("Leave %d sequences with %d counts" % (len(filtered_paths), counts))
        return filtered_paths, filtered_seqs

def filter_by_seqlen(paths, seqs, min_len):
    filtered_paths = []
    filtered_seqs = []
    for p, s in zip(paths, seqs):
        if len(s) >= min_len:
            filtered_paths.append(p)
            filtered_seqs.append(s)
    return filtered_paths, filtered_seqs


if __name__ == "__main__":

    for h in args.heuristics:
        try:
            with open(f"{args.output_dir}/{h}.fa", 'w') as x:
                pass
        except Exception:
            exit(f"Cannot create file {args.output_dir}/{h}.fa")

    # load graph
    sg = SGAFilter.SgaGraph()
    sg.init_graph(args.graph, verbose=args.verbose)
    sg.load_graph(args.counts, verbose=args.verbose)
    sg.normalize_counts(verbose=args.verbose)
    
    for h in args.heuristics:
        if h == "longest":
            contigs = heuristics.take_longest(sg)
        elif h == "longestfc":
            contigs = heuristics.take_longest_minfc(sg, args.fold_change)
        else: # h = bestfc
            contigs =  heuristics.take_best_fc(sg, args.fold_change)

        seqs = get_path_sequences(sg, contigs)

        if args.verbose:
            print(f"{h} before filtering by length")
            stats(sg, contigs, seqs)
        contigs, seqs = filter_by_seqlen(contigs, seqs, args.minimal_length)
        save_paths_to_fasta(sg, contigs, seqs, f"{args.output_dir}/{h}.fa", 
                                                f"{args.output_dir}/{h}_counts.tsv")

        if args.verbose:
            print(f"{h} after filtering by length")
            stats(sg, contigs, seqs)


