# widely used methods methods
import math
import numpy
import time

epsilon = 0.000001

PATH_SEP = '&'
NODE_SEP = '|'
JOINED_NODE_SEP = '~'


def give_time():
    return time.asctime(time.localtime(time.time()))


def compl(s):
    if s == '' or s == []: return ''
    if s is None: return None

    def c(x):
        if x == 'A':
            return 'T'
        elif x == 'C':
            return 'G'
        elif x == 'G':
            return 'C'
        elif x == 'T':
            return 'A'
        return 'N'

    return ''.join(map(c, s))


def cons_pairs(input_list):
    for i in xrange(len(input_list)-1):
        yield (input_list[i], input_list[i+1])


def index_pairs(n):
    for j in xrange(1, n):
        for i in xrange(j):
            yield (i, j)

# def foldchange(n1, n2):
#         if n2 == 0: return float("inf")
#         return 1. * n1 / n2


def foldchange(n1, n2):
    if n2 == 0: n2 = epsilon
    if n1 == 0: n1 = epsilon
    return 1. * n1 / n2


# def abslog2foldchange(n1, n2):
#         return abs(math.log(foldchange(n1, n2), 2)) if n1 != 0 else float("-inf")
def abslog2foldchange(n1, n2):
        return abs(math.log(foldchange(n1, n2), 2))


def foldchange_compare(n1, n2, min_fc):
    # check if 1/min_fc <= n1/n2 <= min_fc
    if n1 == n2 == 0: return False
    if n1 >= n2*min_fc: return True
    if n2 >= n1*min_fc: return True
    return False


# check if foldchange for counts c1 stays enriched in the same direction after adding counts c2
def foldchange_dir(c1, c2, min_fc):
        if c1[0] >= c1[1]*min_fc and c1[0]+c2[0] >= (c1[1]+c2[1])*min_fc:
            return True
        if c1[1] >= c1[0]*min_fc and c1[1]+c2[1] >= (c1[0]+c2[0])*min_fc:
            return True
        return False


def iter_fasta(filename):
    name = None
    seq = ""
    for line in open(filename):
        if line.startswith(">"):
            if name is not None:
                yield name, seq
            name = line.strip()[1:]
            seq = ""
        else:
            seq += line.strip()

    if name is not None:
        yield name, seq


def write_to_fasta(names, sequences, filename):
    with open(filename, 'w+') as f:
        for name, seq in zip(names, sequences):
            f.write(">%s\n%s\n" % (name, seq))


def variance(l):
    #m = mean(l)
    #return sum([(x - m)**2 for x in l])/len(l)
    return numpy.var(l)


def dispersion(l):
    if numpy.mean(l) == 0:
        return 0     # TODO: if all values are 0 there is no dispersion... is it right?
    return variance(l) / numpy.mean(l)


jaccard_header = "Set1 \ set2\tSet1 n Set2\tSet2 \ set1\tJaccard index"


def jaccard_index(set1, set2):
    wspolne = len(set1.intersection(set2))
    idx = 1. * wspolne / len(set1.union(set2))
    #print "Set1\set2=", len(set1) - wspolne, '\t', "Set1 n Set2=", wspolne, '\t', "Set2\set1=", len(set2) - wspolne
    #print "%n\t%d\t%d\t%f" % (len(set1) - wspolne, wspolne, len(set2) - wspolne, idx)
    print '{:0,d}\t{:0,d}\t{:0,d}\t{:0,.4f}'.format(len(set1) - wspolne, wspolne, len(set2) - wspolne, idx)
    return idx