import heapq
from common import *


def take_longest(graph):
    used = {}

    lengths = graph.graph.nodes.data("length")
    paths = []

    heap = [(-v, k) for k, v in lengths]
    heapq.heapify(heap)

    while len(heap) > 0:
        length, nodename = heapq.heappop(heap)
        if nodename in used:
            continue

        path = [nodename]
        used[nodename] = True
        best_fwd = (0, None)
        best_bwd = (0, None)
        first = nodename
        last = nodename

        while True:
            if best_fwd[1] is None:
                for neighb in graph.graph.successors(last):
                    if neighb not in used and lengths[neighb] > best_fwd[0]:
                        best_fwd = (lengths[neighb], neighb)

            if best_bwd[1] is None:
                for neighb in graph.graph.predecessors(first):
                    if neighb not in used and lengths[neighb] > best_bwd[0]:
                        best_bwd = (lengths[neighb], neighb)

            if best_fwd[1] is not None and (best_bwd[1] is None or best_fwd[0] > best_bwd[0]):
                last = best_fwd[1]
                path.append(last)
                best_fwd = (0, None)
                used[last] = True

            elif best_bwd[1] is not None:
                first = best_bwd[1]
                path = [first] + path
                best_bwd = (0, None)
                used[first] = True
            else:
                break

        paths.append(path)
    return paths


def take_longest_minfc(graph, min_fc):
    used = {}

    lengths = graph.graph.nodes.data("length")
    paths = []

    heap = [(-v, k) for k, v in lengths]
    heapq.heapify(heap)
    while len(heap) > 0:
        length, nodename = heapq.heappop(heap)
        if nodename in used:
            continue

        counts = graph.counts[nodename]
        if not foldchange_compare(counts[0], counts[1], min_fc):
            continue

        path = [nodename]
        used[nodename] = True
        first = nodename
        last = nodename

        while True:
            best_bwd = (0, None)
            best_fwd = (0, None)
            for neighb in graph.graph.successors(last):
                if neighb not in used and lengths[neighb] > best_fwd[0]:
                    counts_tmp = graph.counts[neighb]
                    if foldchange_dir(counts,  counts_tmp, min_fc):
                        best_fwd = (lengths[neighb], neighb)

            for neighb in graph.graph.predecessors(first):
                if neighb not in used and lengths[neighb] > best_bwd[0]:
                    counts_tmp = graph.counts[neighb]
                    if foldchange_dir(counts,  counts_tmp, min_fc):
                        best_bwd = (lengths[neighb], neighb)

            if best_fwd[1] is not None and (best_bwd[1] is None or best_fwd[0] > best_bwd[0]):
                last = best_fwd[1]
                path.append(last)
                counts_tmp = graph.counts[last]
                used[last] = True

            elif best_bwd[1] is not None:
                first = best_bwd[1]
                path = [first] + path
                counts_tmp = graph.counts[first]
                used[first] = True
            else:
                break
            counts = [counts[0] + counts_tmp[0], counts[1] + counts_tmp[1]]
        paths.append(path)
    return paths


def take_best_fc(graph, min_fc):

    used = {}
    paths = []

    heap = [(-abslog2foldchange(*v), k) for k, v in graph.counts.items()]
    heapq.heapify(heap)
    while len(heap) > 0:
        lfc, nodename = heapq.heappop(heap)
        if nodename in used:
            continue

        counts = graph.counts[nodename]
        fc = foldchange(*counts)
        if not foldchange_compare(counts[0], counts[1], min_fc):
            continue

        path = [nodename]
        used[nodename] = True
        first = nodename
        last = nodename

        while True:
            best_fwd = (0, None)
            best_bwd = (0, None)
            for neighb in graph.graph.successors(last):
                counts_tmp = graph.counts[neighb]
                if neighb not in used and foldchange_dir(counts, counts_tmp, min_fc):
                    best_fwd = (abslog2foldchange(*counts_tmp), neighb)

            for neighb in graph.graph.predecessors(first):
                counts_tmp = graph.counts[neighb]
                if neighb not in used and foldchange_dir(counts, counts_tmp, min_fc):
                    best_bwd = (abslog2foldchange(*counts_tmp), neighb)

            if best_fwd[1] is not None and (best_bwd[1] is None or
                                            (fc > 1 and best_fwd[0] > best_bwd[0]) or
                                            (fc < 1 and best_fwd[0] < best_bwd[0])):
                last = best_fwd[1]
                path.append(last)
                used[last] = True
                counts_tmp = graph.counts[last]
            elif best_bwd[1] is not None:
                first = best_bwd[1]
                path = [first] + path
                used[first] = True
                counts_tmp = graph.counts[first]
            else:
                break

            counts = [counts[0]+counts_tmp[0], counts[1]+counts_tmp[1]]
            # fc = foldchange(*counts)
        paths.append(path)
    return paths

