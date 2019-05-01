import csv
from RMI import *
from Bio import SeqIO
import struct
import numpy as np
import time

EXPERTS_LEVELS = [
    [100],
    [1000],
    [10, 10],
    [10, 100],
    [10, 1000],
    [100, 100],
    [100, 1000],
    [1000, 1000],
    [1000, 10000],
]

def avg(lst):
    return sum(lst) / len(lst)

class Args:
    pass

def run_train(csvwriter, modelcsvwriter, experts_levels):
    print("Doing  experts level", experts_levels)
    start_time = time.time()

    INT32_SIZE = 4
    NUM_RUNS = 6
    NUM_DROPS = 1

    args = Args()
    args.bwt_file = "../data/chr1.bwt"
    args.ref_seq = "../data/pp_chr1.fa"
    args.query_db = "../data/chr1_small_1k.fa"
    args.query_size = 21
    args.output_file = "output"

    rmi = RMI(experts_levels)

    ref_seq = ""
    for record in SeqIO.parse(args.ref_seq, "fasta"):
        ref_seq = record.seq

    sa_bwt = []

    with open(args.bwt_file, 'rb') as bwt_file:
        for i in range(len(ref_seq) + 1):
            sa_bwt.append(struct.unpack("<L", bwt_file.read(INT32_SIZE))[0])

    queries = []
    for record in SeqIO.parse(args.query_db, "fasta"):
        queries.append(record.seq)

    read_time = time.time()
    print("Done reading")

    query_size = args.query_size

    x = []
    y = []
    nucleo = {}
    nucleo['A'] = 0
    nucleo['C'] = 1
    nucleo['G'] = 2
    nucleo['T'] = 3
    for i in range(len(ref_seq) + 1):
        if sa_bwt[i] + query_size > len(ref_seq):
            continue
        ref_subseq = ref_seq[(sa_bwt[i]):(sa_bwt[i] + query_size)]
        subseq = 0
        for j in range(query_size):
            subseq = subseq << 2 | nucleo[ref_subseq[j]]

        x.append(subseq)
        y.append(i)

    rmi.fit(np.asarray(x).reshape(-1, 1), np.asarray(y))

    fit_time = time.time()
    print("done fitting")

    pred_times = []
    for i in range(NUM_RUNS):
        curr_time = time.time()
        x = []
        for query in queries:
            query_int = 0
            for j in range(query_size):
                query_int = query_int << 2 | nucleo[query[j]]
            x.append(query_int)

        y = rmi.predict(np.asarray(x).reshape(-1, 1))
        predict_time = time.time()

        if i >= NUM_DROPS:
            pred_times.append(predict_time - curr_time)

    print("done predicting")

    def get_ref_int(ind):
    	ref_subseq = get_ref_seq(ind)
    	ref_int = 0
    	for i in range(query_size):
        	ref_int = ref_int << 2 | nucleo[ref_subseq[i]]
    	return ref_int

    def get_ref_seq(ind):
        if sa_bwt[ind] + query_size > len(ref_seq):
            return None
        return ref_seq[sa_bwt[ind]:(sa_bwt[ind] + query_size)]

    def binary_search(query_str, lower, upper, strict):  # strict false means look for lower bound
        in_low = lower
        in_hi = upper
        while get_ref_seq(lower) is None:
            lower += 1
        while get_ref_seq(upper) is None:
            upper -= 1
        if get_ref_seq(lower) > query_str:
            return in_low - 1
        if get_ref_seq(upper) < query_str:
            return in_hi + 1
        if lower == upper:
            return lower
        mid = (lower + upper) // 2
        while get_ref_seq(mid) is None and mid >= lower:
            mid -= 1
        mid_seq = get_ref_seq(mid)
        if mid_seq < query_str or (mid_seq == query_str and strict):
            return binary_search(query_str, mid + 1, upper, strict)
        else:
            return binary_search(query_str, lower, mid, strict)

    def exponential_search(query_str, start_sa):
        lower_bound = None
        upper_bound = None
        curr_seq = get_ref_seq(start_sa)
        while curr_seq is None:
            start_sa += 1
            curr_seq = get_ref_seq(start_sa)
        if curr_seq < query_str:
            lower_bound = start_sa
        if curr_seq > query_str:
            upper_bound = start_sa
        window_size = 1
        ind = start_sa
        if upper_bound is None:
            while start_sa + window_size < len(ref_seq) + 1:
                ind = start_sa + window_size
                window_size = 2 * window_size
                found_seq = get_ref_seq(ind)
                while found_seq is None:
                    ind += 1
                    found_seq = get_ref_seq(ind)
                if found_seq > query_str:
                    upper_bound = ind
                    break
                if found_seq < query_str:
                    lower_bound = ind
        window_size = 1
        ind = start_sa
        if lower_bound is None:
            while start_sa - window_size >= 0:
                ind = start_sa - window_size
                window_size = 2 * window_size
                found_seq = get_ref_seq(ind)
                while found_seq is None:
                    ind -= 1
                    found_seq = get_ref_seq(ind)
                if found_seq < query_str:
                    lower_bound = ind
                    break
                if found_seq > query_str:
                    upper_bound = ind
        if lower_bound is None:
            lower_bound = 0
        if upper_bound is None:
            upper_bound = len(sa_bwt) - 1
        return (binary_search(query_str, lower_bound, upper_bound, False),
                binary_search(query_str, lower_bound, upper_bound, True) + 1)

    search_fn = exponential_search

    pred_dists = []
    with open(args.output_file, "wb") as output_file:
        for i in range(0, len(queries)):
            start_sa = int(y[i])
            if start_sa < 0:
                start_sa = 0
            if start_sa > len(ref_seq):
                start_sa = len(ref_seq)
            (l, u) = search_fn(queries[i], start_sa)
            pred_dists.append(max(l - start_sa, start_sa - u, 0))
            output_file.write(struct.pack("<L", l))
            output_file.write(struct.pack("<L", min(u, len(ref_seq))))

    search_times = []
    for i in range(NUM_RUNS - NUM_DROPS):
        curr_time = time.time()
        for i in range(0, len(queries)):
            start_sa = int(y[i])
            if start_sa < 0:
                start_sa = 0
            if start_sa > len(ref_seq):
                start_sa = len(ref_seq)
            (l, u) = search_fn(queries[i], start_sa)
        search_time = time.time()
        search_times.append(search_time - curr_time)

    print("done searching")

    csvwriter.writerow([experts_levels, sum(pred_dists) / len(pred_dists), max(pred_dists), read_time - start_time, fit_time - read_time, avg(pred_times), avg(search_times)])

    for dist in pred_dists:
        modelcsvwriter.writerow([dist])

    print(experts_levels)
    print("Average prediction distance: ", sum(pred_dists) / len(pred_dists))
    print("Maximum prediction distance: ", max(pred_dists))
    print("I/O Read Time (s): ", read_time - start_time)
    print("Fit Time (s): ", fit_time - read_time)
    print("Predict Time (s): ", avg(pred_times))
    print("Search Time (s): ", avg(search_times))

with open('results.csv', 'w', newline='\n') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')
    csvwriter.writerow(['experts_levels', 'io_read_time', 'fit_time', 'pred_time', 'search_time'])
    for experts_config in EXPERTS_LEVELS:
        with open('results_' + str(experts_config) + '.csv', 'w', newline='\n') as modelcsvfile:
            modelcsvwriter = csv.writer(modelcsvfile)
            modelcsvwriter.writerow(['distance'])
            run_train(csvwriter, modelcsvwriter, experts_config)


