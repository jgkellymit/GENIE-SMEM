from RMI import *
from Bio import SeqIO
import numpy as np
import time
import json

start_time = time.time()

fm_file = "../data/medium_data-FM.json"
ref_seq_file = "../data/medium_data.fa"
query_db = "../data/query10.fa"
query_size = 10
output_file = "output.txt"

rmi = RMI([10, 100])

ref_seq = ""
for record in SeqIO.parse(ref_seq_file, "fasta"):
    ref_seq = record.seq


with open(fm_file, 'r') as fm_f:
    fm_index = json.load(fm_f)
    suffix_array = fm_index["suffix_array"]



queries = []
for record in SeqIO.parse(query_db, "fasta"):
    queries.append(record.seq)

read_time = time.time()


x = []
y = []
nucleo = {}
nucleo['A'] = 0
nucleo['C'] = 1
nucleo['G'] = 2
nucleo['T'] = 3

print(suffix_array)
# TCTGACCTGA GGAGAACTGT GCTCCGCCTT CAGAGTACCA CCGAAATCTG

for i in range(len(ref_seq) + 1):
    if suffix_array[i] + query_size > len(ref_seq):
        continue
    ref_subseq = ref_seq[(suffix_array[i]) - 1:(suffix_array[i] + query_size) - 1]
    subseq = 0
    for j in range(query_size):
        subseq = subseq << 2 | nucleo[ref_subseq[j]]

    x.append(subseq)
    y.append(i)

rmi.fit(np.asarray(x).reshape(-1,1), np.asarray(y))

fit_time = time.time()

x_predict = []


for query in queries:
    print(query)
    query_int = 0
    for j in range(query_size):
        query_int = query_int << 2 | nucleo[query[j]]
    x_predict.append(query_int)


y_predict = rmi.predict(np.asarray(x_predict).reshape(-1, 1))


predict_time = time.time()

def get_ref_int(ind):
    ref_subseq = get_ref_seq(ind)
    ref_int = 0
    for i in range(query_size):
        ref_int = ref_int << 2 | nucleo[ref_subseq[i]]
    return ref_int

def get_ref_seq(ind):
    if suffix_array[ind] + query_size > len(ref_seq):
        return None
    return ref_seq[suffix_array[ind] - 1:(suffix_array[ind] + query_size) - 1]

def binary_search(query_str, lower, upper, strict): # strict false means look for lower bound
    if lower == upper:
        return lower
    if upper - lower == 1:
        l_seq = get_ref_seq(lower)
        u_seq = get_ref_seq(upper)
        if strict:
            if u_seq == query_str:
                return upper
            else:
                return lower
        else:
            if l_seq == query_str:
                return lower
            else:
                return upper
    mid = (lower + upper) // 2
    mid_seq = get_ref_seq(mid)
    while mid_seq is None and mid >= lower:
        mid -= 1
        mid_seq = get_ref_seq(mid)
    if mid_seq < query_str or (mid_seq == query_str and strict):
        return binary_search(query_str, mid, upper, strict)
    else:
        return binary_search(query_str, lower, mid, strict)


def exponential_search(query_str, start_sa):
    lower_bound = None
    upper_bound = None
    curr_seq = get_ref_seq(start_sa)
    while curr_seq is None:
        start_sa += 1
        curr_seq = get_ref_seq(start_sa)

    print(curr_seq)
    print(query_str)

    if curr_seq < query_str:
        lower_bound = start_sa
    elif curr_seq > query_str:
        upper_bound = start_sa

    print(upper_bound)
    print(lower_bound)

    window_size = 1

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
        upper_bound = len(suffix_array) - 1
    return (binary_search(query_str, lower_bound, upper_bound, False), binary_search(query_str, lower_bound, upper_bound, True))


search_fn = exponential_search

pred_dists = []
with open(output_file, "w") as output_file:
    for i in range(0, len(queries)):
        start_sa = y_predict[i]
        print(start_sa)
        (l, u) = search_fn(queries[i], int(start_sa))
        pred_dists.append(max(l - int(start_sa), int(start_sa) - u, 0))
        output_file.write(str(l))
        output_file.write(str((u, len(ref_seq))))

print(l, u)
search_time = time.time()

print("Average prediction distance: ", sum(pred_dists)/len(pred_dists))
print("Maximum prediction distance: ", max(pred_dists))
print("I/O Read Time (s): ", read_time - start_time)
print("Fit Time (s): ", fit_time - read_time)
print("Predict Time (s): ", predict_time - fit_time)
print("Search Time (s): ", search_time - predict_time)

# (9133, 9172)
