from RMI import *
from Bio import SeqIO
import numpy as np
import json
import pickle
import os
from ExactMatch import ExactMatch


class RMI_LUT:

    def __init__(self, RMI_structure, LUT_size, data_file):
        self.nucleo = {"A": 0, "C": 1,"G": 2,"T": 3}
        self.structure = RMI_structure
        self.prediction_size = LUT_size
        self.data_file = data_file

        directory = os.getcwd()
        top_dir = directory.split("SMEM")[0]
        data_dir = os.path.join(top_dir, "SMEM/data")


        self.ref_seq = ""
        for record in SeqIO.parse(os.path.join(data_dir, data_file), "fasta"):
            self.ref_seq = record.seq

        fm_file = data_file.split(".")[0] + "-FM.json"
        with open(os.path.join(data_dir, fm_file), 'r') as fm_f:
            fm_index = json.load(fm_f)
            self.suffix_array = fm_index["suffix_array"]
            self.ref_seq_size = fm_index["ref_size"] - 1  # this originally counts the $

        self.rmi = RMI(RMI_structure)


    def train_RMI(self):
        coded_substring = []
        sa_position = []
        for i in range(len(self.ref_seq) + 1):
            if self.suffix_array[i] - 1 + self.prediction_size > self.ref_seq_size:
                continue
            ref_subseq = self.ref_seq[(self.suffix_array[i]) - 1:(self.suffix_array[i] + self.prediction_size) - 1]
            subseq = 0
            for j in range(self.prediction_size):
                subseq = subseq << 2 | self.nucleo[ref_subseq[j]]

            coded_substring.append(subseq)
            sa_position.append(i)

        self.rmi.fit(np.asarray(coded_substring).reshape(-1, 1), np.asarray(sa_position))


    def rmi_predict(self, query, encoded=False):
        encoded_query = []
        if encoded:
            encoded_query.append(query)
        else:
            query_int = 0
            for j in range(self.prediction_size):
                query_int = query_int << 2 | self.nucleo[query[j]]
            encoded_query.append(query_int)

        return self.rmi.predict(np.asarray(encoded_query).reshape(-1, 1))



    def get_suffix_rmi(self, query, encoded=False):
        # start = datetime.datetime.now()

        start_sa = self.rmi_predict(query, encoded)
        # mid = datetime.datetime.now()
        result = self.exponential_search(query, int(start_sa))
        # end = datetime.datetime.now()

        # predict_time = (mid - start).total_seconds()
        # search_time = (end - mid).total_seconds()

        return result


    def get_ref_int(self, ind):
        ref_subseq = self.get_ref_seq(ind)
        ref_int = 0
        for i in range(self.prediction_size):
            ref_int = ref_int << 2 | self.nucleo[ref_subseq[i]]
        return ref_int


    def get_ref_seq(self, ind):
        if self.suffix_array[ind] - 1 + self.prediction_size > self.ref_seq_size:
            return None
        return self.ref_seq[self.suffix_array[ind] - 1:(self.suffix_array[ind] + self.prediction_size) - 1]


    def binary_search(self, query_str, lower, upper, strict):  # strict false means look for lower bound
        if lower == upper:
            return lower
        if upper - lower == 1:
            if strict:
                u_seq = self.get_ref_seq(upper)
                if u_seq == query_str:
                    return upper
                else:
                    return lower
            else:
                l_seq = self.get_ref_seq(lower)
                if l_seq == query_str:
                    return lower
                else:
                    return upper
        mid = (lower + upper) // 2
        mid_seq = self.get_ref_seq(mid)
        while mid_seq is None and mid > lower:
            mid -= 1
            mid_seq = self.get_ref_seq(mid)

            if mid == lower:  # avoid potential infinite recursion
                if strict:
                    u_seq = self.get_ref_seq(upper)
                    if u_seq == query_str:
                        return upper
                    else:
                        return lower
                else:
                    if mid_seq == query_str:
                        return lower
                    else:
                        return upper

        if mid_seq < query_str or (mid_seq == query_str and strict):
            return self.binary_search(query_str, mid, upper, strict)
        else:
            return self.binary_search(query_str, lower, mid, strict)


    def exponential_search(self, query_str, start_sa):
        lower_bound = None
        upper_bound = None
        curr_seq = self.get_ref_seq(start_sa)
        while curr_seq is None:
            start_sa += 1
            curr_seq = self.get_ref_seq(start_sa)

        if curr_seq < query_str:
            lower_bound = start_sa
        elif curr_seq > query_str:
            upper_bound = start_sa

        window_size = 1

        if upper_bound is None:
            while start_sa + window_size < len(self.ref_seq) + 1:
                ind = start_sa + window_size
                window_size = 2 * window_size
                found_seq = self.get_ref_seq(ind)
                while found_seq is None:
                    ind += 1
                    found_seq = self.get_ref_seq(ind)
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
                found_seq = self.get_ref_seq(ind)
                while found_seq is None:
                    ind -= 1
                    found_seq = self.get_ref_seq(ind)
                if found_seq < query_str:
                    lower_bound = ind
                    break
                if found_seq > query_str:
                    upper_bound = ind
        if lower_bound is None:
            lower_bound = 0
        if upper_bound is None:
            upper_bound = len(self.suffix_array) - 1
        return (self.binary_search(query_str, lower_bound, upper_bound, False),
                self.binary_search(query_str, lower_bound, upper_bound, True))

    def save(self, file):
        this_model = [self.structure, self.prediction_size, self.data_file, self.rmi]
        with open(file, "wb") as wf:
            pickle.dump(this_model, wf)

    @staticmethod
    def load(file):
        with open(file,"rb") as rf:
            model_params = pickle.load(rf)

        new_rmi = RMI_LUT(model_params[0], model_params[1], model_params[2])
        new_rmi.rmi = model_params[3]
        return new_rmi

    def read_query_and_encode(self, query_db="data/query500.fa"):

        queries = []
        for record in SeqIO.parse(query_db, "fasta"):
            queries.append(record.seq)

        for query in queries:
            query_int = 0
            for j in range(len(query)):
                query_int = query_int << 2 | self.nucleo[query[j]]
        return query_int

if __name__ == '__main__':

    import sys
    import datetime
    # expert = [int(sys.argv[1]), int(sys.argv[2])]
    expert = [10, 100]

    q_size = 15

    rmi_lut = RMI_LUT(expert, q_size, "big_data.fa")
    rmi_lut.train_RMI()
    #
    # e_match = ExactMatch("big_data.fa")
    # e_match.load_fm_index()
    #
    # print("Starting query time analysis")
    # search_total = 0
    # predict_total = 0
    # og_total = 0
    # for x in range(100):
    #     q = e_match.create_query(q_size)
    #
    #     suffix1, predict, search = rmi_lut.get_suffix_rmi(q, False)
    #     start = datetime.datetime.now()
    #     suffix2 = e_match.exact_match_back_prop(q)
    #     end = datetime.datetime.now()
    #
    #     og_timed = end - start
    #
    #     predict_total += predict
    #     search_total += search
    #     og_total += og_timed.total_seconds()
    #
    #     if suffix1 != suffix2:
    #         print(suffix1)
    #         print(suffix2)
    #         raise Exception
    #
    # print(predict_total)
    # print(search_total)
    # print(predict_total + search_total)
    #
    # print(og_total)

    rmi_lut.save("rmi_file.pkl")

    # new_rmi = RMI_LUT.load("rmi_file.pkl")

    # print(new_rmi.get_suffix_rmi("AAAACCACTA"))

