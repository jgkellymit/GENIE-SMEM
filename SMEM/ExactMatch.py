from os import path
import json
import random



class ExactMatch:

    def __init__(self, reference_sequence_file: str, query_sequence_file: str=None):
        self.ref_seq_file = reference_sequence_file

        self.ref_sequence = None
        self.ref_size = None
        self.query_sequence = None
        self.fm_file = self.ref_seq_file.split(".")[0] + "-FM.json"
        self.fm_index = {}

        if query_sequence_file is not None:
            self.load_query(query_sequence_file)


    def create_fm_index(self):
        with open(path.join("data", self.ref_seq_file), "r") as refFile:
            refFile.readline()
            self.ref_sequence = ""
            for line in refFile:
                self.ref_sequence += line.strip()

        bwt_array, first_index_array, suffix_array = self.create_bwt_matrix()

        o_matrix = self.get_occurance_matrix(bwt_array)
        c_dic = self.get_count_dic(first_index_array)
        self.fm_index = {"bwt_array": bwt_array, "suffix_array": suffix_array, "occurance_matrix": o_matrix,
                         "count_dic": c_dic}

        with open(path.join("data", self.fm_file), "w") as fm_file:
            fm_file.write(json.dumps(self.fm_index, indent=4, sort_keys=True))


    def create_bwt_matrix(self):
        self.ref_sequence += "$"
        self.ref_size = len(self.ref_sequence)
        bwt_matrix = []

        # Rotate
        for index in range(self.ref_size):
            bwt_matrix.append(self.ref_sequence[index:] + self.ref_sequence[:index])
        bwt_matrix.sort()

        bwt_array = []
        first_index_array = []
        suffix_array = []
        for rotation in bwt_matrix:
            bwt_array.append(rotation[-1])
            first_index_array.append(rotation[0])
            suffix_array.append(self.ref_size - rotation.index("$"))

        return bwt_array, first_index_array, suffix_array

    @staticmethod
    def get_occurance_matrix(bwt_array):
        o_dic = {}
        index = 0
        for item in bwt_array:
            if item not in o_dic:
                o_dic[item] = [0]
            item_list = o_dic[item]
            while len(item_list) <= index:
                item_list.append(item_list[-1])
            item_list[-1] += 1

            index += 1

        # Extend everything to full matrix length
        for key in o_dic:
            item_list = o_dic[key]
            while len(item_list) < index:
                item_list.append(item_list[-1])

        return o_dic

    @staticmethod
    def get_count_dic(first_index_array):
        count_dic = {}
        index = 0
        for char in first_index_array:
            if char not in count_dic:
                count_dic[char] = index
            index += 1
        count_dic[""] = index
        return count_dic


    def load_query(self, query_seq_file):
        with open(path.join("data", query_seq_file), "r") as q_file:
            self.query_sequence = ""
            for line in q_file:
                self.query_sequence += line.strip()


    # Create a query of size "query_size" from the reference sequence
    def create_query(self, query_size, query_output_file=None):
        if self.ref_sequence is None:
            with open(path.join("data", self.ref_seq_file), "r") as refFile:
                refFile.readline()
                self.ref_sequence = ""
                for line in refFile:
                    self.ref_sequence += line.strip()
            self.ref_size = len(self.ref_sequence)

        random.seed(0)
        query_start = random.randint(0, self.ref_size - query_size)
        query = self.ref_sequence[query_start: query_start + query_size]
        if query_output_file is not None:
            with open(path.join("data", query_output_file), "w+") as q_out_f:
                q_out_f.write(">query from '" + self.ref_seq_file + "' zero-index location: " + str(query_start) + "\n")
                q_out_f.write(query + "\n")

        self.query_sequence = query
        return query


    # Do exact match back prop and return the locations in the suffix array of matches
    def exact_match_back_prop(self, query_seq: str):
        if self.fm_index == {}:
            try:
                with open(path.join("data", self.fm_file), "r") as FM_file:
                    self.fm_index = json.load(FM_file)
            except FileNotFoundError:
                print("No FM index file found. Run ExactMatch.createFMIndex to create an FM index.")
                return

        start = 1
        end = self.fm_index["count_dic"][""]
        for char in reversed(query_seq):
            char_count = self.fm_index["count_dic"][char]
            if start - 1 <= 0:
                start = char_count + 1
            else:
                start = char_count + 1 + self.fm_index["occurance_matrix"][char][start - 2]  # minus 2 bc zero index

            end = char_count + self.fm_index["occurance_matrix"][char][end - 1]  # minus 1 bc zero index

            # No match in ref seq
            if start > end:
                return -1

        return start - 1, end - 1  # return zero indexed start and end


    def exact_match(self, query_seq: str=None):
        if query_seq is None:
            if self.query_sequence is None:
                print("No query sequence. Either input sequence, load_file, or generate.")
                return
            query_seq = self.query_sequence

        start, end = self.exact_match_back_prop(query_seq)

        matches = []
        for index in range(start, end + 1):
            matches.append(self.fm_index["suffix_array"][index])

        matches.sort()
        return matches

    def get_position(self, suffix_array_index):
        return self.fm_index["suffix_array"][suffix_array_index]

# e_match = ExactMatch("mississippi.fa")
# e_match.create_fm_index()
# e_match.create_query(3, "query1.fa")
# print(e_match.exact_match())
