from ExactMatch import ExactMatch
from LUT import LUT
from learned_index.RMI_LUT import RMI_LUT
import datetime
import random
import os

class SMEM:
    def __init__(self, matcher: ExactMatch):
        self.matcher = matcher
        self.lut = LUT(self.matcher)
        self.lut.load_lut()
        self.rmi_lut = None


    def get_suffix_index(self, query):
        return self.matcher.exact_match_back_prop(query)


    def get_smems_lut(self, query):

        all_smems = {}
        current_index = 0

        # First smem discovery
        current_sub = query[:self.lut.lut_size]

        encoded_sub = str(self.lut.convert_seq_to_num(current_sub))

        # LUT match found -> forward extend from end of match
        if encoded_sub in self.lut.lut:
            forward_match = self.forward_extension(query, current_index + self.lut.lut_size, current_sub, self.lut.lut[encoded_sub][0])

        # LUT not useful for beginning -> do forward extension from start
        else:
            forward_match = self.forward_extension(query, current_index)

        first_smem = forward_match[1]
        all_smems[first_smem] = forward_match[0][first_smem]


        # Rest of smems (starting point in middle)
        current_index += len(first_smem)

        # Frame syntax: (<string_seq>, <suffix_tuple>, <start_index>, <backward_extend>, <forward_extend>)

        prev_smem_length = len(first_smem)
        end_smem_index = current_index
        while end_smem_index < len(query):
            prev_frame = None
            current_smem_possibility = None
            current_smem_suffix = ()
            current_smem_end_index = -1
            prev_smem_start = end_smem_index - prev_smem_length

            for i in range(self.lut.lut_size):
                if i >= prev_smem_length:  # Within prev smem, dont need to check
                    continue


                current_index = end_smem_index - i
                if current_index + self.lut.lut_size > len(query):
                    continue

                current_sub = query[current_index: current_index + self.lut.lut_size]
                encoded_sub = str(self.lut.convert_seq_to_num(current_sub))
                if encoded_sub in self.lut.lut:  # Case 1,2,4
                    if prev_frame is None:  # First check
                        # print("First with lut match")
                        prev_frame = (current_sub, self.lut.lut[encoded_sub][0], current_index, True, True)
                    elif prev_frame == ():  # Case 4
                        # print("Case 4")
                        prev_frame = (current_sub, self.lut.lut[encoded_sub][0], current_index, True, False)
                    else:
                        if self.check_sequential(self.lut.lut[encoded_sub][1], self.lut.lut[str(self.lut.convert_seq_to_num(prev_frame[0]))][1]):  # Case 1
                            # print("Case 1")
                            if prev_frame[4]:  # forward extend
                                forward_match = self.forward_extension(query, prev_frame[2] + self.lut.lut_size, prev_frame[0], prev_frame[1])
                                if prev_frame[3]:
                                    backward_smems = self.backward_extension(query, prev_frame[2], forward_match[0])
                                    if current_smem_possibility is None or len(backward_smems[0]) >= len(current_smem_possibility):
                                        current_smem_possibility = backward_smems[0]
                                        current_smem_suffix = backward_smems[1]
                                        current_smem_end_index = backward_smems[2]
                                else:
                                    if forward_match[1] == "":
                                        pass
                                    elif current_smem_possibility is None or len(forward_match[1]) >= len(current_smem_possibility):
                                        current_smem_possibility = forward_match[1]
                                        current_smem_suffix = forward_match[0][forward_match[1]]
                                        current_smem_end_index = len(forward_match[1]) + prev_frame[2]
                            elif prev_frame[3]:  # backward extend
                                # This is a theoretical improvement but doesnt seem to happen in practice
                                if current_smem_possibility is not None and (prev_frame[2] - prev_smem_start) + self.lut.lut_size < len(current_smem_possibility):
                                    continue  # this SMEM cannot be longer than previously found SMEM

                                backward_smems = self.backward_extension(query, prev_frame[2], {prev_frame[0]: prev_frame[1]})
                                if current_smem_possibility is None or len(backward_smems[0]) >= len(current_smem_possibility):
                                    current_smem_possibility = backward_smems[0]
                                    current_smem_suffix = backward_smems[1]
                                    current_smem_end_index = backward_smems[2]

                            else:
                                raise Exception("previous frame should either forward or backward extend")

                            prev_frame = (current_sub, self.lut.lut[encoded_sub][0], current_index, True, False)

                        else:  # Case 2
                            # print("Case 2")
                            if prev_frame[4]:  # forward extend
                                forward_match = self.forward_extension(query, prev_frame[2] + self.lut.lut_size, prev_frame[0], prev_frame[1])
                                if forward_match[1] == "":
                                    pass
                                elif current_smem_possibility is None or len(forward_match[1]) >= len(current_smem_possibility):
                                    current_smem_possibility = forward_match[1]
                                    current_smem_suffix = forward_match[0][forward_match[1]]
                                    current_smem_end_index = len(forward_match[1]) + prev_frame[2]
                            else:
                                if current_smem_possibility is None or len(current_sub) >= len(current_smem_possibility):
                                    current_smem_possibility = current_sub
                                    current_smem_suffix = self.lut.lut[encoded_sub][0]
                                    current_smem_end_index = self.lut.lut_size + current_index

                            prev_frame = (current_sub, self.lut.lut[encoded_sub][0], current_index, True, False)
                else:
                    if prev_frame is None or prev_frame == ():  # First check
                        # print("First no lut match")
                        prev_frame = ()
                    else:  # Case 3
                        if prev_frame[4]:  # forward extend
                            forward_match = self.forward_extension(query, prev_frame[2] + self.lut.lut_size, prev_frame[0],
                                                                   prev_frame[1])
                            if forward_match[1] == "":
                                pass
                            elif current_smem_possibility is None or len(forward_match[1]) >= len(current_smem_possibility):
                                current_smem_possibility = forward_match[1]
                                current_smem_suffix = forward_match[0][forward_match[1]]
                                current_smem_end_index = len(forward_match[1]) + prev_frame[2]

                        else:
                            if current_smem_possibility is None or len(prev_frame[0]) >= len(current_smem_possibility):
                                current_smem_possibility = prev_frame[0]
                                current_smem_suffix = prev_frame[1]
                                current_smem_end_index = self.lut.lut_size + prev_frame[2]

                        prev_frame = ()

            # Last frame extension
            if prev_frame is not None and prev_frame != ():
                if prev_frame[4]:  # forward extend
                    forward_match = self.forward_extension(query, prev_frame[2] + self.lut.lut_size, prev_frame[0],
                                                           prev_frame[1])
                    if prev_frame[3]:
                        backward_smems = self.backward_extension(query, prev_frame[2], forward_match[0])
                        if current_smem_possibility is None or len(backward_smems[0]) >= len(current_smem_possibility):
                            current_smem_possibility = backward_smems[0]
                            current_smem_suffix = backward_smems[1]
                            current_smem_end_index = backward_smems[2]
                    else:
                        if forward_match[1] == "":
                            pass
                        elif current_smem_possibility is None or len(forward_match[1]) >= len(current_smem_possibility):
                            current_smem_possibility = forward_match[1]
                            current_smem_suffix = forward_match[0][forward_match[1]]
                            current_smem_end_index = len(forward_match[1]) + prev_frame[2]
                elif prev_frame[3]:  # backward extend
                    backward_smems = self.backward_extension(query, prev_frame[2], {prev_frame[0]: prev_frame[1]})
                    if current_smem_possibility is None or len(backward_smems[0]) >= len(current_smem_possibility):
                        current_smem_possibility = backward_smems[0]
                        current_smem_suffix = backward_smems[1]
                        current_smem_end_index = backward_smems[2]



            if current_smem_possibility is None:  # No LUT match at any partition
                new_smem = self.get_SMEM_at_index(query, end_smem_index)
                all_smems[new_smem[0]] = new_smem[1]
                end_smem_index = new_smem[2]
                prev_smem_length = len(new_smem[0])

                # print("No Match")

            else:
                end_smem_index = current_smem_end_index
                all_smems[current_smem_possibility] = current_smem_suffix
                prev_smem_length = len(current_smem_possibility)

                # print("Matched")

            # print(end_smem_index)

        return all_smems



    @staticmethod
    def check_sequential(list1, list2):
        for item1 in list1:
            for item2 in list2:
                if item1 + 1 == item2:
                    return True
        return False



    def get_smems_rmi(self, query):
        self.rmi_lut = RMI_LUT.load("learned_index/rmi_file.pkl")

        all_smems = {}
        current_index = 0

        # First smem discovery
        current_sub = query[:self.rmi_lut.prediction_size]

        rmi_prediction = self.rmi_lut.get_suffix_rmi(current_sub, encoded=False)
        # Prediction found -> forward extend from end of match
        if rmi_prediction[1] >= rmi_prediction[0]:
            forward_match = self.forward_extension(query, current_index + self.rmi_lut.prediction_size, current_sub, rmi_prediction)

        # LUT not useful for beginning -> do forward extension from start
        else:
            forward_match = self.forward_extension(query, current_index)

        first_smem = forward_match[1]
        all_smems[first_smem] = forward_match[0][first_smem]


        # Rest of smems (starting point in middle)
        current_index += len(first_smem)

        # Frame syntax: (<string_seq>, <suffix_tuple>, <start_index>, <backward_extend>, <forward_extend>)

        prev_smem_length = len(first_smem)
        end_smem_index = current_index
        while end_smem_index < len(query):
            prev_frame = None
            current_smem_possibility = None
            current_smem_suffix = ()
            current_smem_end_index = -1
            prev_smem_start = end_smem_index - prev_smem_length

            for i in range(self.rmi_lut.prediction_size):
                if i >= prev_smem_length:  # Within prev smem, dont need to check
                    continue

                current_index = end_smem_index - i
                if current_index + self.rmi_lut.prediction_size > len(query):
                    continue

                current_sub = query[current_index: current_index + self.rmi_lut.prediction_size]
                rmi_prediction = self.rmi_lut.get_suffix_rmi(current_sub, encoded=False)

                if rmi_prediction[1] >= rmi_prediction[0]:  # Case 1,2,4
                    if prev_frame is None:  # First check
                        # print("First with lut match")
                        prev_frame = (current_sub, rmi_prediction, current_index, True, True)
                    elif prev_frame == ():  # Case 4
                        # print("Case 4")
                        prev_frame = (current_sub, rmi_prediction, current_index, True, False)
                    else:

                        prediction_real_location = self.matcher.get_positions(rmi_prediction[0], rmi_prediction[1])
                        prev_prediction_real_location = self.matcher.get_positions(prev_frame[1][0], prev_frame[1][1])

                        if self.check_sequential(prediction_real_location, prev_prediction_real_location):  # Case 1
                            # print("Case 1")
                            if prev_frame[4]:  # forward extend
                                forward_match = self.forward_extension(query, prev_frame[2] + self.rmi_lut.prediction_size, prev_frame[0], prev_frame[1])
                                if prev_frame[3]:
                                    backward_smems = self.backward_extension(query, prev_frame[2], forward_match[0])
                                    if current_smem_possibility is None or len(backward_smems[0]) >= len(current_smem_possibility):
                                        current_smem_possibility = backward_smems[0]
                                        current_smem_suffix = backward_smems[1]
                                        current_smem_end_index = backward_smems[2]
                                else:
                                    if forward_match[1] == "":
                                        pass
                                    elif current_smem_possibility is None or len(forward_match[1]) >= len(current_smem_possibility):
                                        current_smem_possibility = forward_match[1]
                                        current_smem_suffix = forward_match[0][forward_match[1]]
                                        current_smem_end_index = len(forward_match[1]) + prev_frame[2]
                            elif prev_frame[3]:  # backward extend

                                # This is a theoretical improvement but doesnt seem to happen in practice
                                if current_smem_possibility is not None and (
                                        prev_frame[2] - prev_smem_start) + self.rmi_lut.prediction_size < len(
                                        current_smem_possibility):
                                    continue  # this SMEM cannot be longer than previously found SMEM

                                backward_smems = self.backward_extension(query, prev_frame[2], {prev_frame[0]: prev_frame[1]})
                                if current_smem_possibility is None or len(backward_smems[0]) >= len(current_smem_possibility):
                                    current_smem_possibility = backward_smems[0]
                                    current_smem_suffix = backward_smems[1]
                                    current_smem_end_index = backward_smems[2]

                            else:
                                raise Exception("previous frame should either forward or backward extend")

                            prev_frame = (current_sub, rmi_prediction, current_index, True, False)

                        else:  # Case 2
                            # print("Case 2")
                            if prev_frame[4]:  # forward extend
                                forward_match = self.forward_extension(query, prev_frame[2] + self.rmi_lut.prediction_size, prev_frame[0], prev_frame[1])
                                if forward_match[1] == "":
                                    pass
                                elif current_smem_possibility is None or len(forward_match[1]) >= len(current_smem_possibility):
                                    current_smem_possibility = forward_match[1]
                                    current_smem_suffix = forward_match[0][forward_match[1]]
                                    current_smem_end_index = len(forward_match[1]) + prev_frame[2]
                            else:
                                if current_smem_possibility is None or len(current_sub) >= len(current_smem_possibility):
                                    current_smem_possibility = current_sub
                                    current_smem_suffix = rmi_prediction
                                    current_smem_end_index = self.rmi_lut.prediction_size + current_index

                            prev_frame = (current_sub, rmi_prediction, current_index, True, False)
                else:
                    if prev_frame is None or prev_frame == ():  # First check
                        # print("First no lut match")
                        prev_frame = ()
                    else:  # Case 3
                        if prev_frame[4]:  # forward extend
                            forward_match = self.forward_extension(query, prev_frame[2] + self.rmi_lut.prediction_size, prev_frame[0],
                                                                   prev_frame[1])
                            if forward_match[1] == "":
                                pass
                            elif current_smem_possibility is None or len(forward_match[1]) >= len(current_smem_possibility):
                                current_smem_possibility = forward_match[1]
                                current_smem_suffix = forward_match[0][forward_match[1]]
                                current_smem_end_index = len(forward_match[1]) + prev_frame[2]

                        else:
                            if current_smem_possibility is None or len(prev_frame[0]) >= len(current_smem_possibility):
                                current_smem_possibility = prev_frame[0]
                                current_smem_suffix = prev_frame[1]
                                current_smem_end_index = self.rmi_lut.prediction_size + prev_frame[2]

                        prev_frame = ()

            # Last frame extension
            if prev_frame is not None and prev_frame != ():
                if prev_frame[4]:  # forward extend
                    forward_match = self.forward_extension(query, prev_frame[2] + self.rmi_lut.prediction_size, prev_frame[0],
                                                           prev_frame[1])
                    if prev_frame[3]:
                        backward_smems = self.backward_extension(query, prev_frame[2], forward_match[0])
                        if current_smem_possibility is None or len(backward_smems[0]) >= len(current_smem_possibility):
                            current_smem_possibility = backward_smems[0]
                            current_smem_suffix = backward_smems[1]
                            current_smem_end_index = backward_smems[2]
                    else:
                        if forward_match[1] == "":
                            pass
                        elif current_smem_possibility is None or len(forward_match[1]) >= len(current_smem_possibility):
                            current_smem_possibility = forward_match[1]
                            current_smem_suffix = forward_match[0][forward_match[1]]
                            current_smem_end_index = len(forward_match[1]) + prev_frame[2]
                elif prev_frame[3]:  # backward extend
                    backward_smems = self.backward_extension(query, prev_frame[2], {prev_frame[0]: prev_frame[1]})
                    if current_smem_possibility is None or len(backward_smems[0]) >= len(current_smem_possibility):
                        current_smem_possibility = backward_smems[0]
                        current_smem_suffix = backward_smems[1]
                        current_smem_end_index = backward_smems[2]



            if current_smem_possibility is None:  # No LUT match at any partition
                new_smem = self.get_SMEM_at_index(query, end_smem_index)
                all_smems[new_smem[0]] = new_smem[1]
                end_smem_index = new_smem[2]
                prev_smem_length = len(new_smem[0])
                # print("No Match SMEM: " + new_smem[0])


            else:
                end_smem_index = current_smem_end_index
                all_smems[current_smem_possibility] = current_smem_suffix
                prev_smem_length = len(current_smem_possibility)
                # print("Matched SMEM: " + current_smem_possibility)

            # print(end_smem_index)

        return all_smems




    def backward_extension(self, query, start_index, forward_matches):
        largest = ""
        suffix_of_largest = None
        end_index = -1

        largest_forward = ""

        for key in forward_matches:
            if len(key) > len(largest_forward):
                largest_forward = key

            suffix_tuple = None

            for i in range(start_index-1, -1, -1):

                currentSearch = query[i: start_index] + key
                if suffix_tuple is None:
                    suffix_tuple = self.get_suffix_index(currentSearch)
                else:
                    suffix_tuple = self.matcher.exact_match_back_prop_add_one(currentSearch[0], suffix_tuple)

                if suffix_tuple == -1:
                    break
                else:
                    if len(currentSearch) > len(largest):
                        largest = currentSearch
                        suffix_of_largest = suffix_tuple
                        end_index = start_index + len(key)

        if len(largest_forward) > len(largest):
            largest = largest_forward
            suffix_of_largest = forward_matches[largest_forward]
            end_index = start_index + len(largest_forward)

        return largest, suffix_of_largest, end_index

    def forward_extension(self, query, start_index, largest="", suffix_tuple=None):
        forward_matches = {}
        if suffix_tuple is not None:
            forward_matches[largest] = suffix_tuple

        currentSearch = largest
        for i in range(start_index+1, len(query)+1):
            currentSearch = largest + query[start_index:i]

            # TODO -- Faster to narrow down using the suffix indices?
            suffix_tuple = self.get_suffix_index(currentSearch)

            if suffix_tuple == -1:
                return forward_matches, currentSearch[:-1]  # Don't include last index
            else:
                forward_matches[currentSearch] = suffix_tuple

        # Hit end of query
        return forward_matches, currentSearch


    """
        Get the SMEMs for every position in the query seqeunce

        Args: query (str)
              minimum_length (int) - minimum length of the SMEMs

        Returns: dictionary, keys are the start index of an SMEM in
                 the query, values are the SMEM in the form:
                 [sequence (str), (start_suffix_index (int)), end_suffix_index(int))]
    """
    def get_SMEMS(self, query, minimum_length):

        currentIndex = 0
        smems = {}

        while currentIndex < len(query):
            smem = self.get_SMEM_at_index(query, currentIndex)
            if len(smem[0]) >= minimum_length:
                smems[smem[0]] = smem[1]
            currentIndex = smem[2]

        return smems

    def get_SMEM_at_index(self, query, start_index):

        # Forward extend
        forward_extension = self.forward_extension(query, start_index)

        #print("Forward Extension Matches: " + str(forward_matches))

        largest_backward = self.backward_extension(query, start_index, forward_extension[0])

        #print("Backward Extension Matches: " + str(backward_matches))

        # Get SMEM from matches
        if len(forward_extension[1]) > len(largest_backward[0]):
            return [forward_extension[1], forward_extension[0][forward_extension[1]], len(forward_extension[1]) + start_index]
        else:
            return [largest_backward[0], largest_backward[1], largest_backward[2]]




def create_random_query(query_size):
    q = ""
    for x in range(query_size):
        q += random.choice(["A", "G", "C", "T"])
    return q


def create_query_from_ref(ref_seq, query_size):
    ref_size = len(ref_seq)
    query = ""
    while len(query) < query_size:
        position = random.randint(0, ref_size)
        size = random.randint(1, 30)
        if size + position > ref_size:
            continue
        query += ref_seq[position: position + size]
    return query[:query_size]


if __name__ == '__main__':
    match = ExactMatch("medium_data.fa")
    match.load_fm_index()
    match.load_ref_sequence()

    # q = create_random_query(200)
    q = create_query_from_ref(match.ref_sequence, 200)
    print(q)
    # q = "CGCAATACCTTAATCGTTTAGTCTTTCATTCCAACTGGAAAGAATTCCCGCACTTCTAGACGCGAGTTCTGCAAAGCCTCAGTTTTACTTCGTATTGCCCGTAGGCAAACGTCCTCTCATTTACCTTCAGTGGCAATGTCTGTGTCAACTGCTACACTTACCTCGATGTAAGAGATGACTACGGCTCTGGAAGATACCCT"

    smem = SMEM(match)


    start = datetime.datetime.now()
    a = smem.get_smems_lut(q)
    mid = datetime.datetime.now()

    b = smem.get_SMEMS(q, 1)

    sec = datetime.datetime.now()

    c = smem.get_smems_rmi(q)

    end = datetime.datetime.now()
    print("Lut time:")
    print(mid - start)

    print("RMI time")
    print(end - sec)

    print("OG time")
    print(sec - mid)

    # print(a)
    # print(b)
    # print(c)
    #
    print(len(a))
    print(len(b))
    print(len(c))
    # if len(a) != len(b):
    #     print("SMEMS DONT MATCH")
