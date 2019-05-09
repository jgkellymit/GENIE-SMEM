from ExactMatch import ExactMatch
from LUT import LUT
import datetime

class SMEM:
    def __init__(self, matcher: ExactMatch):
        self.matcher = matcher
        self.lut = LUT(self.matcher)
        self.lut.load_lut()

    def get_suffix_index(self, query):
        return self.matcher.exact_match_back_prop(query)



    def get_smems_lut(self, query):

        all_smems = {}
        current_index = 0

        # First smem discovery
        current_sub = query[:self.lut.lut_size]

        # LUT match found -> forward extend from end of match
        if current_sub in self.lut.lut:
            forward_match = self.forward_extension(query, current_index + self.lut.lut_size, current_sub, self.lut.lut[current_sub][0])

        # LUT not useful for beginning -> do forward extension from start
        else:
            forward_match = self.forward_extension(query, current_index)

        first_smem = forward_match[1]
        all_smems[first_smem] = forward_match[0][first_smem]


        # Rest of smems (starting point in middle)
        current_index += len(first_smem)

        # Frame syntax: (<string_seq>, <suffix_tuple>, <start_index>, <backward_extend>, <forward_extend>)


        end_smem_index = current_index
        while end_smem_index < len(query):
            prev_frame = None
            current_smem_possibility = None
            current_smem_suffix = ()
            current_smem_end_index = -1

            for i in range(self.lut.lut_size):
                current_index = end_smem_index - i
                if current_index + self.lut.lut_size > len(query):
                    continue

                current_sub = query[current_index: current_index + self.lut.lut_size]
                if current_sub in self.lut.lut:  # Case 1,2,4
                    if prev_frame is None:  # First check
                        # print("First with lut match")
                        prev_frame = (current_sub, self.lut.lut[current_sub][0], current_index, True, True)
                    elif prev_frame == ():  # Case 4
                        # print("Case 4")
                        prev_frame = (current_sub, self.lut.lut[current_sub][0], current_index, True, False)
                    else:
                        if self.check_sequential(self.lut.lut[current_sub][1], self.lut.lut[prev_frame[0]][1]):  # Case 1
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
                                backward_smems = self.backward_extension(query, prev_frame[2], {prev_frame[0]: prev_frame[1]})
                                if current_smem_possibility is None or len(backward_smems[0]) >= len(current_smem_possibility):
                                    current_smem_possibility = backward_smems[0]
                                    current_smem_suffix = backward_smems[1]
                                    current_smem_end_index = backward_smems[2]

                            else:
                                raise Exception("previous frame should either forward or backward extend")

                            prev_frame = (current_sub, self.lut.lut[current_sub][0], current_index, True, False)

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
                                    current_smem_suffix = self.lut.lut[current_sub][0]
                                    current_smem_end_index = self.lut.lut_size + current_index

                            prev_frame = (current_sub, self.lut.lut[current_sub][0], current_index, True, False)
                else:
                    if prev_frame is None or prev_frame == ():  # First check
                        # print("First no lut match")
                        prev_frame = ()
                    else:  # Case 3
                        # print("Case 3")
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
                # print("No Match")

            else:
                end_smem_index = current_smem_end_index
                all_smems[current_smem_possibility] = current_smem_suffix
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







    def get_SMEMS_with_lut(self, query, minimum_smem_length):

        ref = self.matcher.ref_sequence[:-1]

        curr_SMEM_start = 0
        curr_sub_start = 0
        smems = {}

        forward_match = None
        smem_prev_indices = None
        # scroll scross the query finding matches in the lut
        # expand the match backwards when it no longer matches
        while curr_sub_start + self.lut.lut_size <= len(query):
            sub = query[curr_sub_start: curr_sub_start+self.lut.lut_size]

            # check if substring is in lut
            if sub in self.lut.lut:
                if forward_match == None:
                    forward_match = sub

                else:
                    # TODO: add a check to make sure that the substrings are in sequence
                    # in the reference. This should be done by adding the start index of a substring
                    # to the lut so that we can check it in O(1).
                    forward_match += sub[-1]

                if curr_sub_start + self.lut.lut_size == len(query):
                    smems[forward_match] = self.get_suffix_index(forward_match)

                curr_sub_start += 1

            else:
                # if there is no forward match yet you cant back extend, so continue
                if forward_match == None:
                    curr_sub_start += 1
                    curr_SMEM_start += 1
                else:
                    #back extend the forward match and then add it as an smem
                    backward_match = self.backward_extension(query, curr_SMEM_start, {forward_match:None})
                    if len(backward_match[0]) == 0:
                        smems[forward_match] = self.get_suffix_index(forward_match)
                    else:
                        smems[backward_match] = backward_match[1]

                    # move pointers and reset forward match
                    curr_sub_start += 1
                    curr_SMEM_start = curr_sub_start

                    forward_match = None

        return smems



    def get_SMEMS_with_lut_simple(self, query, minimum_smem_length= None):
        smems = {}

        start_index = 0
        end_index = start_index + self.lut.lut_size

        while start_index < len(query):
            if end_index >= len(query):
                end_index = len(query)

            sub = query[start_index: end_index]
            if sub in self.lut.lut:
                #forward extend from end of sub
                suffix_tuple = self.lut.lut(sub)[0]
                forward_matches = {sub: suffix_tuple}
                largest_forward = sub
                currentSearch = sub

                for i in range(end_index, len(query)):
                    currentSearch += query[i]
                    suffix_tuple = self.get_suffix_index(currentSearch)

                    if suffix_tuple == -1:
                        #end_index = i
                        break
                    else:
                        forward_matches[currentSearch] = suffix_tuple
                        if len(currentSearch) > len(largest_forward):
                            largest_forward = currentSearch

                #backward extend from start of sub
                largest_backward = self.backward_extension(query, start_index, forward_matches)

                #get the largest from the forward and backward extensions
                if len(largest_forward) > len(largest_backward[0]):
                    smems[largest_forward] = forward_matches[largest_forward]
                else:
                    smems[largest_backward[0]] = largest_backward[1]
                #move pointers
                start_index = end_index + 1
                end_index = start_index + self.lut.lut_size

            else:
                #standard SMEM search along indices in the substring that didnt appear in the LUT
                current_index = start_index
                while current_index < end_index:
                    smem = self.get_SMEM_at_index(query, current_index)

                    if len(smem[0]) >= minimum_smem_length:
                        smems[smem[0]] = smem[1]

                    current_index += 1

                #move pointers
                start_index = end_index + 1
                end_index = start_index + self.lut.lut_size

        return smems


    def backward_extension(self, query, start_index, forward_matches):
        largest = ""
        suffix_of_largest = None
        end_index = -1

        largest_forward = ""

        for key in forward_matches:
            if len(key) > len(largest_forward):
                largest_forward = key

            for i in range(start_index-1, -1, -1):

                currentSearch = query[i: start_index] + key

                # TODO -- this can be sped way up since we already have the previous suffix indices
                suffix_tuple = self.get_suffix_index(currentSearch)

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
    import random
    q = ""
    for x in range(query_size):
        q += random.choice(["A", "G", "C", "T"])
    return q

if __name__ == '__main__':
    match = ExactMatch("medium_data.fa")
    match.load_fm_index()

    # q = "CCTAACCCTAACCGCGCAGGCGCAGAGACACATGCTGGCTGTGATACGTGGCCGGCCC" \
    #     "TCGCTCCAGCAGCTGGACCCCAGAACACAGTGGCGCAGGCTGGGTGGAGCCGTCCCCCCTTGTGCCACTTCTGGATGCTAGGGTTACACTG"

    q = create_random_query(200)
    # q = "TTGACGTGGCTATTGTTTCGTCGCAGATCTTGCGCGAGGGAACACGTGCCGTGCCAAAACCAATGGCCCATGTTTAAGCAAGGAACCCTCCGCAGTCCCGCGCTCACTAAGTGTGTATGGTGAAGTGACATGCAAGGATCGGGTGAATTAGTGTTCATTCTCAAGTGTTCCCAGTGGATACTCTTCGTTCCTTCACTGAGCAACGAAGCTACACCCGGAGTTGCGCTTAGGTTACGAAACTGCATGTATGTATCTGCTCACAGCAAAGGTGCCATTTGCCGCATATCACTCTTGCTATGC"

    # q = "TCGGCTGTCGCGGCGCTTGAAGCAAGAGCCGGGTGTAGCTCGTCTAATTGAAGGTAGAGTTCATTGATTAGTCGTGGATTAGCTTATAACGAGTGTGCGA"
    # q = "GTTGTATGTGTACCATCGTT"
    print(q)

    smem = SMEM(match)


    start = datetime.datetime.now()
    a = smem.get_smems_lut(q)
    mid = datetime.datetime.now()

    b = smem.get_SMEMS(q, 1)

    end = datetime.datetime.now()
    print("Lut time:")
    print(mid - start)

    print("OG time")
    print(end -  mid)

    print(a)
    print(b)

    print(len(a))
    print(len(b))
    if len(a) != len(b):
        print("SMEMS DONT MATCH")

    #
    # query = match.create_query(100)


    # # query="GCGCAGGCGCAGGGGCGTGTGTGCCTGTTTCTCCAC"
    # smem = SMEM(match)
    #
    # minimum_smem_length = 10
    # start = datetime.datetime.now()
    # standard = smem.get_SMEMS(query, minimum_smem_length)
    # end1 = datetime.datetime.now()
    # simple_lut = smem.get_SMEMS_with_lut_simple(query, minimum_smem_length)
    # end2 = datetime.datetime.now()
    # lut = smem.get_SMEMS_with_lut(query, minimum_smem_length)
    # end3 = datetime.datetime.now()
    # print("\n")
    # print("Query Sequence: ")
    # print(query)
    # print("___________________________")
    # print("Standard SMEM Result:")
    # print(standard)
    # print("\n")
    # print("Simple LUT SMEM Result:")
    # print(simple_lut)
    # print("\n")
    # print("LUT SMEM Result:")
    # print(lut)
    #
    #
    # print("___________________________")
    # print("Standard SMEM Time: " + str(end1-start))
    # print("Simple LUT SMEM Time: " + str(end3-end2))
    # print("LUT SMEM Time: " + str(end2-end1))
    # print("\n")
