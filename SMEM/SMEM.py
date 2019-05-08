from ExactMatch import ExactMatch
from LUT import LUT
import datetime

class SMEM:
    def __init__(self, matcher: ExactMatch):
        self.matcher = matcher
        self.lut = LUT(self.matcher)
        self.lut.load_lut()

    #@staticmethod
    def get_suffix_index(self, query):
        return self.matcher.exact_match_back_prop(query)

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
                    if len(backward_match) == 0:
                        smems[forward_match] = self.get_suffix_index(forward_match)
                    else:
                        smems[backward_match] = self.get_suffix_index(backward_match)

                    # move pointers and reset forward match
                    curr_sub_start += 1
                    curr_SMEM_start = curr_sub_start

                    forward_match = None

        return smems



    def get_SMEMS_with_lut_simple(self, query, minimum_smem_length= None):
        smems = {}

        start_index = 0
        end_index = start_index + self.lut.lut_size

        while end_index <= len(query):
            sub = query[start_index: end_index]
            if sub in self.lut.lut:
                #forward extend from end of sub
                suffix_tuple = self.get_suffix_index(sub)
                forward_matches = {sub: suffix_tuple}
                largest_forward = ''
                currentSearch = sub

                for i in range(end_index+1, len(query)+1):
                    currentSearch += query[i]
                    suffix_tuple = self.get_suffix_index(currentSearch)

                    if suffix_tuple == -1:
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
                for i in range(self.lut.lut_size):
                    smem = self.get_SMEM_at_index(query, start_index)

                    if len(smem[0]) >= minimum_smem_length:
                        smems[smem[0]] = smem[1]

                #move pointers
                start_index = end_index + 1
                end_index = start_index + self.lut.lut_size

        return smems


    def backward_extension(self, query, start_index, forward_matches):
        largest = ''
        suffix_of_largest = None

        for key in forward_matches:

            for i in range(start_index-1, -1, -1):

                currentSearch = query[i: start_index] + key

                suffix_tuple = self.get_suffix_index(currentSearch)

                if suffix_tuple == -1:
                    break
                else:
                    if len(currentSearch) > len(largest):
                        largest = currentSearch
                        suffix_of_largest = suffix_tuple

        return largest, suffix_of_largest

    def forward_extension(self, query, start_index):
        forward_matches = {}
        largest = ''

        for i in range(start_index+1, len(query)+1):
            currentSearch = query[start_index:i]

            # TODO -- Faster to narrow down using the suffix indices?
            suffix_tuple = self.get_suffix_index(currentSearch)

            if suffix_tuple == -1:
                break
            else:
                forward_matches[currentSearch] = suffix_tuple
                if len(currentSearch) > len(largest):
                    largest = currentSearch

        return forward_matches, largest


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
                smems[currentIndex] = smem
            currentIndex += len(smem[0])

        return smems

    def get_SMEM_at_index(self, query, start_index):

        # Forward extend
        forward_extension = self.forward_extension(query, start_index)

        #print("Forward Extension Matches: " + str(forward_matches))

        largest_backward = self.backward_extension(query, start_index, forward_extension[0])

        #print("Backward Extension Matches: " + str(backward_matches))

        # Get SMEM from matches
        if len(forward_extension[1]) > len(largest_backward[0]):
            return [forward_extension[1], forward_extension[0][forward_extension[1]]]
        else:
            return [largest_backward[0], largest_backward[1]]



if __name__ == '__main__':
    matcher = ExactMatch("medium_data.fa")
    matcher.load_fm_index()

    #query = matcher.create_query(100)
    query="CCCAACCCCAACCCTAAGGGTGGAGGCGTGGCGCAGGCGCAGAG"
    print("Query Sequence: " + query)
    smem = SMEM(matcher)

    # start = datetime.datetime.now()
    print(smem.get_SMEMS(query, 2))
    # end1 = datetime.datetime.now()
    # print(smem.get_SMEMS_with_lut(query, None))
    # end2 = datetime.datetime.now()

    print(smem.get_SMEMS_with_lut_simple(query, 2))

    # print("Standard SMEM Time: " + str(end1-start))
    # print("LUT SMEM Time: " + str(end2-end1))
