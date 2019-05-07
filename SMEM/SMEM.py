from ExactMatch import ExactMatch
from LUT import LUT
import datetime

class SMEM:
    def __init__(self, matcher: ExactMatch):
        self.matcher = matcher

    #@staticmethod
    def get_suffix_index(self, query):
        return self.matcher.exact_match_back_prop(query)

    def get_SMEMS_with_lut(self, query, lut_size, minimum_smem_length):
        lut = LUT(self.matcher)
        lut.generate_lut(lut_size)
        #print("LUT: " + str(lut.LUT))
        ref = self.matcher.ref_sequence[:self.matcher.ref_size-1]

        curr_SMEM_start = 0
        curr_sub_start = 0
        smems = {}

        forward_match = None

        #scroll scross the query finding matches in the lut
        #expand the match backwards when it no longer matches
        while curr_sub_start + lut_size <= len(query):
            sub = query[curr_sub_start: curr_sub_start+lut_size]

            #check if substring is in lut
            if sub in lut.LUT:
                if forward_match == None:
                    forward_match = sub
                else:
                    #TODO: add a check to make sure that the substrings are in sequence
                    #in the reference. This should be done by adding the start index of a substring
                    #to the lut so that we can check it in O(1).
                    forward_match += sub[-1]

                if curr_sub_start + lut_size == len(query):
                    smems[forward_match] = self.get_suffix_index(forward_match)

                curr_sub_start += 1

            else:
                #if there is no forward match yet you cant back extend, so continue
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

                    #move pointers and reset forward match
                    curr_sub_start += 1
                    curr_SMEM_start = curr_sub_start

                    forward_match = None

        return smems

    def backward_extension(self, query, start_index, forward_matches):
        backward_matches = {}
        for key in forward_matches:

            for i in range(start_index-1, -1, -1):

                currentSearch = query[i: start_index] + key

                suffix_tuple = self.get_suffix_index(currentSearch)

                if suffix_tuple == -1:
                    break
                else:
                    backward_matches[currentSearch] = suffix_tuple
        return backward_matches


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
        forward_matches = {}

        #forward extend
        for i in range(start_index+1, len(query)+1):
            currentSearch = query[start_index:i]

            suffix_tuple = self.get_suffix_index(currentSearch)

            if suffix_tuple == -1:
                break
            else:
                forward_matches[currentSearch] = suffix_tuple

        #print("Forward Extension Matches: " + str(forward_matches))

        backward_matches = self.backward_extension(query, start_index, forward_matches)

        #print("Backward Extension Matches: " + str(backward_matches))

        #get SMEM from matches
        largest = ''

        if len(backward_matches) == 0:
            for match in forward_matches:
                if len(match) > len(largest):
                    largest = match
            return [largest, forward_matches[largest]]
        else:
            for match in backward_matches:
                if len(match) > len(largest):
                    largest = match

            return [largest, backward_matches[largest]]



if __name__ == '__main__':
    matcher = ExactMatch("full_data.fa")
    matcher.create_fm_index()
    query = matcher.create_query(10000)
    #query = ""
    print("Query Sequence: " + query)
    smem = SMEM(matcher)

    time1 = datetime.datetime.now()
    smem.get_SMEMS(query, 3)
    time2 = datetime.datetime.now()
    smem.get_SMEMS_with_lut(query, 1000, None)
    end = datetime.datetime.now()

    print("No LUT time: " + str(time2-time1))
    print("LUT time: " +  str(end - time2))
