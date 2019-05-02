from ExactMatch import ExactMatch


class SMEM:
    def __init__(self, matcher: ExactMatch):
        self.matcher = matcher

    def get_suffix_index(self, query):
        return self.matcher.exact_match_back_prop(query)

    """
        Get the SMEMs for every position in the query seqeunce

        Args: query (str)
        Returns: dictionary, keys are the start index of an SMEM in
                 the query, values are the SMEM in the form:
                 [sequence (str), (start_suffix_index (int)), end_suffix_index(int))]
    """
    def get_SMEMS(self, query):

        currentIndex = 0
        smems = {}

        while currentIndex <= len(query)-1:
            smem = self.get_SMEM_at_index(query, currentIndex)
            smems[currentIndex] = smem
            currentIndex += len(smem[0])
        return smems

    def get_SMEM_at_index(self, query, start_index):
        forward_matches = {}
        backward_matches = {}

        #forward extend
        for i in range(start_index+1, len(query)+1):
            currentSearch = query[start_index:i]

            suffix_tuple = self.get_suffix_index(currentSearch)

            if suffix_tuple == -1:
                break
            else:
                forward_matches[currentSearch] = suffix_tuple

        print("Forward Extension Matches: " + str(forward_matches))

        #backward extend
        for key in forward_matches:

            for i in range(start_index-1, -1, -1):

                currentSearch = query[i: start_index] + key

                suffix_tuple = self.get_suffix_index(currentSearch)

                if suffix_tuple == -1:
                    break
                else:
                    backward_matches[currentSearch] = suffix_tuple

        print("Backward Extension Matches: " + str(backward_matches))

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
    matcher = ExactMatch("mississippi.fa")
    #query = matcher.create_query(3)
    query = "ipspspi"
    print("Query Sequence: " + query)
    smem = SMEM(matcher)
    print("SMEM Results: " + str(smem.get_SMEMS(query)))
