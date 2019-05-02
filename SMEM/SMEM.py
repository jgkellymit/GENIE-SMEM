from ExactMatch import ExactMatch


class SMEM:
    def __init__(self, matcher: ExactMatch):
        self.matcher = matcher

    def get_suffix_index(self, query):
        return self.matcher.exact_match_back_prop(query)


    def get_SMEMS(self, query):


        currentIndex = 0
        smem_table = {}

        while currentIndex <= len(query):
            smem = self.get_SMEM_at_index(query, currentIndex)


            sequence = smem[0]
            start_index_suffix = smem[1][0]
            end_index_suffix = smem[1][1]

            for i in range(start_index_suffix, end_index_suffix):
            if i in smem_table:
                if len(sequence) > len(smem_table[i]):
                    smem_table[i] = sequence

            currentIndex = end_index_suffix + 1

        return smem_table


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

        #get SMEM from backward_matches
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
    query = matcher.create_query(2)
    print(query)
    smem = SMEM(matcher)
    print(smem.get_SMEMS(query))
