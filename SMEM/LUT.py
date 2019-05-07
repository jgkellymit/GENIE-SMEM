from ExactMatch import ExactMatch

class LUT:

    def __init__(self, matcher: ExactMatch):
        self.matcher = matcher
        self.LUT = None

    def generate_lut(self, size):
        self.LUT = {}
        ref = self.matcher.ref_sequence[:self.matcher.ref_size-1]

        for pos in range(self.matcher.ref_size):

            for i in range(size):
                if self.matcher.ref_size - size > pos-i >= 0:
                    substring = ref[pos-i:pos+size-i]

                    if substring not in self.LUT:
                        self.LUT[substring] = self.matcher.exact_match_back_prop(substring)


# matcher = ExactMatch("mississippi.fa")
# matcher.create_fm_index()
#lut = LUT(matcher)
#lut.generate_lut(4)
#print("LUT: " + str(lut.LUT))
