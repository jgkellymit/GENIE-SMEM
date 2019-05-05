from ExactMatch import ExactMatch

class LUT:

    def __init__(self, matcher: ExactMatch):
        self.matcher = matcher
        self.LUT = None

    def generate_lut(self, size):
        self.LUT = {}
        ref = self.matcher.ref_sequence[:self.matcher.ref_size-1]

        for pos in range(self.matcher.ref_size):
            substrings = []

            for i in range(size):
                if self.matcher.ref_size - size > pos-i >= 0:
                    substrings.append(ref[pos-i:pos+size-i])

            if len(substrings) > 0:
                self.LUT[pos] = substrings


matcher = ExactMatch("mississippi.fa")
matcher.create_fm_index()
lut = LUT(matcher)
lut.generate_lut(4)

print(lut.LUT)
