from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException

class ParentFormula(Formula):
    def calculate_relation(self, raw_values):
        if len(raw_values) < 3:
            raise LineFormatException()

        parent_alleles = self.split_sat(raw_values.pop())
        child_alleles = self.split_sat(raw_values.pop())
        locus = ' '.join(raw_values)  # for loci names contain space

        if len(parent_alleles) != 2 or len(child_alleles) != 2:
            raise AllelesException()

        parent_set = set(parent_alleles)  # unique parent alleles
        child_set = set(child_alleles)  # unique child alleles

        relation = 0
        freq_dict = self.get_frequencies(locus, list(parent_set & child_set))
        # i, j     i, j
        if len(freq_dict) == 2:
            divider = 4
            dividend = 0
            for key in freq_dict:
                divider = divider * freq_dict[key]
                dividend = dividend + freq_dict[key]

            if divider != 0:
                # (f1 + f2) / 4 * f1 * f2
                relation = dividend / divider
        elif len(freq_dict) == 1:
            # only 1 allele in common
            # 1 / {1,2}*{1,2}*f
            freq = list(freq_dict.values())[0]
            relation = 1 / (len(parent_set)*len(child_set) * freq)

        return self.make_result(locus, '/'.join(parent_alleles), '/'.join(child_alleles), relation)
