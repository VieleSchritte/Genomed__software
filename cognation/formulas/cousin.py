from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException

class CousinFormula(Formula):
    def calculate_relation(self, raw_values):
        lr = 0
        if len(raw_values) < 3:
            raise LineFormatException()

        cousin2_alleles = self.split_sat(raw_values.pop())
        cousin1_alleles = self.split_sat(raw_values.pop())
        locus = ''.join(raw_values)

        if len(cousin1_alleles) != 2 or len(cousin2_alleles) != 2:
            raise AllelesException()

        cousin1_set = set(cousin1_alleles)
        cousin2_set = set(cousin2_alleles)
        intersection_list = list(cousin1_set & cousin2_set)
        freq_dict = self.get_frequencies(locus, intersection_list)

        # case1 - no common alleles
        if len(intersection_list) == 0:
            lr = 0.75

        # case2 - two common alleles, both are heterozygous
        elif len(intersection_list) == 2:
            freq1 = freq_dict[intersection_list[0]]
            freq2 = freq_dict[intersection_list[1]]
            lr = 0.75 + 0.125*(freq1 + freq2)/(2 * freq1 * freq2)

        # only one intersection - one common allele
        elif len(intersection_list) == 1:
            freq = freq_dict[intersection_list[0]]

            # case 3 - both are homozygous
            if len(cousin1_set) == len(cousin2_set) == 1:
                lr = 0.75 + 0.25/freq

            # case4 - one of them is heterozygous
            elif len(cousin1_set) == 1 and len(cousin2_set) == 2:
                lr = 0.75 + 0.125/freq

            # case 5 - both are the same heterozyhous
            elif len(cousin1_set) == len(cousin2_set) == 2:
                lr = 0.75 + 0.125/(2 * freq)

            return self.make_result(locus, '/'.join(cousin1_alleles), '/'.join(cousin2_alleles), lr)
