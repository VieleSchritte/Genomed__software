from __future__ import unicode_literals
from .base import Formula


class CousinFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 2)
        cousin1_set, cousin2_set = sets
        intersection = intersections[0]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        freq_dict = self.get_frequencies(locus, intersection)
        lr = 0.75

        # no common alleles
        if len(intersection) == 0:
            return self.make_result(locus, lr, dict_make_result)

        # both are heterozygous, two common alleles
        if len(intersection) == 2:
            freq1 = freq_dict[list(intersection)[0]]
            freq2 = freq_dict[list(intersection)[1]]
            lr += 0.125 * (freq1 + freq2) / (2 * freq1 * freq2)
            return self.make_result(locus, lr, dict_make_result)

        freq = freq_dict[list(intersection)[0]]

        #  both are homozygous, one common allele
        if len(cousin1_set) == len(cousin2_set) == 1:
            lr += 0.25 / freq
            return self.make_result(locus, lr, dict_make_result)

        #  both are heterozygous, one common allele
        if len(cousin1_set) == len(cousin2_set) == 2:
            lr += 0.125 / (2 * freq)
            return self.make_result(locus, lr, dict_make_result)

        # one of them is homozygous, one common allele
        lr += 0.125 / freq

        return self.make_result(locus, lr, dict_make_result)
