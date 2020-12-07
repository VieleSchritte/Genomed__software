from __future__ import unicode_literals
from .base import Formula


class UncleFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 2)
        nephew_set, uncle_set = sets
        intersection = intersections[0]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        freq_dict = self.get_frequencies(locus, intersection)

        lr = 0.5
        if len(intersection) == 0:
            return self.make_result(locus, lr, dict_make_result)

        if len(intersection) == 2:
            freq1 = freq_dict[list(intersection)[0]]
            freq2 = freq_dict[list(intersection)[1]]
            lr += 0.25 * (freq1 + freq2) / (2 * (freq1 * freq2))
            return self.make_result(locus, lr, dict_make_result)

        freq = freq_dict[list(intersection)[0]]

        if len(uncle_set) == len(nephew_set) == 2:
            lr += 0.25 / (2 * freq)

        if len(uncle_set) == len(nephew_set) == 1:
            lr += 0.5 / freq

        if len(uncle_set) != len(nephew_set):
            lr += 0.25 / freq

        return self.make_result(locus, lr, dict_make_result)
