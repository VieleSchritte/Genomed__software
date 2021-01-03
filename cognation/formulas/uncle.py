from __future__ import unicode_literals
from .base import Formula, Calculations


class UncleFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 2)
        nephew_set, uncle_set = sets
        intersection = intersections[0]
        c = Calculations()

        if self.is_gender_specific(locus):
            return self.result_gender_specific(locus, dict_make_result)
        freq_dict = self.get_frequencies(locus, intersection)
        if len(intersection) == 0:
            return self.make_result(locus, 0.5, dict_make_result)
        if len(intersection) == 2:
            freq1, freq2 = freq_dict[list(intersection)[0]], freq_dict[list(intersection)[1]]
            lr = 0.5 + 0.25 * (freq1 + freq2) / (2 * (freq1 * freq2))
            return self.make_result(locus, lr, dict_make_result)
        freq = freq_dict[list(intersection)[0]]
        answers = {
            len(uncle_set) == len(nephew_set) == 2: 0.5 + 0.25 / (2 * freq),
            len(uncle_set) == len(nephew_set) == 1: 0.5 + 0.5 / freq,
            len(uncle_set) != len(nephew_set): 0.5 + 0.25 / freq
        }
        lr = c.get_lr_from_cond_dict_short(answers)
        return self.make_result(locus, lr, dict_make_result)
