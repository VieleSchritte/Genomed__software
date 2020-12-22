from __future__ import unicode_literals
from .base import Formula


class StepbrotherFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 2)
        insp_set, step_bro = sets
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

        if len(insp_set) == len(step_bro) == 2:
            lr += 0.25 / (2 * freq)

        if len(insp_set) == len(step_bro) == 1:
            lr += 0.5 / freq

        if len(insp_set) != len(step_bro):
            lr += 0.25 / freq

        return self.make_result(locus, lr, dict_make_result)
