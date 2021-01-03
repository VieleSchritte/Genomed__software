from __future__ import unicode_literals
from .base import Formula, Calculations


class CoupleFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        child_alleles, child_set = alleles[2], sets[2]

        if self.is_gender_specific(locus):
            return self.result_gender_specific(locus, dict_make_result)

        freq_dict = self.get_frequencies(locus, child_alleles)
        freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
        c = Calculations()
        for i in range(1, 3):
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)
        if len(child_set) == 1:  # aa an an
            lr = (c.F(freq1)) ** 2
            return self.make_result(locus, 1 / lr, dict_make_result)

        lr = 2 * c.F(freq1) * c.F(freq2) - (2 * freq1 * freq2) ** 2  # ab an bn
        return self.make_result(locus, 1 / lr, dict_make_result)
