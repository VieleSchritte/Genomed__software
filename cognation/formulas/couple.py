from __future__ import unicode_literals
from .base import Formula, Calculations


class CoupleFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        child_alleles, child_set = alleles[0], sets[0]

        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        freq_dict = self.get_frequencies(locus, child_alleles)
        freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
        c = Calculations()
        lr = 0

        # If there are no intersections between child and couple, return lr = 0 and start counting mutations
        for i in range(0, 2):
            if len(intersections[i]) == 0:
                return self.make_result(locus, lr, dict_make_result)

        # aa an an
        if len(child_set) == 1:
            print('aa an an')
            lr = (c.F(freq1)) ** 2
            return self.make_result(locus, lr, dict_make_result)

        # ab an bn
        print('ab an bn')
        lr = 2 * c.F(freq1) * c.F(freq2) - (2 * freq1 * freq2) ** 2
        return self.make_result(locus, lr, dict_make_result)
