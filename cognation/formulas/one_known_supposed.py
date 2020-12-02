from __future__ import unicode_literals
from .base import Formula, Calculations


class OneKnownSupposedFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        known_alleles, child_alleles, sup_alleles = alleles
        known_set, child_set, sup_set = sets
        kch_inter, sk_inter, sch_inter = intersections

        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        if locus == 'AMEL':
            return self.make_result(locus, 1, dict_make_result)

        lr = 0
        if len(sch_inter) == 0 or len(kch_inter) == 0:
            return self.make_result(locus, lr, dict_make_result)

        freq_dict = self.get_frequencies(locus, child_alleles + known_alleles + sup_alleles)
        print(locus)
        print('child, known, sup: ', child_alleles, known_alleles, sup_alleles)
        print(freq_dict)
        freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]

        # cases ab ab an/bn
        if child_set == known_set and len(child_set) == 2:
            lr = (freq1 + freq2) * (2 - (freq1 + freq2))
            print(lr)
            print()
            return self.make_result(locus, lr, dict_make_result)

        c = Calculations()
        # case ab an bn
        if len(child_set) == 2:
            freq = freq_dict[list(child_set - known_set)[0]]
            lr = c.F(freq)
            print(lr)
            print()
            return self.make_result(locus, lr, dict_make_result)

        freq = freq_dict[list(sch_inter)[0]]
        lr = c.F(freq)
        print(lr)
        print()
        return self.make_result(locus, lr, dict_make_result)
