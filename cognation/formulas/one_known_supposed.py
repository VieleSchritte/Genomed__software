from __future__ import unicode_literals
from .base import Formula, Calculations


class OneKnownSupposedFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        sup_alleles, known_alleles, child_alleles = alleles
        sup_set, known_set, child_set = sets
        sk_inter, sch_inter, kch_inter = intersections

        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        lr = 0
        if len(sch_inter) == 0 or len(kch_inter) == 0:
            return self.make_result(locus, lr, dict_make_result)

        freq_dict = self.get_frequencies(locus, child_alleles + known_alleles + sup_alleles)
        # cases ab ab an/bn
        if len(child_set) == len(known_set) == 2 and len(sch_inter) == 1:
            freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
            lr = (freq1 + freq2) * (2 - (freq1 + freq2))
            return self.make_result(locus, lr, dict_make_result)

        freq = freq_dict[list(sch_inter)[0]]
        c = Calculations()
        lr = c.F(freq)
        return self.make_result(locus, lr, dict_make_result)
