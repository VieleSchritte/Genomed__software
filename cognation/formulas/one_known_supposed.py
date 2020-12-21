from __future__ import unicode_literals
from .base import Formula, Calculations


class OneKnownSupposedFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        child_alleles = alleles[2]
        known_set, child_set = sets[0], sets[2]
        kch_inter, sch_inter = intersections[1], intersections[2]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        lr = 0
        for i in range(1, 3):
            if len(intersections[i]) == 0:
                return self.make_result(locus, lr, dict_make_result)

        freq_dict = self.get_frequencies(locus, child_alleles)
        freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
        c = Calculations()

        pos_dict = {
            len(child_set) == 1: c.F(freq1),  # aa an an
            child_set == known_set and len(child_set) == 2: (freq1 + freq2) * (2 - (freq1 + freq2)),  # ab ab an/bn
            child_set != known_set and len(child_set) == 2: c.F(freq_dict[list(sch_inter)[0]])  # ab an bn
        }
        for key in pos_dict.keys():
            if key:
                lr = pos_dict[key]
                return self.make_result(locus, 1 / lr, dict_make_result)
