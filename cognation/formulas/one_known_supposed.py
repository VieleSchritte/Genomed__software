from __future__ import unicode_literals
from .base import Formula, Calculations


class OneKnownSupposedFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        child_alleles = alleles[2]
        known_set, child_set = sets[0], sets[2]
        kch_inter, sch_inter = intersections[1], intersections[2]

        if self.is_gender_specific(locus):
            return self.result_gender_specific(locus, dict_make_result)
        for i in range(1, 3):
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)
        freq_dict = self.get_frequencies(locus, child_alleles)
        freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
        c = Calculations()

        answers = {
            len(child_set) == 1: c.F(freq1),  # aa an an
            child_set == known_set and len(child_set) == 2: (freq1 + freq2) * (2 - (freq1 + freq2)),  # ab ab an/bn
            child_set != known_set and len(child_set) == 2: c.F(freq_dict[list(sch_inter)[0]])  # ab an bn
        }
        lr = c.get_lr_from_cond_dict_short(answers)
        return self.make_result(locus, 1 / lr, dict_make_result)
