from __future__ import unicode_literals
from .base import Formula, Calculations


class CoupleFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        child_alleles, mother_alleles, father_alleles = alleles
        child_set = sets[0]
        cm_inter, cf_inter = intersections[0:2]

        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        freq_dict = self.get_frequencies(locus, child_alleles + mother_alleles + father_alleles)
        c = Calculations()
        lr = 0

        if len(cm_inter) != 0 and len(cf_inter) != 0:
            if len(child_set) == 1:
                freq = freq_dict[child_alleles[0]]
                lr = (c.F(freq)) ** 2

            else:
                freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
                lr = 2 * c.F(freq1) * c.F(freq2) - (2 * freq1 * freq2) ** 2

        return self.make_result(locus, lr, dict_make_result)
