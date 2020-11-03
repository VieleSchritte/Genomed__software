from __future__ import unicode_literals
from .base import Formula
from .base import Calculations


class TwoParentsFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        child_alleles, mother_alleles, father_alleles = alleles
        child_set, parent1_set, parent2_set = sets
        intersection1, intersection2 = intersections[0], intersections[1]

        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        freq_dict = self.get_frequencies(locus, child_alleles + mother_alleles + father_alleles)
        c = Calculations()
        lr = 0

        if len(intersection1) != 0 and len(intersection2) != 0:
            if len(child_set) == 1:
                print('homo')
                print()
                freq = freq_dict[child_alleles[0]]
                lr = (c.F(freq)) ** 2

            else:
                print('hetero')
                print()
                freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
                lr = 2 * c.F(freq1) * c.F(freq2) - (2 * freq1 * freq2) ** 2

        return self.make_result(locus, lr, dict_make_result)
