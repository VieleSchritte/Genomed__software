from __future__ import unicode_literals
from .base import Formula, Calculations


# FORMULA_TYPE_GRANDPARENT
class GrandParentFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 2)
        grandchild_alleles, grandparent_alleles = alleles
        grandchild_set, grandparent_set = sets
        intersection = intersections[0]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)




class Confirmations:
    @staticmethod
    def homo_gc_confirmation():
        confirmation = 1
        return confirmation

    @staticmethod
    def hetero_gc_confirmation():
        confirmation = 1
        return confirmation