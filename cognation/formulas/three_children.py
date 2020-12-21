from __future__ import unicode_literals
from .base import Formula, Calculations
from .two_children import TwoChildrenFormula
from .parent import ParentFormula


class ThreeChildrenFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 4)
        parent_alleles, children_genotypes = alleles[0], alleles[1:]
        parent_set = sets[0]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        for i in range(3):  # no intersections between children and parent => return lr = 0
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)

        raw_values = [locus, '/'.join(parent_alleles)]
        formulas = [ParentFormula(Formula), TwoChildrenFormula(Formula)]
        c = Calculations()
        lr = c.get_repeatable_lr(raw_values, children_genotypes, formulas)
        if lr:
            return self.make_result(locus, lr, dict_make_result)
        else:
            children_alleles = c.get_overall_alleles(children_genotypes)
            freq_dict = self.get_frequencies(locus, children_alleles)
            possible_parent_genotypes = c.get_possible_genotypes(children_alleles, children_genotypes, [parent_set, 'supposed'])
            lr = c.get_lr_from_possible(possible_parent_genotypes, freq_dict)
            return self.make_result(locus, 1 / lr, dict_make_result)
