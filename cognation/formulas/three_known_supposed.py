from __future__ import unicode_literals
from .base import Formula, Calculations
from .one_known_supposed import OneKnownSupposedFormula
from .two_known_supposed import TwoKnownSupposedFormula


class ThreeKnownSupposed(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 5)
        known_alleles, supposed_alleles, children_genotypes = alleles[0], alleles[1], alleles[2:]
        known_set, supposed_set, children_sets = sets[0], sets[1], sets[2:]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        for i in range(1, 7):
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)

        c = Calculations()
        raw_values = [locus, '/'.join(known_alleles), '/'.join(supposed_alleles)]
        formulas = [OneKnownSupposedFormula(Formula), TwoKnownSupposedFormula(Formula)]
        lr = c.get_repeatable_lr(raw_values, children_genotypes, formulas)
        if lr:
            return self.make_result(locus, lr, dict_make_result)
        else:
            common = c.get_overall_alleles(alleles)
            freq_dict = self.get_frequencies(locus, common)
            possible_parent_genotypes = c.get_possible_genotypes(children_sets, children_genotypes, known_set)
            print(possible_parent_genotypes)
