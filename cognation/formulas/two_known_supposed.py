from __future__ import unicode_literals
from .base import Formula, Calculations
from .one_known_supposed import OneKnownSupposedFormula


class TwoKnownSupposedFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 4)
        known_alleles, supposed_alleles, children_genotypes = alleles[0], alleles[1], alleles[2:]
        known_set, children_sets = sets[0], sets[2:]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        # If there are no intersections between children and parents, return lr = 0 and start counting mutations
        for i in range(1, 3):
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)

        # If children's genotypes are same, use OneKnownSupposedFormula
        raw_values = [locus, '/'.join(known_alleles), '/'.join(supposed_alleles)]
        c = Calculations()
        lr = c.get_repeatable_lr(raw_values, children_genotypes, [OneKnownSupposedFormula(Formula)])
        if lr:
            return self.make_result(locus, lr, dict_make_result)
        else:
            children_alleles = c.get_overall_alleles(children_genotypes)
            freq_dict = self.get_frequencies(locus, children_alleles)
            possible_parents_genotypes = c.get_possible_genotypes(children_alleles, children_genotypes, [known_set, 'known'])
            lr = c.get_lr_from_possible(possible_parents_genotypes, freq_dict)
            return self.make_result(locus, 1 / lr, dict_make_result)
