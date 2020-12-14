from __future__ import unicode_literals
from .base import Formula, Calculations
from .one_known_supposed import OneKnownSupposedFormula


class TwoKnownSupposedFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 4)
        known_alleles, supposed_alleles, children_genotypes = alleles[0], alleles[1], alleles[2:]
        known_set, supposed_set, child1_set, child2_set = sets
        ch1ch2_inter = intersections[-1]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        # If there are no intersections between children and parents, return lr = 0 and start counting mutations
        for i in range(1, 5):
            if i == 2 or i == 3:
                continue
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)
        common = []
        for genotype in alleles:
            common += genotype
        common_set, freq_dict = set(common), self.get_frequencies(locus, common)
        c = Calculations()

        # If children's genotypes are same, use OneKnownSupposedFormula
        raw_values = [locus, '/'.join(known_alleles), '/'.join(supposed_alleles)]
        lr = c.get_repeatable_lr(raw_values, children_genotypes, [OneKnownSupposedFormula(Formula)])
        if lr:
            return self.make_result(locus, lr, dict_make_result)


