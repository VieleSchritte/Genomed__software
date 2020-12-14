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

        for i in range(len(intersections)):
            exceptions_list = [2, 4, 5]
            if i in exceptions_list:
                continue
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)

        c = Calculations()
        formulas = [OneKnownSupposedFormula(Formula), TwoKnownSupposedFormula(Formula)]
        raw_values = [locus, '/'.join(known_alleles), '/'.join(supposed_alleles)]
        lr = c.get_repeatable_lr(raw_values, children_genotypes, formulas)
        if lr:
            return self.make_result(locus, lr, dict_make_result)
        else:
            common = []
            for genotype in alleles:
                common += genotype
            common_set = set(common)
            freq_dict = self.get_frequencies(locus, common)

            lr = 2
            homo_counter = 0
            for child_set in children_sets:
                if len(child_set) == 1:
                    homo_counter += 1
            if homo_counter == 2:  # aa ab bb ab ab
                for allele in freq_dict.keys():
                    lr *= freq_dict[allele]
                return self.make_result(locus, lr, dict_make_result)
            """homo_dict = {
                1: {
                    (1, 1, 1): 2 * freq1 * freq3,
                    (1, 0, 1): [2 * freq1 * freq3, 2 * freq1 * freq2],
                    (0, 1, 1): [2 * freq1 * freq3, 2 * freq1 * freq2],
                    (1, 1, 0): [2 * freq1 * freq3, 2 * freq1 * freq2],
                }
                0: {
                    (1, 1, 1): 2 * freq3 * (freq1 + freq2),
                    (1, 1, 0): 2 * freq2 * freq4,
                    (1, 0, 1): 2 * freq2 * freq4,
                    (0, 1, 1): 2 * freq2 * freq4,
                }
            }"""
