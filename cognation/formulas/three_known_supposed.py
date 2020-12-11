from __future__ import unicode_literals
from .base import Formula, Calculations
from .one_known_supposed import OneKnownSupposedFormula
from .two_known_supposed import TwoKnownSupposedFormula


class ThreeKnownSupposed(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 5)
        known_alleles, supposed_alleles, children_genotypes = alleles[0], alleles[1], alleles[2:]
        known_set, child1_set, child2_set, child3_set, supposed_set = sets

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

            # ab ac bc ab ac/bc
            if len(common_set) == 3 and len(child1_set) == 2:
                freq1, freq2 = freq_dict[children_genotypes[0][0]], freq_dict[children_genotypes[0][1]]
                freq3 = freq_dict[list(child2_set - child1_set)[0]]
                lr = 2 * freq3 * (freq1 + freq2)
                return self.make_result(locus, lr, dict_make_result)

            # default describes all other cases
            freq1, freq2 = freq_dict[supposed_alleles[0]], freq_dict[supposed_alleles[1]]
            lr = 2 * freq1 * freq2
            return self.make_result(locus, lr, dict_make_result)
