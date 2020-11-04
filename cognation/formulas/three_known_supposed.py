from __future__ import unicode_literals
from .base import Formula
from .one_known_supposed import OneKnownSupposedFormula
from .two_known_supposed import TwoKnownSupposedFormula


class ThreeKnownSupposed(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        supposed_alleles, known_alleles, child3_alleles, child2_alleles, child1_alleles = alleles
        supposed_set, known_set, child3_set, child2_set, child1_set = sets
        sk_inter, sch3_inter, sch2_inter, sch1_inter, kch3_inter, kch2_inter, kch1_inter, ch3ch2_inter, ch3ch1_inter, ch2ch1_inter = intersections

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr = '-'
        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        children_alleles = [child1_alleles, child2_alleles, child3_alleles]

        # If child alleles are the same, use OneKnownSupposed
        if child1_set == child2_set == child3_set:
            raw_values = [locus, '/'.join(child1_alleles), '/'.join(known_alleles), '/'.join(supposed_alleles)]
            result = OneKnownSupposedFormula(Formula).calculate_relation(raw_values)
            result['part4'] = '/'.join(child2_alleles)
            result['part5'] = '/'.join(child3_alleles)
            return result

        unique_genotype = []
        repeat_genotype = []

        for i in range(len(children_alleles)):
            for j in range(len(children_alleles)):
                for k in range(len(children_alleles)):

                    if j > i and children_alleles[i] == children_alleles[j]:
                        if child_alleles[k] != children_alleles[i]:
                            unique_genotype = children_alleles[k]
                            repeat_genotype = children_alleles[i]

        # If there are two same child alleles, use TwoKnownSupposed

        if len(unique_genotype) != 0 and len(repeat_genotype) != 0:
            raw_values = [locus, '/'.join(supposed_alleles), '/'.join(known_alleles), '/'.join(unique_genotype), '/'.join(repeat_genotype)]
            result = TwoKnownSupposedFormula(Formula).calculate_relation(raw_values)
            result['part5'] = '/'.join(child3_alleles)
            return result
