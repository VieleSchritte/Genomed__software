from __future__ import unicode_literals
from cognation.formulas.base import Formula, AllelesException
from .base import Calculations


class SiblingFormula(Formula):
    def calculate_relation(self, raw_values):

        (locus, alleles_list, sets_list, inter_list) = self.getting_alleles_locus(raw_values, 3)

        child_alleles, parent_alleles, sibling_alleles = alleles_list
        child_set, parent_set, sibling_set = sets_list

        #  cp = child and parent, cs = child and sibling, ps = parent and sibling
        intersection_cp, intersection_cs, intersection_ps = inter_list

        if self.is_gender_specific(locus):
            return self.make_result(locus, '/'.join(child_alleles), '/'.join(parent_alleles), '-', '/'.join(sibling_alleles))

        #  if there's no common alleles between parent and child or between parent and sibling, there's no relation
        if len(intersection_cp) == 0 or len(intersection_ps) == 0:
            return self.make_result(locus, '/'.join(child_alleles), '/'.join(parent_alleles), 0, '/'.join(sibling_alleles))

        if len(child_alleles) != 2 or len(parent_alleles) != 2 or len(sibling_alleles) != 2:
            raise AllelesException()

        lr = 1
        calc = Calculations()

        if len(child_set) == 1:

            freq_dict = self.get_frequencies(locus, child_set)
            refutation = calc.homo_refutation(freq_dict[0])

            freq1, freq2 = self.get_frequencies(locus, sibling_set)

            #  if both child and sibling are homozygous with the common allele then return lr = 1
            if len(sibling_set) == 1 and len(intersection_cs) == 1:
                return self.make_result(locus, '/'.join(child_alleles), '/'.join(parent_alleles), lr, '/'.join(sibling_alleles))

            
