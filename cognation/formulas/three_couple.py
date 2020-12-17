from __future__ import unicode_literals
from .base import Formula, Calculations
from .couple import CoupleFormula
from .two_couple import TwoCoupleFormula


class ThreeCoupleFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 5)
        father_alleles, mother_alleles, children_genotypes = alleles[0], alleles[1], alleles[2:]
        children_sets = sets[2:]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)
        for i in range(1, 7):
            if intersections[i] == 0:
                return self.make_result(locus, 0, dict_make_result)

        c = Calculations()
        raw_values = [locus, '/'.join(father_alleles), '/'.join(mother_alleles)]
        formulas = [CoupleFormula(Formula), TwoCoupleFormula(Formula)]
        lr = c.get_repeatable_lr(raw_values, children_genotypes, formulas)
        if lr:
            return self.make_result(locus, lr, dict_make_result)
        else:
            hetero_counter = c.hetero_counter(children_sets)
            children_alleles = c.get_overall_alleles(children_genotypes)
            freq_dict = self.get_frequencies(locus, children_alleles)
            freq1, freq2, freq3 = c.get_correct_frequency_order_couple(children_sets, freq_dict, children_alleles)
            freq_sum = 0
            for allele in children_alleles:
                freq_sum += freq_dict[allele]
            overall_dict = {
                3: {
                    4: c.multiply_lr_on_children_allele(8, children_alleles, freq_dict),
                    3: c.multiply_lr_on_children_allele(8, children_alleles, freq_dict) * freq_sum
                },
                2: 8 * freq1 ** 2 * freq2 * freq3,
                1: c.multiply_lr_on_children_allele(2, children_alleles, freq_dict) ** 2
            }
            lr = c.get_lr_from_dict_couple(overall_dict, hetero_counter, len(children_alleles))
            return self.make_result(locus, 1 / lr, dict_make_result)
