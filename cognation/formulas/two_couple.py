from __future__ import unicode_literals
from .base import Formula
from .couple import CoupleFormula
from .base import Calculations


class TwoCoupleFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 4)
        father_alleles, mother_alleles, children_genotypes = alleles[0], alleles[1], alleles[2:]
        children_sets = sets[2:]

        if self.is_gender_specific(locus):
            return self.result_gender_specific(locus, dict_make_result)
        for i in range(1, 5):
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)

        raw_values = [locus, '/'.join(father_alleles), '/'.join(mother_alleles)]
        c = Calculations()
        lr = c.get_repeatable_lr(raw_values, children_genotypes, [CoupleFormula(Formula)])
        if lr:
            return self.make_result(locus, lr, dict_make_result)
        else:
            children_alleles = c.get_overall_alleles(children_genotypes)
            freq_dict = self.get_frequencies(locus, children_alleles)
            hetero_counter = c.hetero_counter(children_sets)
            freq1, freq2, freq3 = 0, 0, 0
            if len(children_alleles) < 4:
                freq1, freq2, freq3 = c.get_correct_frequency_order_couple(children_sets, freq_dict, children_alleles)
            overall_dict = {
                2: {
                    3: c.multiply_lr_on_children_allele(4, children_alleles, freq_dict) * (2 + freq1),  # ab ac an bc
                    4: c.multiply_lr_on_children_allele(16, children_alleles, freq_dict)  # ab cd ad/ac bc/bd
                },
                1: {
                    2: 4 * freq2 * (2 - freq1 - freq2) * freq1 ** 2,  # aa ab an ab
                    3: 8 * freq2 * freq3 * freq1 ** 2
                },
                0: (c.multiply_lr_on_children_allele(2, children_alleles, freq_dict)) ** 2  # aa bb ab ab
            }
            lr = c.get_lr_from_dict_couple(overall_dict, hetero_counter, len(children_alleles))
            return self.make_result(locus, 1 / lr, dict_make_result)
