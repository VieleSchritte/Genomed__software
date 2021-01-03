from __future__ import unicode_literals
from .base import Formula
from .base import Calculations


class CousinFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 2)
        inspected_genotype, intersection, inspected_set = alleles[0], intersections[0], sets[0]

        if self.is_gender_specific(locus):
            return self.result_gender_specific(locus, dict_make_result)
        freq_dict = self.get_frequencies(locus, inspected_genotype)
        c = Calculations()
        alleles_number = len(c.get_overall_alleles(alleles))
        freq1, freq2 = freq_dict[inspected_genotype[0]], freq_dict[inspected_genotype[1]]
        if len(inspected_set) == 2 and len(intersection) == 1:
            freq1, freq2 = freq_dict[list(intersection)[0]], freq_dict[list(inspected_set - intersection)[0]]
        answer_dict = {
            (0, 2): 0.75,
            (0, 3): 0.75,
            (0, 4): 0.75,
            (1, 1): 0.75 + 0.25 / freq1,
            (1, 2): 0.75 + 0.125 / freq1,
            (1, 3): 0.75 + 0.125 / (2 * freq1),
            (2, 2): 0.75 + 0.125 * (freq1 + freq2) / (2 * (freq1 * freq2))
        }
        if (len(intersection), alleles_number) in answer_dict.keys():
            lr = answer_dict[(len(intersection), alleles_number)]
            return self.make_result(locus, lr, dict_make_result)
