from __future__ import unicode_literals
from .base import Formula, AllelesException

class StepbrotherFormula(Formula):
    def calculate_relation(self, raw_values):
        (stb1_alleles, stb2_alleles, locus, stb1_set, stb2_set, intersection) = self.getting_alleles_locus(raw_values)

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr
        if self.is_gender_specific(locus):
            return self.make_result(locus, '/'.join(stb2_alleles), '/'.join(stb1_alleles), '-')

        # Skip line with warning if there's incorrect number of alleles
        if len(stb1_alleles) != 2 or len(stb2_alleles) != 2:
            raise AllelesException()

        freq_dict = self.get_frequencies(locus, intersection)
        inter_list, uncle_set_list, nephew_set_list = list(intersection), list(stb2_set), list(stb1_set)

        lr = 0.5
        if len(inter_list) == 0:
            return self.make_result(locus, '/'.join(stb1_alleles), '/'.join(stb2_alleles), lr)

        if len(inter_list) == 2:
            freq1 = freq_dict[inter_list[0]]
            freq2 = freq_dict[inter_list[1]]
            lr += 0.25 * (freq1 + freq2) / (2 * (freq1 * freq2))
            return self.make_result(locus, '/'.join(stb1_alleles), '/'.join(stb2_alleles), lr)

        freq = freq_dict[inter_list[0]]

        if len(uncle_set_list) == len(nephew_set_list) == 2:
            lr += 0.25 / (2 * freq)

        if len(uncle_set_list) == len(nephew_set_list) == 1:
            lr += 0.5 / freq

        if len(uncle_set_list) != len(nephew_set_list):
            lr += 0.25 / freq

        return self.make_result(locus, '/'.join(stb1_alleles), '/'.join(stb2_alleles), lr)