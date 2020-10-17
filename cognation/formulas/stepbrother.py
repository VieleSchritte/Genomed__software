from __future__ import unicode_literals
from .base import Formula


class StepbrotherFormula(Formula):
    def calculate_relation(self, raw_values):
        (stb1_alleles, stb2_alleles, locus, stb1_set, stb2_set, intersection) = self.getting_alleles_locus(raw_values, 2)

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr = '-'
        if self.is_gender_specific(locus):
            return self.make_result2(locus, '/'.join(stb2_alleles), '/'.join(stb1_alleles), '-')

        freq_dict = self.get_frequencies(locus, intersection)

        lr = 0.5
        if len(intersection) == 0:
            return self.make_result2(locus, '/'.join(stb1_alleles), '/'.join(stb2_alleles), lr)

        if len(intersection) == 2:
            freq1 = freq_dict[list(intersection)[0]]
            freq2 = freq_dict[list(intersection)[1]]
            lr += 0.25 * (freq1 + freq2) / (2 * (freq1 * freq2))
            return self.make_result2(locus, '/'.join(stb1_alleles), '/'.join(stb2_alleles), lr)

        freq = freq_dict[list(intersection)[0]]

        if len(stb1_set) == len(stb2_set) == 2:
            lr += 0.25 / (2 * freq)

        if len(stb1_set) == len(stb2_set) == 1:
            lr += 0.5 / freq

        if len(stb1_set) != len(stb2_set):
            lr += 0.25 / freq

        return self.make_result2(locus, '/'.join(stb1_alleles), '/'.join(stb2_alleles), lr)
