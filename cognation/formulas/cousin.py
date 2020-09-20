from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException

class CousinFormula(Formula):
    def calculate_relation(self, raw_values):
        lr = 0
        if len(raw_values) < 3:
            raise LineFormatException()

        cousin2_alleles = self.split_sat(raw_values.pop())
        cousin1_alleles = self.split_sat(raw_values.pop())
        locus = ''.join(raw_values)

        if len(cousin1_alleles) != 2 or len(cousin2_alleles) != 2:
            raise AllelesException()

        cousin1_set = set(cousin1_alleles)  # unique cousin 1 alleles
        cousin2_set = set(cousin2_alleles)  # unique cousin 2 alleles
        intersection = list(cousin1_set & cousin2_set)  # unique common alleles

        cousin1_set_len = len(cousin1_set)
        cousin2_set_len = len(cousin2_set)
        inter_len = len(intersection)

        freq_dict = self.get_frequencies(locus, intersection)

        if inter_len == 0:
            # no common alleles
            lr = 0.75
        elif inter_len == 2:
            # borh are heterozygous, two common alleles
            freq1 = freq_dict[intersection[0]]
            freq2 = freq_dict[intersection[1]]
            lr = 0.75 + 0.125 * (freq1 + freq2) / (2 * freq1 * freq2)
        elif inter_len == 1:
            freq = freq_dict[intersection[0]]
            if cousin1_set_len == cousin2_set_len == 1:
                lr = 0.75 + 0.25 / freq
            elif cousin1_set_len == 1 and cousin2_set_len == 2:
                lr = 0.75 + 0.125 / freq
            elif cousin1_set_len == cousin2_set_len == 2:
                lr = 0.75 + 0.125 / (2 * freq)

        return self.make_result (locus, '/'.join(cousin1_alleles), '/'.join(cousin2_alleles), lr)
