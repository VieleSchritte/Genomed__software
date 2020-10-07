from __future__ import unicode_literals
from .base import Formula, AllelesException


class CousinFormula(Formula):
    def calculate_relation(self, raw_values):
        (cousin1_alleles, cousin2_alleles, locus, cousin1_set, cousin2_set, intersection) = self.getting_alleles_locus(raw_values)

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr
        if self.is_gender_specific(locus):
            return self.make_result(locus, '/'.join(cousin2_alleles), '/'.join(cousin1_alleles), '-')

        if len(cousin1_alleles) != 2 or len(cousin2_alleles) != 2:
            raise AllelesException()

        freq_dict = self.get_frequencies(locus, intersection)
        lr = 0.75

        inter_len = len(list(intersection))
        cousin1_set_len = len(list(cousin1_set))
        cousin2_set_len = len(list(cousin2_set))

        # no common alleles
        if inter_len == 0:
            return self.make_result(locus, '/'.join(cousin2_alleles), '/'.join(cousin1_alleles), lr)

        # both are heterozygous, two common alleles
        if inter_len == 2:
            freq1 = freq_dict[list(intersection)[0]]
            freq2 = freq_dict[list(intersection)[1]]
            lr += 0.125 * (freq1 + freq2) / (2 * freq1 * freq2)
            return self.make_result(locus, '/'.join(cousin2_alleles), '/'.join(cousin1_alleles), lr)

        freq = freq_dict[list(intersection)[0]]

        #  both are homozygous, one common allele
        if cousin1_set_len == cousin2_set_len == 1:
            lr += 0.25 / freq
            return self.make_result(locus, '/'.join(cousin2_alleles), '/'.join(cousin1_alleles), lr)

        #  both are heterozygous, one common allele
        if cousin1_set_len == cousin2_set_len == 2:
            lr += 0.125 / (2 * freq)
            return self.make_result(locus, '/'.join(cousin2_alleles), '/'.join(cousin1_alleles), lr)

        # one of them is homozygous, one common allele
        lr += 0.125 / freq

        return self.make_result(locus, '/'.join(cousin2_alleles), '/'.join(cousin1_alleles), lr)
