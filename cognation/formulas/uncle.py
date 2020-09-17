from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException

class UncleFormula(Formula):
    def calculate_relation(self, raw_values):
        lr = 0
        if len(raw_values) < 3:
            raise LineFormatException()

        nephew_alleles = self.split_sat(raw_values.pop())
        uncle_alleles = self.split_sat(raw_values.pop())
        locus = ''.join(raw_values)

        if len(nephew_alleles) != 2 or len(uncle_alleles) != 2:
            raise AllelesException()

        uncle_set = set(uncle_alleles)  # unique uncle alleles
        nephew_set = set(nephew_alleles)  # unique nephew alleles
        intersection_list = list(uncle_set & nephew_set)
        freq_dict = self.get_frequencies(locus, intersection_list)

        # case -1 - no common alleles
        if len(intersection_list) == 0:
            lr = 0.5

        # case 0 - two common alleles, both are heterozigous
        elif len(intersection_list) == 2:
            freq1 = freq_dict[intersection_list[0]]
            freq2 = freq_dict[intersection_list[1]]
            lr = 0.5 + (freq1 + freq2)/(8 * freq1 * freq2)

        # one common allele
        elif len(intersection_list) ==1:
            freq = freq_dict[intersection_list[0]]

            # case 1 - both are homozigous
            if len(uncle_set) == len(nephew_set) == 1:
                lr = 0.5 + 1/(2 * freq)

            # case 2 - uncle is homozigous, nephew is heterozigous
            elif len(uncle_set) == 1 and len(nephew_set) == 2:
                lr = 0.5 + 1/(4 * freq)

            # case 3 - both are heterozigous
            elif len(uncle_set) == len(nephew_set) == 2:
                lr = 0.5 + 1/(8 * freq)

        return self.make_result(locus, '/'.join(uncle_alleles), '/'.join(nephew_alleles), lr)

