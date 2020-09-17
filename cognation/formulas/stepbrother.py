from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException

# !!! The same formula for uncle !!!
class StepbrotherFormula(Formula):
    def calculate_relation(self, raw_values):
        lr = 0
        if len(raw_values) < 3:
            raise LineFormatException()

        stepbrother1_alleles = self.split_sat(raw_values.pop())
        stepbrother2_alleles = self.split_sat(raw_values.pop())
        locus = ''.join(raw_values)

        if len(stepbrother1_alleles) != 2 or len(stepbrother2_alleles) != 2:
            raise AllelesException()

        stepbrother2_set = set(stepbrother2_alleles)  # unique uncle alleles
        stepbrother1_set = set(stepbrother1_alleles)  # unique nephew alleles
        intersection_list = list(stepbrother2_set & stepbrother1_set)
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
            if len(stepbrother2_set) == len(stepbrother1_set) == 1:
                lr = 0.5 + 1/(2 * freq)

            # case 2 - uncle is homozigous, nephew is heterozigous
            elif len(stepbrother2_set) == 1 and len(stepbrother1_set) == 2:
                lr = 0.5 + 1/(4 * freq)

            # case 3 - both are heterozigous
            elif len(stepbrother2_set) == len(stepbrother1_set) == 2:
                lr = 0.5 + 1/(8 * freq)

        return self.make_result(locus, '/'.join(stepbrother2_alleles), '/'.join(stepbrother1_alleles), lr)


