from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException

class ParentFormula(Formula):
    def calculate_relation(self, raw_values):
        if len(raw_values) < 3:
            raise LineFormatException()

        alleles_locus_tuple = self.getting_alleles_locus(raw_values)

        child_alleles = alleles_locus_tuple[0]
        parent_alleles = alleles_locus_tuple[1]
        locus = alleles_locus_tuple[2]
        child_set = alleles_locus_tuple[3]
        parent_set = alleles_locus_tuple[4]
        intersection = alleles_locus_tuple[5]

        # Function in base.py for checking out if the locus is gender-specific; if yes it returnes lr
        gender_check = self.gender_specific(locus)

        if gender_check != None:
            relation = gender_check
        else:

            if len(parent_alleles) != 2 or len(child_alleles) != 2:
                raise AllelesException()

            relation = 0
            freq_dict = self.get_frequencies(locus, intersection)
            # i, j     i, j
            if len(freq_dict) == 2:
                divider = 4
                dividend = 0
                for key in freq_dict:
                    divider = divider * freq_dict[key]
                    dividend = dividend + freq_dict[key]

                if divider != 0:
                    # (f1 + f2) / 4 * f1 * f2
                    relation = dividend / divider
            elif len(freq_dict) == 1:
                # only 1 allele in common
                # 1 / {1,2}*{1,2}*f
                freq = list(freq_dict.values())[0]
                relation = 1 / (len(parent_set)*len(child_set) * freq)

        return self.make_result(locus, '/'.join(parent_alleles), '/'.join(child_alleles), relation)
