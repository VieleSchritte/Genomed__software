from __future__ import unicode_literals
from .base import Formula, Calculations
from .brother import BrotherFormula


class TwoBrothersFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        brother2_alleles, brother1_alleles, inspected_alleles = alleles
        brother2_set, brother1_set, inspected_set = sets
        br2br1_inter, br2insp_inter, br1insp_inter = intersections

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr = '-'
        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        c = Calculations()
        common_set = set(inspected_alleles + brother1_alleles + brother2_alleles)
        freq_dict = self.get_frequencies(locus, list(common_set))
        confirmation = 1

        # special cases aa aa any, ab ab any, ab aa bn => confirmation = 1
        if brother1_set == inspected_set:
            lr = self.get_division_lr(locus, inspected_set, inspected_alleles, confirmation)
            return self.make_result(locus, lr, dict_make_result)

        # aa ab cc/cd cases => lr = 0
        zero_case1 = len(br2br1_inter) == len(br2insp_inter) == 0
        # ab not an/bn not an/bn => lr = 0
        zero_case2 = len(br1insp_inter) == len(br2insp_inter) == 0
        # ab ac ad => lr = 0
        zero_case3 = len(common_set) == 4 and br1insp_inter == br2br1_inter == br2insp_inter
        if zero_case1 or zero_case2 or zero_case3:
            return self.make_result(locus, 0, dict_make_result)

        # brother1_set = brother2_set => call BrotherFormula
        if brother1_set == brother2_set:
            raw_values = [locus, '/'.join(inspected_alleles), '/'.join(brother1_alleles)]
            result = BrotherFormula(Formula).calculate_relation(raw_values)
            lr = result['lr']
            return self.make_result(locus, lr, dict_make_result)

        # Homozygous inspected person
        if len(inspected_set) == 1:
            freq1, freq2 = freq_dict[inspected_alleles[0]], freq_dict[list(brother1_set - inspected_set)[0]]
            refutation = c.homo_refutation(freq1)

            # aa ab ac
            if br1insp_inter == br2insp_inter:
                confirmation = 2 * freq1 / (2 - freq1 + 2 * freq2)
                return self.make_result(locus, confirmation / refutation, dict_make_result)

            # aa ab bc
            confirmation = 2 * freq1 / (2 + freq2)
            return self.make_result(locus, confirmation / refutation, dict_make_result)

        # Heterozygous inspected person
        else:
            freq1, freq2 = freq_dict[inspected_alleles[0]], freq_dict[inspected_alleles[1]]
            refutation = c.hetero_refutation(freq1, freq2)

            # ab aa bn
            if len(brother1_set) == 1 and len(br2br1_inter) == 0:
                return self.make_result(locus, confirmation / refutation, dict_make_result)

            # ab aa ac
            if len(brother1_set) == 1 and len(br1insp_inter) != 0:
                freq1, freq2 = freq_dict[brother1_alleles[0]], freq_dict[list(inspected_set - brother1_set)[0]]
                confirmation = freq2 / (2 - freq1)
                return self.make_result(locus, confirmation / refutation, dict_make_result)

            # ab ac bd
            if len(br2br1_inter) == 0 and len(br1insp_inter) != 0:
                return self.make_result(locus, 0.5 / refutation, dict_make_result)

            # ab ac cc
            if len(brother2_set) == 1 and br1insp_inter != br2br1_inter:
                freq1, freq3 = freq_dict[list(br1insp_inter)[0]], freq_dict[brother2_alleles[0]]
                freq2 = freq_dict[list(inspected_set - brother1_set)[0]]
                confirmation = 2 * freq2 * freq3 / c.F(freq1)
                return self.make_result(locus, confirmation / refutation, dict_make_result)

            # ab ac bc
            if len(common_set) == 3:
                freq1, freq2 = freq_dict[list(br1insp_inter)[0]], freq_dict[list(brother1_set - inspected_set)[0]]
                confirmation = 2 * freq1 / (2 + freq2)
                return self.make_result(locus, confirmation / refutation, dict_make_result)

            # ab ac cd
            freq1, freq2 = freq_dict[list(inspected_set - brother1_set)[0]], freq_dict[list(brother1_set - inspected_set)[0]]
            confirmation = 2 * freq1 / (2 + freq2)
            return self.make_result(locus, confirmation / refutation, dict_make_result)
