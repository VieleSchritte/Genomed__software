from __future__ import unicode_literals
from .base import Formula, Calculations


class TwoBrothersFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        brother2_alleles, brother1_alleles, inspected_alleles = alleles
        brother2_set, brother1_set, inspected_set = sets
        br2br1_inter, br2insp_inter, br1insp_inter = intersections

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr = '-'
        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        lr = 0
        if br2br1_inter == 0 or br2insp_inter == 0:
            return self.make_result(locus, lr, dict_make_result)

        c = Calculations()
        freq_dict = self.get_frequencies(locus, inspected_alleles + brother1_alleles + brother2_alleles)
        confirmation = 1
        common_set = set(inspected_alleles + brother1_alleles + brother2_alleles)

        # aa ab bc, ab ac bc, ab ac cd
        if len(br1insp_inter) == len(br2br1_inter) == 1 and br1insp_inter != br2br1_inter and len(brother2_set) == 2:
            freq1, freq2 = freq_dict[list(br1insp_inter)[0]], freq_dict[list(br2br1_inter)[0]]
            confirmation = freq1 / (2 + freq2)

        # Homozygous inspected person
        if len(inspected_set) == 1:
            freq = freq_dict[inspected_alleles[0]]
            refutation = c.homo_refutation(freq)

            # aa ab ac
            if len(brother1_set) == len(brother2_set) == 2:
                freq1, freq2 = freq, freq_dict[list(brother1_set - inspected_set)[0]]
                confirmation = 2 * freq1 / (2 - freq1 + 2 * freq2)
                lr = confirmation / refutation
                return self.make_result(locus, lr, dict_make_result)

            # aa aa any => confirmation = 1
            lr = confirmation / refutation
            return self.make_result(locus, lr, dict_make_result)

        # Heterozygous inspected person
        else:
            freq1, freq2 = freq_dict[inspected_alleles[0]], freq_dict[inspected_alleles[1]]
            refutation = c.hetero_refutation(freq1, freq2)

            # ab aa ac
            if len(brother1_set) == 1 and len(brother2_set) == 2 and len(common_set) == 2:
                freq1, freq2 = freq_dict[brother1_alleles[0]], freq_dict[list(inspected_set - brother1_set)[0]]
                confirmation = freq2 / (2 - freq1)
                lr = confirmation / refutation
                return self.make_result(locus, lr, dict_make_result)

            # ab ac cc
            if len(brother1_set) != len(brother2_set) and br2insp_inter != br1insp_inter:
                freq1, freq2, freq3 = freq_dict[list(br1insp_inter)[0]], freq_dict[list(inspected_set - brother1_set)[0]], freq_dict[brother2_alleles[0]]
                confirmation = 2 * freq2 * freq3 / c.F(freq1)
                lr = confirmation / refutation
                return self.make_result(locus, lr, dict_make_result)

            # ab ac bd
            if len(br2br1_inter) == 0 and len(common_set) == 4:
                confirmation = 0.5
                lr = confirmation / refutation
                return self.make_result(locus, lr, dict_make_result)

            # ab aa bn, ab ab any
            lr = confirmation / refutation
            return self.make_result(locus, lr, dict_make_result)
