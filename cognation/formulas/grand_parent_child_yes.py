from __future__ import unicode_literals
from cognation.formulas.base import Formula, Calculations


class ParentGrandChildYes(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        parent_alleles, grandparent_alleles, child_alleles = alleles
        parent_set, grandparent_set, child_set = sets
        pg_inter, pch_inter, gch_inter = intersections

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr = '-'
        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        c = Calculations()
        freq_dict = self.get_frequencies(locus, parent_alleles + grandparent_alleles + child_alleles)
        common_set = set(parent_alleles + grandparent_alleles + child_alleles)

        # Homozygous child
        if len(child_set) == 1:
            # aa an an
            if gch_inter == pch_inter:
                freq = freq_dict[child_alleles[0]]
                lr = c.F(freq)
                return self.make_result(locus, lr, dict_make_result)

            # aa bb ab
            if len(grandparent_set) == 1 and pg_inter != gch_inter:
                freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[grandparent_alleles[0]]
                lr = 2 * freq1 * freq2
                return self.make_result(locus, lr, dict_make_result)

            # aa bc ab/ac
            if len(common_set) == 3 and len(gch_inter) == 0:
                freq1, freq2, freq3 = freq_dict[child_alleles[0]], freq_dict[grandparent_alleles[0]], freq_dict[grandparent_alleles[1]]
                lr = 2 * freq1 * (freq2 + freq3)
                return self.make_result(locus, lr, dict_make_result)

        # Heterozygous child
        else:
            # ab aa an
            if gch_inter == pg_inter and len(grandparent_set) == 1:
                freq = freq_dict[grandparent_alleles[0]]
                lr = c.F(freq)
                return self.make_result(locus, lr, dict_make_result)

            # ab ab an/bn
            if child_set == grandparent_set:
                freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
                lr = (freq1 + freq2) * (2 - (freq1 + freq2))
                return self.make_result(locus, lr, dict_make_result)

            # ab ac an/bc
            if len(child_set) == len(grandparent_set) and child_set != grandparent_set:
                freq1 = freq_dict[list(gch_inter)[0]]
                freq2, freq3 = freq_dict[list(child_set - grandparent_set)[0]], freq_dict[list(grandparent_set - child_set)[0]]
                lr = c.F(freq1) + 2 * freq2 * freq3
                return self.make_result(locus, lr, dict_make_result)

            # ab cc ac/bc
            if len(grandparent_set) == 1 and len(gch_inter) == 0:
                freq1, freq2, freq3 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]], freq_dict[grandparent_alleles[0]]
                lr = 2 * freq3 * (freq1 + freq2)
                return self.make_result(locus, lr, dict_make_result)

            # ab cd ac/ad/bc/bd
            if len(common_set) == 4 and len(gch_inter) == 0:
                freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
                freq3, freq4 = freq_dict[grandparent_alleles[0]], freq_dict[grandparent_alleles[1]]
                lr = 2 * (freq1 + freq2) * (freq3 + freq4)
                return self.make_result(locus, lr, dict_make_result)