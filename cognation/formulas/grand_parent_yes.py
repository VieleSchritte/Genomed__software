from __future__ import unicode_literals
from cognation.formulas.base import Formula, Calculations


class YesGrandParent(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        child_alleles, grandparent_alleles = alleles[0], alleles[1]
        child_set, grandparent_set = sets[0], sets[1]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)
        for i in range(1, 3):
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)
        c = Calculations()
        freq_dict = self.get_frequencies(locus, c.get_overall_alleles(alleles))

        if len(child_set) == 1:  # Homozygous child
            freq1, freq2, freq3 = freq_dict[child_alleles[0]], freq_dict[grandparent_alleles[0]], freq_dict[grandparent_alleles[1]]
            answers = {
                len(intersections[0]) != 0: c.F(freq1),  # aa an an
                len(intersections[0]) == 0 and len(grandparent_set) == 1: 2 * freq1 * freq2,  # aa bb ab
                len(intersections[0]) == 0 and len(grandparent_set) == 2: 2 * freq1 * (freq2 + freq3)  # aa bc ab/ac
            }
            lr = c.get_lr_from_cond_dict_short(answers)
            return self.make_result(locus, 1 / lr, dict_make_result)

        else:  # Heterozygous child
            cg_inter = intersections[0]
            freq1, freq2, freq3, freq4 = self.get_freq_order(cg_inter, grandparent_set, freq_dict, alleles, child_set)
            answers = {
                len(cg_inter) != 0 and len(grandparent_set) == 1: c.F(freq1),  # ab aa an
                len(cg_inter) != 0 and grandparent_set == child_set: (freq1 + freq2) * (2 - (freq1 + freq2)),  # ab ab an/bn
                len(cg_inter) != 0 and len(grandparent_set) == 2 and grandparent_set != child_set: c.F(freq1) + 2 * freq2 * freq3,  # ab ac an/bc
                len(cg_inter) == 0 and len(grandparent_set) == 1: 2 * freq3 * (freq1 + freq2),  # ab cc ac/bc
                len(cg_inter) == 0 and len(grandparent_set) != 1: 2 * (freq1 + freq2) * (freq3 + freq4)  # ab cd ac/ad/bc/bd
            }
            lr = c.get_lr_from_cond_dict_short(answers)
            return self.make_result(locus, 1 / lr, dict_make_result)

    @staticmethod
    def get_freq_order(cg_inter, grandparent_set, freq_dict, alleles, child_set):
        child_alleles, grandparent_alleles = alleles[0], alleles[1]
        if len(cg_inter) != 0:
            if len(grandparent_set) == 1:
                freq1, freq2 = freq_dict[grandparent_alleles[0]], freq_dict[list(child_set - grandparent_set)[0]]
                return freq1, freq2, 1, 1
            elif grandparent_set == child_set:
                freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
                return freq1, freq2, 1, 1
            else:
                freq1, freq2 = freq_dict[list(cg_inter)[0]], freq_dict[list(child_set - cg_inter)[0]]
                freq3 = freq_dict[list(grandparent_set - cg_inter)[0]]
                return freq1, freq2, freq3, 1
        else:
            freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
            freq3, freq4 = freq_dict[grandparent_alleles[0]], freq_dict[grandparent_alleles[1]]
            return freq1, freq2, freq3, freq4
