from __future__ import unicode_literals
from cognation.formulas.base import Formula, Calculations


class YesGrandParent(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        child_alleles, grandparent_alleles = alleles[0], alleles[1]
        child_set, grandparent_set = sets[0], sets[1]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)
        c = Calculations()
        freq_dict = self.get_frequencies(locus, c.get_overall_alleles(alleles))

        if len(child_set) == 1:  # Homozygous child
            freq1, freq2, freq3 = freq_dict[child_alleles[0]], freq_dict[grandparent_alleles[0]], freq_dict[grandparent_alleles[1]]
            answers = {
                len(intersections[0]) != 0: c.F(freq1),  # aa an an
                len(intersections[0]) == 0 and len(grandparent_set) == 1: 2 * freq1 * freq2,  # aa bb ab
                len(intersections[0]) == 0 and len(grandparent_set) == 2: 2 * freq1 * (freq2 + freq3)  # aa bc ab/ac
            }
            for condition in answers.keys():
                if condition:
                    return self.make_result(locus, 1 / answers[condition], dict_make_result)
            return self.make_result(locus, 0, dict_make_result)

        else:  # Heterozygous child
            freq_list = [1, 1, 1, 1]
            cg_inter = intersections[0]
            if len(cg_inter) != 0:
                if len(grandparent_set) == 1:
                    freq_list[0], freq_list[1] = freq_dict[grandparent_alleles[0]], freq_dict[list(child_set - grandparent_set)[0]]
                elif grandparent_set == child_set:
                    freq_list[0], freq_list[1] = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
                else:
                    freq_list[0], freq_list[1] = freq_dict[list(cg_inter)[0]], freq_dict[list(child_set - cg_inter)[0]]
                    freq_list[2] = freq_dict[list(grandparent_set - cg_inter)[0]]
            else:
                freq_list[0], freq_list[1] = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
                freq_list[2], freq_list[3] = freq_dict[grandparent_alleles[0]], freq_dict[grandparent_alleles[1]]
            freq1, freq2, freq3, freq4 = freq_list
            answers = {
                len(cg_inter) != 0 and len(grandparent_set) == 1: c.F(freq1),  # ab aa an
                len(cg_inter) != 0 and grandparent_set == child_set: (freq1 + freq2) * (2 - (freq1 + freq2)),  # ab ab an/bn
                len(cg_inter) != 0 and len(grandparent_set) == 2 and grandparent_set != child_set: c.F(freq1) + 2 * freq2 * freq3,  # ab ac an/bc
                len(cg_inter) == 0 and len(grandparent_set) == 1: 2 * freq3 * (freq1 + freq2),  # ab cc ac/bc
                len(cg_inter) == 0 and len(grandparent_set) != 1: 2 * (freq1 + freq2) * (freq3 + freq4)  # ab cd ac/ad/bc/bd
            }
            for key in answers.keys():
                if key:
                    return self.make_result(locus, 1 / answers[key], dict_make_result)
            return self.make_result(locus, 0, dict_make_result)
