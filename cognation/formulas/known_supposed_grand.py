from __future__ import unicode_literals
from .base import Formula, Calculations


class KnownSupposedGrand(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 4)
        child_alleles, grandparent_alleles = alleles[0], alleles[2]
        child_set, known_set, grandparent_set = sets[0:3]
        target_sets, target_alleles, target_inters = sets[0:3], [alleles[0], alleles[2]], intersections[0:2]
        ch_known_inter, ch_grand_inter = intersections[0], intersections[1]

        if len(intersections[0]) == 0 or len(intersections[2]) == 0:
            return self.make_result(locus, 0, dict_make_result)
        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)
        c = Calculations()
        all_alleles = c.get_overall_alleles(target_alleles)
        freq_dict = self.get_frequencies(locus, all_alleles)
        freq1, freq2, freq3, freq4 = self.get_freqs_order(target_sets, target_alleles, freq_dict, target_inters)

        if len(child_set) == 1:  # Homozygous child
            answers = {
                len(ch_grand_inter) == 0 and len(grandparent_set) == 1: 2 * freq1 * freq2,  # aa an bb ab
                len(ch_grand_inter) == 0 and len(grandparent_set) == 2: 2 * freq1 * (freq2 + freq3),  # aa an bc b/ac
                len(ch_grand_inter) != 0: c.F(freq1)  # aa an an an
            }
            lr = c.get_lr_from_cond_dict(answers)
            return self.make_result(locus, 1 / lr, dict_make_result)

        if known_set == child_set:
            answers = {
                len(ch_grand_inter) == 0 and len(grandparent_set) == 1: 2 * freq3 * (freq1 + freq2),  # ab ab cc ac/bc
                len(ch_grand_inter) == 0 and len(grandparent_set) == 2: 2 * (freq1 + freq2) * (freq3 + freq4),  # ab ab cd ac/ad/bc/bd
                len(ch_grand_inter) == 2: (freq1 + freq2) * (2 - (freq1 + freq2)),  # ab ab ab an/bn
                len(ch_grand_inter) == 1 and len(grandparent_set) == 2: c.F(freq1) + 2 * freq2 * freq3,  # ab ab ac an/bc
                len(ch_grand_inter) == 1 and len(grandparent_set) == 1: c.F(freq1)  # ab ab aa an
            }
            lr = c.get_lr_from_cond_dict(answers)
            return self.make_result(locus, 1 / lr, dict_make_result)

        child_b_allele = list(child_set - ch_known_inter)
        answers = {
            len(ch_grand_inter) == 0 and len(grandparent_set) == 1: 2 * freq2 * freq3,
            len(ch_grand_inter) == 0 and len(grandparent_set) == 2: 2 * freq2 * (freq3 + freq4),
            len(ch_grand_inter) != 0 and child_b_allele in grandparent_alleles: c.F(child_b_allele),
            len(ch_grand_inter) != 0 and child_b_allele not in grandparent_alleles and len(grandparent_set) == 1: 2 * freq1 * freq2,
            len(ch_grand_inter) != 0 and child_b_allele not in grandparent_alleles and len(grandparent_set) == 2: 2 * freq2 * (freq1 + freq3)
        }
        lr = c.get_lr_from_cond_dict(answers)
        return self.make_result(locus, 1 / lr, dict_make_result)

    @staticmethod
    def get_freqs_order(target_sets, target_alleles, freq_dict, target_inters):
        child_set, known_set, grandparent_set = target_sets
        child_alleles, grandparent_alleles = target_alleles
        ch_known_inter, ch_grand_inter = target_inters
        if len(child_set) == 1:  # Homozygous child
            freq1, freq2, freq3 = freq_dict[child_alleles[0]], freq_dict[grandparent_alleles[0]], freq_dict[grandparent_alleles[1]]
            return freq1, freq2, freq3
        if known_set == child_set:
            freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
            if len(ch_grand_inter) == 0:
                freq3, freq4 = freq_dict[grandparent_alleles[0]], freq_dict[grandparent_alleles[1]]
                return freq1, freq2, freq3, freq4
            if len(ch_grand_inter) == 1:
                freq1, freq2 = freq_dict[list(ch_grand_inter)], freq_dict[list(child_set - ch_grand_inter)]
                freq3 = freq_dict[list(grandparent_set - ch_grand_inter)]
                return freq1, freq2, freq3, 1
        freq1, freq2 = freq_dict[list(ch_known_inter)], freq_dict[list(child_set - ch_known_inter)]
        freq3, freq4 = freq_dict[grandparent_alleles[0]], freq_dict[grandparent_alleles[1]]
        if len(grandparent_set) == 2 and len(ch_grand_inter) == 1:
            freq3 = freq_dict[list(grandparent_set - child_set)]
            return freq1, freq2, freq3, freq4
