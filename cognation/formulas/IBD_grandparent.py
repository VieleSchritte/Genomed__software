from __future__ import unicode_literals
from .base import Formula


# FORMULA_TYPE_GRANDPARENT
class IBDGrandParent(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 2)
        gc_alleles = alleles[0]
        gc_set, gp_set = sets
        intersection, len_inter = intersections[0], len(intersections[0])

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        freq_dict = self.get_frequencies(locus, gc_set)
        freq1, freq2 = freq_dict[gc_alleles[0]], freq_dict[gc_alleles[1]]
        if len(gc_set) == 2 and len_inter == 1:
            freq1, freq2 = freq_dict[list(intersection)[0]], freq_dict[list(gc_set - intersection)[0]]

        lrs_dict = {
            1: {
                len_inter != 0 and gp_set == gc_set: 0.5 + 0.5 / freq1,  # aa aa
                len_inter != 0 and gp_set != gc_set: 0.5 + 0.25 / freq1,  # aa ab
                len_inter == 0: 0.5  # aa bb/bc
            },
            2: {
                len_inter != 0 and gp_set == gc_set: 0.5 + 0.25 * (freq1 + freq2) / (2 * (freq1 * freq2)),  # ab ab
                len_inter != 0 and gp_set != gc_set: 0.5 + 0.25 / (2 * freq1),  # ab ac
                len_inter == 0: 0.5  # ab cd
            }
        }
        for key in lrs_dict.keys():
            if key == len(gc_set):
                target_dict = lrs_dict[key]
                for condition in target_dict.keys():
                    if condition:
                        return self.make_result(locus, target_dict[condition], dict_make_result)
