from __future__ import unicode_literals
from cognation.formulas.base import Formula, Calculations


class SiblingFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        child_alleles = alleles[-1]
        parent_set, sibling_set, child_set = sets
        sp_inter, cp_inter, sc_inter = intersections

        if self.is_gender_specific(locus):
            return self.result_gender_specific(locus, dict_make_result)
        if len(sp_inter) == 0 or len(cp_inter) == 0:
            return self.make_result(locus, 0, dict_make_result)
        c = Calculations()
        all_alleles, set_all, homo_counter = c.get_overall_alleles(alleles), set(c.get_overall_alleles(alleles)), c.homo_counter(sets)
        freq_dict = self.get_frequencies(locus, all_alleles)
        if len(child_set) == 1:  # Homozygous child
            freq1, freq2, freq3 = self.get_freq_order(child_set, all_alleles, freq_dict, child_alleles, parent_set, set_all, cp_inter)
            refutation = c.homo_refutation(freq1)
            conf_dict = {
                (1, 3, 1): 1,
                (2, 2, 2): c.M(freq2, freq1),
                (2, 2, 1): 1,
                (2, 1, 2): c.F(freq1) / (c.F(freq1) + c.F(freq2) - 2 * freq1 * freq2),
                (3, 1, 2): c.M(freq3, freq1)
            }
            if (len(all_alleles), homo_counter, len(sibling_set)) in conf_dict.keys():
                lr = conf_dict[(len(all_alleles), homo_counter, len(sibling_set))] / refutation
                return self.make_result(locus, lr, dict_make_result)

        else:  # Heterozygous child
            freq1, freq2, freq3, freq4 = self.get_freq_order(child_set, all_alleles, freq_dict, child_alleles, parent_set, set_all, cp_inter)
            refutation = c.hetero_refutation(freq1, freq2)
            answers = {
                (2, 2, 1, 1): c.M(freq1, freq2),  # ab aa aa
                (2, 1, 2, 2): 1,  # ab aa ab
                (3, 1, 2, 1): c.M(freq3, freq1),  # ab aa ac
                (2, 1, 1, 1): 1,  # ab ab aa
                (2, 0, 2, 2): 1,  # ab ab ab
                (3, 0, 2, 1): {
                    parent_set == child_set: c.M(freq3, freq1) + c.M(freq3, freq2),  # ab ab ac
                    parent_set == sibling_set: 2 * freq2 / (2 - freq1 - freq3),  # ab ac ac
                    parent_set != child_set != sibling_set: 1  # ab ac bc
                },
                (3, 0, 2, 2): 1,  # ab ac ab
                (3, 1, 1, 1): c.M(freq1, freq3),  # ab ac aa
                (4, 0, 2, 1): c.M(freq4, freq2),  # ab ac ad
                (3, 1, 1, 0): c.M(freq3, freq2),  # ab ac cc
                (4, 0, 2, 0): c.M(freq4, freq2)  # ab ac cd
            }
            for key in answers.keys():
                if (len(all_alleles), homo_counter, len(sibling_set), len(sc_inter)) == key:
                    if type(answers[key]) == dict:
                        target_dict = answers[key]
                        confirmation = c.get_lr_from_cond_dict_short(target_dict)
                        return self.make_result(locus, confirmation / refutation, dict_make_result)
                    return self.make_result(locus, answers[key] / refutation, dict_make_result)

    @staticmethod
    def get_freq_order(child_set, all_alleles, freq_dict, child_alleles, parent_set, set_all, cp_inter):
        c = Calculations()
        minus_child, minus_parent = list(set_all - child_set)[0], list(set_all - parent_set)[0]
        minus_both = list(set_all - parent_set - child_set)[0]
        answers = {
            1: {
                len(all_alleles) == 2: (freq_dict[child_alleles[0]], freq_dict[minus_child], 1),
                len(all_alleles) == 3: (freq_dict[child_alleles[0]], freq_dict[list(parent_set - child_set)[0]], freq_dict[minus_parent])
            },
            2:  {
                2: {
                    len(cp_inter) == 1: (freq_dict[list(cp_inter)[0]], freq_dict[list(set_all - cp_inter)[0]], 1, 1),
                    len(cp_inter) == 2: (freq_dict[child_alleles[0]], freq_dict[child_alleles[1]], 1, 1)
                },
                3: {
                    len(cp_inter) == 1: (freq_dict[list(cp_inter)[0]], freq_dict[list(child_set - cp_inter)[0]], freq_dict[minus_child]),
                    len(cp_inter) == 2: (freq_dict[child_alleles[0]], freq_dict[child_alleles[1]], freq_dict[list(minus_child)[0]])
                }
            }
        }
        for key in answers.keys():
            if len(child_set) == key:
                target_dict = answers[key]
                if type(list(target_dict.keys())[0]) != int:
                    return c.get_lr_from_cond_dict_short(target_dict)
                for target_key in target_dict.keys():
                    if len(all_alleles) == key:
                        return c.get_lr_from_cond_dict_short(target_dict[target_key])

        freq1 = freq_dict[child_alleles[0]]
        freq2, freq4 = freq_dict[list(child_set - cp_inter)[0]], freq_dict[minus_both]
        return freq1, freq2, 1, freq4
