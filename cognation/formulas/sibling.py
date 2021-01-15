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
        print(freq_dict)
        if len(child_set) == 1:  # Homozygous child
            refutation = c.homo_refutation(freq_dict[list(child_set)[0]])
            freq1, freq2, freq3 = self.get_freq_order(sets, freq_dict, set_all)
            conf_dict = {
                (1, 3, 1): 1,  # aa aa aa
                (2, 2, 2): c.M(freq2, freq1),  # aa aa ab
                (2, 2, 1): 1,  # aa ab aa
                (2, 1, 2): c.F(freq1) / (c.F(freq1) + c.F(freq2) - 2 * freq1 * freq2),  # aa ab ab
                (3, 1, 2): c.M(freq3, freq1)  # aa ab ac/bc
            }
            if (len(all_alleles), homo_counter, len(sibling_set)) in conf_dict.keys():
                print((len(all_alleles), homo_counter, len(sibling_set)))
                lr = conf_dict[(len(all_alleles), homo_counter, len(sibling_set))] / refutation
                return self.make_result(locus, lr, dict_make_result)

        else:  # Heterozygous child
            refutation = c.hetero_refutation(freq_dict[child_alleles[0]], freq_dict[child_alleles[1]])
            freq1, freq2, freq3, freq4 = self.get_freq_order(sets, freq_dict, set_all)
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
            trial_key = (len(all_alleles), homo_counter, len(sibling_set), len(sc_inter))
            print(freq_dict)
            print(trial_key, freq2, freq4)
            print('child, parent, sibling: ', child_set, parent_set, sibling_set)
            if trial_key in answers.keys():
                target = answers[trial_key]
                if type(target) == dict:
                    confirmation = c.get_lr_from_cond_dict_short(target)
                    return self.make_result(locus, confirmation / refutation, dict_make_result)
                return self.make_result(locus, answers[trial_key] / refutation, dict_make_result)
            else:
                return self.make_result(locus, 0, dict_make_result)

    @staticmethod
    def get_freq_order(sets, freq_dict, set_all):
        parent_set, sibling_set, child_set = sets
        if len(child_set) == 1:  # DON'T CHANGE THESE NUMBERS!!!!
            freq1, freq2, freq3 = 3, 1, 1
            if len(parent_set) == 1:
                if len(sibling_set) == 2:  # only case aa aa ab is important: freqs order
                    freq1, freq2 = freq_dict[list(child_set)[0]], freq_dict[list(sibling_set - child_set)[0]]
                    return freq1, freq2, freq3
                return freq1, freq2, freq3  # case aa aa aa: no matter what freqs are: conf = 1
            if parent_set == sibling_set:  # aa ab ab
                freq1, freq2 = freq_dict[list(child_set)[0]], freq_dict[list(parent_set - child_set)[0]]
                return freq1, freq2, freq3
            if len(sibling_set) == 2:  # aa ab ac/bc
                freq1, freq2 = freq_dict[list(child_set)[0]], freq_dict[list(parent_set - child_set)[0]]
                freq3 = freq_dict[list(set_all - parent_set)[0]]
                return freq1, freq2, freq3
            return freq1, freq2, freq3  # case aa ab aa: no matter what freqs are: conf = 1

        # Heterozygous child
        freq1, freq2, freq3, freq4 = 1, 3, 3, 1  # DON'T CHANGE THESE NUMBERS!!!!
        if len(parent_set) == 1:
            if sibling_set == parent_set:  # ab aa aa
                freq1, freq2 = freq_dict[list(parent_set)[0]], freq_dict[list(child_set - parent_set)[0]]
                return freq1, freq2, freq3, freq4
            if sibling_set == child_set:  # ab aa ab
                return freq1, freq2, freq3, freq4
            freq1, freq2 = freq_dict[list(parent_set)[0]], freq_dict[list(child_set - parent_set)[0]]
            freq3 = freq_dict[list(set_all - child_set)[0]]
            return freq1, freq2, freq3, freq4  # ab aa ac
        if parent_set == child_set:
            if len(set_all) == 3:  # ab ab ac
                child_genotype = list(child_set)  # no matter what order of alleles 1 and 2 is
                freq1, freq2 = freq_dict[child_genotype[0]], freq_dict[child_genotype[1]]
                freq3 = freq_dict[list(set_all - child_set)[0]]
                return freq1, freq2, freq3, freq4
            return freq1, freq2, freq3, freq4  # in other cases no matter what freqs order is
        freq1, freq2 = freq_dict[list(child_set & parent_set)[0]], freq_dict[list(child_set - parent_set)[0]]
        freq3 = freq_dict[list(parent_set - child_set)[0]]
        if len(set_all) == 3:
            return freq1, freq2, freq3, freq4
        if len(set_all) == 4:  # cases ab ac ad/cd
            freq4 = freq_dict[list(set_all - child_set - parent_set)[0]]
            return freq1, freq2, freq3, freq4
        return freq1, freq2, freq3, freq4  # in other cases we just need this function to return anything
