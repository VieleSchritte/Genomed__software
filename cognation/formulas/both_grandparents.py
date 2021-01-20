from __future__ import unicode_literals
from .base import Formula, Calculations


# FORMULA_TYPE_GRANDPARENT
class BothGrandparents(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        grandchild_alleles = alleles[0]
        grandchild_set, grandparents_sets = sets[0], [sets[1], sets[2]]
        if self.is_gender_specific(locus):
            return self.result_gender_specific(locus, dict_make_result)

        calc, c = Calculations(), Confirmations()
        freq_dict = self.get_frequencies(locus, grandchild_alleles)
        freq1, freq2 = freq_dict[grandchild_alleles[0]], freq_dict[grandchild_alleles[1]]

        if len(grandchild_set) == 1:
            refutation = calc.homo_refutation(freq1)
            confirmation = c.homo_gc_confirmation(freq1, grandchild_set, grandparents_sets, intersections)
        else:
            for grandparent_set in grandparents_sets:
                if len(grandparent_set & grandchild_set) == 1:
                    freq1 = freq_dict[list(grandparent_set & grandchild_set)[0]]
                    freq2 = freq_dict[list(grandchild_set - grandparent_set)[0]]
            inters_lens = []
            for intersection in intersections:
                inters_lens.append(len(intersection))
            refutation = calc.hetero_refutation(freq1, freq2)
            confirmation = c.hetero_gc_confirmation(grandparents_sets, grandchild_set, inters_lens, freq1, freq2)
        lr = confirmation / refutation
        return self.make_result(locus, lr, dict_make_result)


class Confirmations:
    @staticmethod
    def homo_gc_confirmation(freq, grandchild_set, grandparents_sets, intersections):
        calc = Calculations()
        # aa aa any
        counter = 0
        for grandparent_set in grandparents_sets:
            if grandparent_set == grandchild_set:
                return calc.F(freq)
            else:
                counter += 1
        intersections_counter = 0
        for i in range(2):
            if len(intersections[i]) != 0:
                intersections_counter += 1
        if counter == 2:
            # aa ab an, n != a
            int_counts_dict = {
                2: 0.75 * calc.F(freq),
                1: 0.5 * calc.F(freq),
                0: 0
            }
            for key in int_counts_dict.keys():
                if intersections_counter == key:
                    return int_counts_dict[key]
                else:
                    return 0

    @staticmethod
    def hetero_gc_confirmation(grandparents_sets, grandchild_set, inters_lens, freq1, freq2):
        calc = Calculations()
        homo_counter = 0
        if inters_lens[0:2] == [0, 0]:
            return 0

        for grandparent_set in grandparents_sets:
            if len(grandparent_set) == 1:
                homo_counter += 1

        h_counter_form_dict = {
            0: {
                (2, 2, 2): 0.75 * (calc.F(freq1) + calc.F(freq2)) - freq1 * freq2,  # ab ab ab
                (2, 0, 0): 0.5 * (calc.F(freq1) + calc.F(freq2)),  # ab ab any != an, bn
                (0, 2, 0): 0.5 * (calc.F(freq1) + calc.F(freq2)),  # ab ab any != an, bn
                (1, 1, 1): 0.5 * (calc.F(freq1) + calc.F(freq2) - freq1 * freq2),  # ab ac bc
                (2, 1, 1): 0.75 * calc.F(freq2) + 0.5 * (calc.F(freq1) - 2 * freq1 * freq2),  # ab ab ac
                (1, 2, 1): 0.75 * calc.F(freq2) + 0.5 * (calc.F(freq1) - 2 * freq1 * freq2),  # ab ab ac
                (1, 1, 2): 0.75 * calc.F(freq2),  # ab ac ac
                (1, 0, 1): 0.5 * calc.F(freq2),  # ab ac any != an, bn
                (0, 1, 1): 0.5 * calc.F(freq2),  # ab ac any != an, bn
                (1, 0, 0): 0.5 * calc.F(freq2),  # ab ac any != an, bn
                (0, 1, 0): 0.5 * calc.F(freq2)  # ab ac any != an, bn
            },
            1: {
                (1, 2, 1): calc.F(freq2) + 0.5 * (calc.F(freq1) - 2 * freq1 * freq2),  # ab aa bn, n != b
                (2, 1, 1): calc.F(freq2) + 0.5 * (calc.F(freq1) - 2 * freq1 * freq2),  # ab aa bn, n != b
                (1, 1, 0): calc.F(freq2) + 0.5 * (calc.F(freq1) - 2 * freq1 * freq2),  # ab aa bn, n != b
                (2, 0, 0): 0.5 * (calc.F(freq1) + calc.F(freq2)),  # ab ab any != an, bn
                (0, 2, 0): 0.5 * (calc.F(freq1) + calc.F(freq2)),  # ab ab any != an, bn
                (1, 1, 1): calc.F(freq2),  # ab aa an/cn, n != b
                (1, 0, 0): calc.F(freq2),  # ab aa an/cn, n != b
                (0, 1, 0): calc.F(freq2),  # ab aa an/cn, n != b
                (1, 0, 1): 0.5 * calc.F(freq2),  # ab ac any != an, bn
                (0, 1, 1): 0.5 * calc.F(freq2),  # ab ac any != an, bn
            },
            2: {
                (1, 1, 0): (freq1 + freq2) * (2 - (freq1 + freq2)),  # ab aa bb
                (0, 1, 0): calc.F(freq2),  # ab aa an/cn, n != b
                (1, 0, 0): calc.F(freq2),  # ab aa an/cn, n != b
                (1, 1, 1): calc.F(freq2)  # ab aa an/cn, n != b
            }
        }

        for h_counter in h_counter_form_dict.keys():
            if homo_counter == h_counter:
                target_dict = h_counter_form_dict[h_counter]
                for key in target_dict.keys():
                    moot_tuples = [(1, 0, 0), (0, 1, 0)]
                    if homo_counter == 1 and tuple(inters_lens) in moot_tuples:
                        for grandparent_set in grandparents_sets:
                            if len(grandparent_set) == 1 and len(grandparent_set & grandchild_set) == 0:
                                return 0.5 * calc.F(freq2)
                            else:
                                return 0
                    if key == tuple(inters_lens):
                        return target_dict[key]
                    else:
                        return 0
            else:
                return 0
