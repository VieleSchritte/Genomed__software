from __future__ import unicode_literals
from cognation.formulas.base import Formula, Calculations


class GrandParentNo(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        child_alleles, grandparent_alleles = alleles[0], alleles[2]
        child_set, parent_set, grandparent_set = sets
        cp_inter, cg_inter = intersections[0], intersections[1]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)
        calc, c = Calculations(), Confirmations
        freq_dict = self.get_frequencies(locus, child_alleles)
        print(locus, freq_dict)
        freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
        if len(child_set) == 1:
            refutation = calc.homo_refutation(freq1)
            confirmation = c.homo_gc_confirmation(grandparent_set, cg_inter, freq1)
        else:
            refutation = calc.hetero_refutation(freq1, freq2)
            if len(cp_inter) == 1:
                freq1, freq2 = freq_dict[list(cp_inter)[0]], freq_dict[list(child_set - parent_set)[0]]
            confirmation = c.hetero_gc_confirmation(sets, intersections, freq1, freq2, grandparent_alleles)
            print(confirmation)
            print()
        lr = confirmation / refutation
        return self.make_result(locus, lr, dict_make_result)


class Confirmations:
    @staticmethod
    def homo_gc_confirmation(grandparent_set, cg_inter, freq):
        if len(grandparent_set) == 2 and len(cg_inter) == 1:
            return 0.5 * (1 + freq)
        else:
            return 1

    @staticmethod
    def hetero_gc_confirmation(sets, intersections, freq1, freq2, grandparent_alleles):
        child_set, parent_set, grandparent_set = sets
        cp_inter, cg_inter, pg_inter = intersections
        inters_lens = [len(cp_inter), len(cg_inter), len(pg_inter)]
        c = Calculations()
        homo_counter = c.homo_counter(sets)
        print(child_set, parent_set, grandparent_set)

        if parent_set == child_set:
            if len(cg_inter) == len(pg_inter) == 0:  # ab ab cc/cd
                return freq1 + freq2
            if len(grandparent_set) == 2 and len(pg_inter) == 1:  # ab ab ac
                return 0.5 * (1 + freq1 + freq2)
            return 1  # ab ab aa/ab
        else:
            print('went to else')
            if inters_lens == [1, 1, 0] and len(grandparent_set) == 1:  # ab an bb
                return 1
            if inters_lens == [1, 1, 1]:
                if homo_counter == 0:
                    return 0.5 * (1 + freq2)  # ab an bn n != b
                else:
                    return freq2
            other_answers = [
                [[1, 2, 1], [1, 1, 0], [2, 1, 0], [1, 2, 1]],
                [[1, 1, 2], [1, 0, 0]]
            ]
            if inters_lens in other_answers[0]:
                return 0.5 * (1 + freq2)  # ab an bn n != b
            if inters_lens in other_answers[1]:
                return freq2  # ab an any != bn

