from __future__ import unicode_literals
from cognation.formulas.base import Formula, Calculations


class GrandParentNo(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        child_alleles = alleles[0]
        child_set, parent_set, grandparent_set = sets
        cp_inter, cg_inter = intersections[0], intersections[1]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)
        calc, c = Calculations(), Confirmations
        freq_dict = self.get_frequencies(locus, child_alleles)
        freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
        if len(child_set) == 1:
            refutation = calc.homo_refutation(freq1)
            confirmation = c.homo_gc_confirmation(grandparent_set, cg_inter, freq1)
        else:
            refutation = calc.hetero_refutation(freq1, freq2)
            if len(cp_inter) == 1:
                freq1, freq2 = freq_dict[list(cp_inter)[0]], freq_dict[list(child_set - parent_set)[0]]
            confirmation = c.hetero_gc_confirmation(intersections, sets, freq1, freq2)
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
    def hetero_gc_confirmation(intersections, sets, freq1, freq2):
        child_set, parent_set, grandparent_set = sets
        if child_set == parent_set:
            if len(intersections[1]) == 0:  # ab ab cc/cd
                return freq1 + freq2
            if grandparent_set == parent_set or len(grandparent_set) == 1:  # ab ab aa/ab
                return 1
            return 0.5 * (1 + freq1 + freq2)

        target_allele = child_set - parent_set
        if grandparent_set == target_allele:  # ab an bb
            return 1
        if list(target_allele)[0] in list(grandparent_set):  # ab an bn, n != b
            return 0.5 * (1 + freq2)
        return freq2  # ab an any != bn
