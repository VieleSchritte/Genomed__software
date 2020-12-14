from __future__ import unicode_literals
from .base import Formula, Calculations
from .parent import ParentFormula


class TwoChildrenFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, part_alleles, part_sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        parent_alleles, child1_alleles, child2_alleles = part_alleles
        child1_set, child2_set = part_sets[1:]
        ch1p_inter, ch2p_inter, ch1ch2_inter = intersections

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        # If children's genotypes are same, use ParentFormula
        if child1_set == child2_set:
            raw_values = [locus, '/'.join(parent_alleles), '/'.join(child1_alleles)]
            result = ParentFormula(Formula).calculate_relation(raw_values)
            lr = result['lr']
            return self.make_result(locus, lr, dict_make_result)

        c = Calculations()
        common_set = set(child1_alleles + child2_alleles + parent_alleles)
        freq_dict = self.get_frequencies(locus, list(common_set))
        lr = 0

        # If there are no intersections, return lr = 0 and start counting mutations
        for i in range(1, 2):
            if len(intersections[i]) == 0:
                return self.make_result(locus, lr, dict_make_result)

        if len(child1_set) == len(child2_set) == 1:
            # aa an an
            if ch1p_inter == ch2p_inter:
                lr = c.F(freq_dict[list(ch1p_inter)[0]])
                return self.make_result(locus, 1 / lr, dict_make_result)
            # aa bb ab
            else:
                lr = 2 * freq_dict[child1_alleles[0]] * freq_dict[child2_alleles[0]]
                return self.make_result(locus, 1 / lr, dict_make_result)

        # aa bc ab/ac and ab cc ac/bc
        if len(child1_set) != len(child2_set):
            children_sets = part_sets[1:]
            homo_allele, hetero_alleles = [], []
            for i in range(len(children_sets)):
                child_set = children_sets[i]
                if len(child_set) == 1:
                    homo_allele = list(child_set)
                else:
                    hetero_alleles = list(child_set)
            lr = 2 * freq_dict[homo_allele[0]] * (freq_dict[hetero_alleles[0]] + freq_dict[hetero_alleles[1]])
            return self.make_result(locus, 1 / lr, dict_make_result)

        # ab ac an/bc
        if len(child1_set) == len(child2_set) == 2 and len(ch1ch2_inter) != 0:
            freq1 = freq_dict[list(ch1ch2_inter)[0]]
            freq2, freq3 = freq_dict[list(child1_set - child2_set)[0]], freq_dict[list(child2_set - child1_set)[0]]
            lr = c.F(freq1) + 2 * freq2 * freq3
            return self.make_result(locus, 1 / lr, dict_make_result)

        freq1, freq2 = freq_dict[child1_alleles[0]], freq_dict[child1_alleles[1]]
        freq3, freq4 = freq_dict[child2_alleles[0]], freq_dict[child2_alleles[1]]
        lr = 2 * (freq1 + freq2) * (freq3 + freq4)
        return self.make_result(locus, 1 / lr, dict_make_result)
