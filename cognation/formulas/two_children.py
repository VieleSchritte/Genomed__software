from __future__ import unicode_literals
from .base import Formula, Calculations
from .parent import ParentFormula


class TwoChildrenFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, part_alleles, part_sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        child2_alleles, child1_alleles, parent_alleles = part_alleles
        child2_set, child1_set, parent_set = part_sets
        ch1ch2_inter, ch2p_inter, ch1p_inter = intersections

        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

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

        freq1, freq2 = freq_dict[child1_alleles[0]], freq_dict[child1_alleles[1]]
        freq3, freq4 = freq_dict[child2_alleles[0]], freq_dict[child2_alleles[1]]
        # Homozygous 1st child
        if len(child1_set) == 1:
            if len(ch1ch2_inter) == 0:
                # aa bb ab
                if len(child2_set) == 1:
                    lr = 2 * freq1 * freq2
                    return self.make_result(locus, lr, dict_make_result)

                # aa bc ab/ac
                else:
                    lr = 2 * freq1 * (freq2 + freq3)
                    return self.make_result(locus, lr, dict_make_result)

            lr = c.F(freq1)
            return self.make_result(locus, lr, dict_make_result)

        # Heterozygous 1st child
        else:
            if len(ch1ch2_inter) == 0:
                # case ab cc ac/bc
                if len(child2_set) == 1:
                    lr = 2 * freq3 * (freq1 + freq2)
                    return self.make_result(locus, lr, dict_make_result)

                # case ab cd ac/ad/bc/bd
                else:
                    lr = 2 * (freq3 + freq4) * (freq1 + freq2)
                    return self.make_result(locus, lr, dict_make_result)

            # case ab ac an/bc
            else:
                freq1 = freq_dict[list(ch1ch2_inter)[0]]
                freq2, freq3 = freq_dict[list(child1_set - child2_set)[0]], freq_dict[list(child2_set - child1_set)[0]]
                lr = c.F(freq1) + 2 * freq2 * freq3
                return self.make_result(locus, lr, dict_make_result)
