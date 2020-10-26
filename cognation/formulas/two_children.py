from __future__ import unicode_literals
from .base import Formula
from .base import Calculations
from .parent import ParentFormula


class TwoChildren(Formula):
    def calculate_relation(self, raw_values):
        (locus, part_alleles, part_sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        child1_alleles, child2_alleles, parent_alleles = part_alleles
        child1_set, child2_set, parent_set = part_sets
        ch1p_intersection, ch2p_intersection, ch1ch2_intersection = intersections

        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        if child1_set == child2_set:
            return ParentFormula(Formula).calculate_relation(raw_values)

        c = Calculations()
        freq_dict = self.get_frequencies(locus, child1_alleles + child2_alleles + parent_alleles)
        lr = 0

        if len(ch1p_intersection) >= 1 and len(ch2p_intersection) >= 1:
            # Homozygous 1st child
            if len(child1_set) == 1:
                # case aa an an
                if ch2p_intersection == ch1p_intersection:
                    lr = c.F(freq_dict[child1_alleles[0]])
                    return self.make_result(locus, lr, dict_make_result)

                common_set = set(child1_alleles + child2_alleles + parent_set)
                # case aa bb ab
                if len(common_set) == 2:
                    freq1, freq2 = freq_dict[list(common_set)[0]], freq_dict[list(common_set)[1]]
                    lr = 2 * freq1 * freq2
                    return self.make_result(locus, lr, dict_make_result)

                # case aa bc ab/ac
                if len(common_set) == 3:
                    freq1 = freq_dict[list(ch1p_intersection)]
                    unique_allele1 = c.get_unique_allele(parent_alleles, child1_alleles)
                    freq2 = freq_dict[unique_allele1]
                    unique_allele2 = c.get_unique_allele(list(common_set), [unique_allele1, list(ch1p_intersection)])
                    freq3 = freq_dict[unique_allele2]
                    lr = 2 * freq1 * (freq2 + freq3)
                    return self.make_result(locus, lr, dict_make_result)

            # Heterozygous 1st child
            else:
                if len(child2_set) == 1:
                    
