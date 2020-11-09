from __future__ import unicode_literals
from .base import Formula
from .couple import CoupleFormula


class TwoCoupleFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 4)
        child2_alleles, child1_alleles, mother_alleles, father_alleles = alleles
        child2_set, child1_set, mother_set, father_set = sets
        ch1ch2_inter, ch2m_inter, ch2f_inter, ch1m_inter, ch1f_inter = intersections[0:5]

        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        common_set = set(child1_alleles + child2_alleles + mother_alleles + father_alleles)
        freq_dict = self.get_frequencies(locus, list(common_set))
        lr = 0

        # If there are no intersections between children and couple, return lr = 0 and start counting mutations
        for i in range(1, 5):
            if intersections[i] == 0:
                return self.make_result(locus, lr, dict_make_result)

        # both children genotypes are same, use CoupleFormula
        if child1_set == child2_set:
            raw_values = [locus, '/'.join(father_alleles), '/'.join(mother_alleles), '/'.join(child1_alleles)]
            result = CoupleFormula(Formula).calculate_relation(raw_values)
            lr = result['lr']
            return self.make_result(locus, lr, dict_make_result)

        # special case aa bb ab ab
        if len(common_set) == 2:
            freq1, freq2 = freq_dict[list(common_set)[0]], freq_dict[list(common_set)[1]]
            lr = (2 * freq1 * freq2) ** 2
            return self.make_result(locus, lr, dict_make_result)

        # Homozygous child1
        if len(child1_set) == 1:
            # aa ab an ab
            if len(ch1ch2_inter) != 0:
                freq1, freq2 = freq_dict[list(ch1m_inter)[0]], freq_dict[list(child2_set - child1_set)[0]]
                lr = 4 * (freq1 ** 2) * freq2 * (2 - freq1 - freq2)
                return self.make_result(locus, lr, dict_make_result)

            # aa bc ab ac
            freq1, freq2, freq3 = freq_dict[child1_alleles[0]], freq_dict[child2_alleles[0]], freq_dict[child2_alleles[1]]
            lr = 8 * (freq1 ** 2) * freq2 * freq3
            return self.make_result(locus, lr, dict_make_result)

        # Heterozygous child1
        else:
            # case ab ac an/ab bc/ac
            if len(ch1ch2_inter) != 0:
                freq1 = freq_dict[list(ch1ch2_inter)[0]]
                freq2, freq3 = freq_dict[list(child1_set - child2_set)[0]], freq_dict[list(child2_set - child1_set)[0]]
                lr = 4 * freq1 * freq2 * freq3 * (2 + freq1)
                return self.make_result(locus, lr, dict_make_result)

            # case ab cd ad/ac bc/bd
            lr = 16
            for i in range(len(list(common_set))):
                freq = freq_dict[list(common_set)[i]]
                lr *= freq
            return self.make_result(locus, lr, dict_make_result)
