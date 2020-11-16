from __future__ import unicode_literals
from .base import Formula, Calculations
from .couple import CoupleFormula
from .two_couple import TwoCoupleFormula


class ThreeCoupleFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 5)
        child3_alleles, child2_alleles, child1_alleles, mother_alleles, father_alleles = alleles
        child3_set, child2_set, child1_set, mother_set, father_set = sets
        mother_father_inter = intersections[-1]

        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        for i in range(2, 9):
            if i == 4:
                continue
            if intersections[i] == 0:
                return self.make_result(locus, 0, dict_make_result)

        c = Calculations()
        unique_genotype, repeat_genotype = c.get_repeat_unique(child1_alleles, child2_alleles, child3_alleles)
        if len(repeat_genotype) != 0:
            raw_values = [locus, '/'.join(father_alleles), '/'.join(mother_alleles), '/'.join(repeat_genotype)]

            if len(unique_genotype) == 0:
                result = CoupleFormula(Formula).calculate_relation(raw_values)
                lr = result['lr']
                return self.make_result(locus, lr, dict_make_result)

            raw_values.append('/'.join(unique_genotype))
            result = TwoCoupleFormula(Formula).calculate_relation(raw_values)
            lr = result['lr']
            return self.make_result(locus, lr, dict_make_result)

        common_set = set(child1_alleles + child2_alleles + mother_alleles + father_alleles)
        alleles_list = list(common_set)
        freq_dict = self.get_frequencies(locus, list(common_set))

        if len(common_set) == 2:
            freq1, freq2 = freq_dict[list(common_set)[0]], freq_dict[list(common_set)[1]]
            lr = (2 * freq1 * freq2) ** 2
            return self.make_result(locus, lr, dict_make_result)

        if len(common_set) == 4:
            lr = 8
            for i in range(len(alleles_list)):
                lr *= alleles_list[i]
            return self.make_result(locus, lr, dict_make_result)

        if len(child1_set) == 1:
            freq1 = freq_dict[list(mother_father_inter)[0]]
            freq2, freq3 = freq_dict[list(mother_set - father_set)[0]], freq_dict[list(father_set - mother_set)[0]]
            lr = 8 * freq1 ** 2 * freq2 * freq3
            return self.make_result(locus, lr, dict_make_result)

        composition = 1
        sum_fr = 0
        for i in range(len(alleles_list)):
            composition *= alleles_list[i]
            sum_fr += alleles_list[i]
        lr = 8 * composition * sum_fr
        return self.make_result(locus, lr, dict_make_result)
