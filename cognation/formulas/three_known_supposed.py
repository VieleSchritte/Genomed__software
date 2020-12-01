from __future__ import unicode_literals
from .base import Formula, Calculations
from .one_known_supposed import OneKnownSupposedFormula
from .two_known_supposed import TwoKnownSupposedFormula


class ThreeKnownSupposed(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 5)
        supposed_alleles, child3_alleles, child2_alleles, child1_alleles, known_alleles = alleles
        supposed_set, child3_set, child2_set, child1_set, known_set = sets

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr = '-'
        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        if locus == 'AMEL':
            return self.make_result(locus, 1, dict_make_result)

        # If there are no intersections between children and parents, return lr = 0
        for i in range(0, 3):
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)

        c = Calculations()
        unique_genotype, repeat_genotype = c.get_repeat_unique(child1_alleles, child2_alleles, child3_alleles)
        if len(repeat_genotype) != 0:
            raw_values = [locus, '/'.join(known_alleles), '/'.join(repeat_genotype), '/'.join(supposed_alleles)]

            # all children genotypes are same, use OneKnownSupposed
            if len(unique_genotype) == 0:
                result = OneKnownSupposedFormula(Formula).calculate_relation(raw_values)
                lr = result['lr']
                return self.make_result(locus, lr, dict_make_result)

            # If there are two same child alleles, use TwoKnownSupposed
            raw_values.append('/'.join(unique_genotype))
            raw_values[-1], raw_values[-2] = raw_values[-2], raw_values[-1]
            result = TwoKnownSupposedFormula(Formula).calculate_relation(raw_values)
            lr = result['lr']
            return self.make_result(locus, lr, dict_make_result)

        common_set = set(child1_alleles + child2_alleles + child3_alleles + supposed_alleles + known_alleles)
        freq_dict = self.get_frequencies(locus, list(common_set))

        # ab ac bc ab ac/bc
        if len(common_set) == 3 and len(child1_set) == 2:
            freq1, freq2 = freq_dict[child1_alleles[0]], freq_dict[child1_alleles[1]]
            freq3 = freq_dict[list(child2_set - child1_set)[0]]
            lr = 2 * freq3 * (freq1 + freq2)
            return self.make_result(locus, lr, dict_make_result)

        # default describes all other cases
        freq1, freq2 = freq_dict[supposed_alleles[0]], freq_dict[supposed_alleles[1]]
        lr = 2 * freq1 * freq2
        return self.make_result(locus, lr, dict_make_result)
