from __future__ import unicode_literals
from .base import Formula, Calculations
from .one_known_supposed import OneKnownSupposedFormula
from .two_known_supposed import TwoKnownSupposedFormula


class ThreeKnownSupposed(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 5)
        supposed_alleles, known_alleles, child3_alleles, child2_alleles, child1_alleles = alleles
        supposed_set, known_set, child3_set, child2_set, child1_set = sets

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr = '-'
        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        # If there are no intersections between children and parents, return lr = 0
        for i in range(1, 7):
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)

        c = Calculations()
        unique_genotype, repeat_genotype_raw, repeat_genotype_join = c.get_child13_genotypes(child1_alleles, child2_alleles, child3_alleles)

        # all children genotypes are same, use OneKnownSupposed
        if len(unique_genotype) == 0 and len(repeat_genotype_join) != 0:
            raw_values = [locus, '/'.join(child1_alleles), '/'.join(known_alleles), '/'.join(supposed_alleles)]
            result = OneKnownSupposedFormula(Formula).calculate_relation(raw_values)
            result['part4'] = '/'.join(child2_alleles)
            result['part5'] = '/'.join(child3_alleles)
            return result

        # If there are two same child alleles, use TwoKnownSupposed
        if len(unique_genotype) != 0 and len(repeat_genotype_raw) != 0:
            raw_values = [locus, '/'.join(supposed_alleles), '/'.join(known_alleles), '/'.join(unique_genotype), '/'.join(repeat_genotype_raw)]
            result = TwoKnownSupposedFormula(Formula).calculate_relation(raw_values)
            result['part5'] = '/'.join(repeat_genotype_join)
            return result

        common_set = set(child1_alleles + child2_alleles + child3_alleles + supposed_alleles + known_alleles)
        freq_dict = self.get_frequencies(locus, list(common_set))

        # special case aa ab bb ab ab
        if len(common_set) == 2:
            freq1, freq2 = freq_dict[list(common_set)[0]], freq_dict[list(common_set)[1]]
            lr = 2 * freq1 * freq2
            return self.make_result(locus, lr, dict_make_result)

        # special case ab ac bd ac bc
        if len(common_set) == 4:
            freq1, freq2 = freq_dict[child3_alleles[0]], freq_dict[child3_alleles[1]]
            lr = 2 * freq1 * freq2
            return self.make_result(locus, lr, dict_make_result)

        # Homozygous first child
        if len(child1_set) == 1:
            freq1, freq2 = freq_dict[supposed_alleles[0]], freq_dict[supposed_alleles[1]]
            lr = 2 * freq1 * freq2
            return self.make_result(locus, lr, dict_make_result)

        freq1, freq2, freq3 = freq_dict[child1_alleles[0]], freq_dict[child1_alleles[1]], freq_dict[list(common_set - child1_set)[0]]
        lr = 2 * freq3 * (freq1 + freq2)
        return self.make_result(locus, lr, dict_make_result)
