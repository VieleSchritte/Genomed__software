from __future__ import unicode_literals
from .base import Formula, Calculations
from .one_known_supposed import OneKnownSupposedFormula


class TwoKnownSupposedFormula(Formula):
    def calculate_relation(self, raw_values):
        print('Two known supposed')
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 4)
        known_alleles, supposed_alleles, children_genotypes = alleles[0], alleles[1], alleles[2:]
        known_set, supposed_set, child1_set, child2_set = sets
        ch1ch2_inter = intersections[3]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        # If there are no intersections between children and parents, return lr = 0 and start counting mutations
        for i in range(1, 5):
            if i == 2 or i == 3:
                continue
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)
        common = []
        for genotype in alleles:
            common += genotype
        common_set, freq_dict = set(common), self.get_frequencies(locus, common)
        c = Calculations()

        # If children's genotypes are same, use OneKnownSupposedFormula
        raw_values = [locus, '/'.join(known_alleles), '/'.join(supposed_alleles)]
        lr = c.get_repeatable_lr(raw_values, children_genotypes, [OneKnownSupposedFormula(Formula)])
        if lr:
            return self.make_result(locus, lr, dict_make_result)

        # homozygous first child
        if len(child1_set) == 1:
            # aa ab ab an
            if (len(common_set) == 2 or len(common_set) == 3) and child2_set == known_set:
                freq = freq_dict[children_genotypes[0]]
                lr = c.F(freq)
                return self.make_result(locus, lr, dict_make_result)

            freq1, freq2 = freq_dict[supposed_alleles[0]], freq_dict[supposed_alleles[1]]
            lr = 2 * freq1 * freq2
            return self.make_result(locus, lr, dict_make_result)

        # Heterozygous first child
        else:
            # ab cd bd ac, ab ac an (n != b) bc
            if len(ch1ch2_inter) == 0 or list(child1_set - child2_set)[0] not in known_alleles:
                freq1, freq2 = freq_dict[supposed_alleles[0]], freq_dict[supposed_alleles[1]]
                lr = 2 * freq1 * freq2
                return self.make_result(locus, lr, dict_make_result)

            # ab ac ab ac/bc
            if list(ch1ch2_inter)[0] in known_alleles and child1_set == known_set:
                freq1, freq2 = freq_dict[known_alleles[0]], freq_dict[known_alleles[1]]
                freq3 = freq_dict[list(child2_set - child1_set)[0]]
                lr = 2 * freq3 * (freq1 + freq2)
                return self.make_result(locus, lr, dict_make_result)

            # ab ac bc an
            freq = freq_dict[list(ch1ch2_inter)[0]]
            lr = c.F(freq)
            return self.make_result(locus, lr, dict_make_result)
