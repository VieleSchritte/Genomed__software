from __future__ import unicode_literals
from .base import Formula, Calculations
from .one_known_supposed import OneKnownSupposedFormula


class TwoKnownSupposedFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 4)
        known_alleles, supposed_alleles, children_genotypes = alleles[0], alleles[1], alleles[2:]
        known_set, supposed_set, children_sets = sets[0], sets[1], sets[2:]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        # If there are no intersections between children and parents, return lr = 0 and start counting mutations
        for child_set in children_sets:
            if child_set & known_set == 0 or child_set & supposed_set == 0:
                return self.make_result(locus, 0, dict_make_result)

        # If children's genotypes are same, use OneKnownSupposedFormula
        raw_values = [locus, '/'.join(known_alleles), '/'.join(supposed_alleles)]
        c = Calculations()
        lr = c.get_repeatable_lr(raw_values, children_genotypes, [OneKnownSupposedFormula(Formula)])
        if lr:
            return self.make_result(locus, lr, dict_make_result)
        else:
            possible_genotypes = c.get_possible_genotypes(c.get_children_set(children_genotypes), children_genotypes, known_set)
            common_pos = []
            for genotype in possible_genotypes:
                common_pos += genotype
            freq_dict = self.get_frequencies(locus, common_pos)

            if len(possible_genotypes) == 0:
                for child_set in children_sets:
                    if len(child_set) == 1:
                        freq = freq_dict[list(child_set)[0]]
                        return c.F(freq)
            if len(possible_genotypes) == 1:
                freq1, freq2 = freq_dict[possible_genotypes[0][0]], freq_dict[possible_genotypes[0][1]]
                return 2 * freq1 * freq2

            freq3 = freq_dict[set(possible_genotypes[0]) & set(possible_genotypes[1])]
            other_freqs = []
            for i in range(len(possible_genotypes)):
                other_freqs.append(freq_dict[list(set(possible_genotypes[i]) - set(possible_genotypes[i - 1]))[0]])
            return 2 * freq3 * (other_freqs[0] + other_freqs[1])
