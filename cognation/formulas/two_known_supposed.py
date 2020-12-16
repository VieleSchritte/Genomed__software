from __future__ import unicode_literals
from .base import Formula, Calculations
from .one_known_supposed import OneKnownSupposedFormula


class TwoKnownSupposedFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 4)
        known_alleles, supposed_alleles, children_genotypes = alleles[0], alleles[1], alleles[2:]
        known_set, children_sets = sets[0], sets[2:]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        # If there are no intersections between children and parents, return lr = 0 and start counting mutations
        for i in range(1, 3):
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)

        # If children's genotypes are same, use OneKnownSupposedFormula
        raw_values = [locus, '/'.join(known_alleles), '/'.join(supposed_alleles)]
        c = Calculations()
        lr = c.get_repeatable_lr(raw_values, children_genotypes, [OneKnownSupposedFormula(Formula)])
        if lr:
            return self.make_result(locus, lr, dict_make_result)
        else:
            children_alleles = c.get_overall_alleles(children_genotypes)
            freq_dict = self.get_frequencies(locus, children_alleles)
            possible_parents_sets = c.get_possible_genotypes(children_alleles, children_genotypes, known_set)
            if type(possible_parents_sets) == set:  # cases where lr = c.F(Pa)
                lr = c.F(freq_dict[list(possible_parents_sets)[0]])
                return self.make_result(locus, 1 / lr, dict_make_result)

            for child_set in children_sets:  # the rest of homozygosity is default cases
                if len(child_set) == 1:
                    possible_genotype = list(possible_parents_sets[0])
                    freq1, freq2 = freq_dict[possible_genotype[0]], freq_dict[possible_genotype[1]]
                    lr = 2 * freq1 * freq2
                    return self.make_result(locus, 1 / lr, dict_make_result)
                else:
                    if len(possible_parents_sets) == 1:
                        possible_genotype = list(possible_parents_sets[0])
                        freq1, freq2 = freq_dict[possible_genotype[0]], freq_dict[possible_genotype[1]]
                        lr = 2 * freq1 * freq2
                        return self.make_result(locus, 1 / lr, dict_make_result)
                    if len(possible_parents_sets) == 2:
                        allele3 = possible_parents_sets[0] & possible_parents_sets[1]
                        allele1, allele2 = possible_parents_sets[0] - allele3, possible_parents_sets[1] - allele3
                        freq3 = freq_dict[list(allele3)[0]]
                        freq1, freq2 = freq_dict[list(allele1)[0]], freq_dict[list(allele2)[0]]
                        lr = 2 * freq3 * (freq1 + freq2)
                        return self.make_result(locus, 1 / lr, dict_make_result)
