from __future__ import unicode_literals
from .base import Formula, Calculations
from .brother import BrotherFormula


class TwoBrothersFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        inspected_alleles, brothers_genotypes = alleles[0], alleles[1:]
        inspected_set, brothers_sets = sets[0], sets[1:]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)
        c = Calculations()
        raw_values = [locus, '/'.join(inspected_alleles)]
        lr = c.get_repeatable_lr(raw_values, brothers_genotypes, [BrotherFormula])
        if lr:
            return self.make_result(locus, lr, dict_make_result)
        else:
            all_alleles = c.get_overall_alleles(alleles)
            freq_dict = self.get_frequencies(locus, all_alleles)
            if len(inspected_set) == 1:
                freq1, freq2, freq3 = freq_dict[inspected_alleles[0]], 1, 1
                homo_refutation = c.homo_refutation(freq_dict[inspected_alleles[0]])
                if inspected_set in brothers_sets:  # aa aa any
                    return self.make_result(locus, 1 / homo_refutation, dict_make_result)
                for bro_set in brothers_sets:
                    if len(bro_set & inspected_set) != 0:
                        bro_genotype = list(bro_set)
                        freq2, freq3 = freq_dict[bro_genotype[0]], freq_dict[bro_genotype[1]]

                homo_dict = {
                    [1, 1, 1]: 2 * freq1 / (2 - freq1 + 2 * freq2),
                    [1, 0, 1]: 2 * freq1 / (2 + freq2)
                }
                if intersections in homo_dict.keys():
                    return self.make_result(locus, homo_dict[intersections] / homo_refutation, dict_make_result)
                else:
                    return self.make_result(locus, 0, dict_make_result)

            hetero_refutation = c.hetero_refutation(freq_dict.keys()[0], freq_dict.keys()[1])
            zero_conditions = [
                intersections == [1, 0, 0],
                intersections[0] == 0,
                intersections == [1, 1, 1] and len(all_alleles) == 4
            ]
            for condition in zero_conditions:
                if condition:
                    return self.make_result(locus, 0, dict_make_result)
            hetero_counter = c.hetero_counter(sets)
            freq1, freq2, freq3 = 1, 1, 1
            for bro_set in brothers_sets:
                if len(bro_set & inspected_set) != 0:
                    freq1 = freq_dict[list(bro_set & inspected_set)[0]]
                    freq2 = freq_dict[list(inspected_set - bro_set & inspected_set)[0]]
            if len(intersections[-1]) == 1:
                freq3 = freq_dict[list(intersections[-1])[0]]

            hetero_dict = {
                (1, 1, 1): {
                    2: freq2 / (2 - freq1),  # ab aa ac
                    3: 2 * freq1 / (2 + freq3)  # ab ac bc
                },
                (1, 1, 0): {
                    1: 1,  # ab aa bn
                    2: 1,  # ab aa bn
                    3: 0.5  # ab ac bd
                },
                (1, 0, 1): {
                    2: 2 * freq2 * freq3 / c.F(freq1),  # ab ac cc
                    3: 2 * freq2 / (2 + freq3)  # ab ac cd
                },
                2: 1  # ab ab any / ab aa bn (case ab aa ab)
            }
            for key in hetero_dict.keys():
                if type(key) == int:
                    if intersections[0] == 2:
                        confirmation = hetero_dict[key]
                if intersections == key:
                    target_dict = hetero_dict[key]
                    if hetero_counter in target_dict.keys():
                        confirmation = target_dict[hetero_counter]
            lr = confirmation / hetero_refutation
            return self.make_result(locus, lr, dict_make_result)
