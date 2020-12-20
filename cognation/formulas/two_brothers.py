from __future__ import unicode_literals
from .base import Formula, Calculations
from .brother import BrotherFormula


class TwoBrothersFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        inspected_alleles, brothers_genotypes = alleles[0], alleles[1:]
        inspected_set, brothers_sets = sets[0], sets[1:]
        inters_lens = []
        for intersection in intersections:
            inters_lens.append(len(intersection))
        inters_lens = tuple(inters_lens)
        c = Calculations()
        all_alleles = c.get_overall_alleles(alleles)

        print(locus)
        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)
        zero_conditions = [
            inters_lens == (1, 0, 0) and inspected_set not in brothers_sets,  # aa ab cc/cd, ab aa cc/cd
            inters_lens == (0, 1, 0) and inspected_set not in brothers_sets,  # aa ab cc/cd
            inters_lens == (1, 1, 1) and len(all_alleles) == 4,  # ab ac ad
            inters_lens[0] == 0 and inspected_set not in brothers_sets  # ab any != an/bn any != an/bn
        ]
        for condition in zero_conditions:
            if condition:
                return self.make_result(locus, 0, dict_make_result)

        raw_values = [locus, '/'.join(inspected_alleles)]
        lr = c.get_repeatable_lr(raw_values, brothers_genotypes, [BrotherFormula(Formula)])
        if lr:
            return self.make_result(locus, lr, dict_make_result)
        else:
            freq_dict = self.get_frequencies(locus, all_alleles)
            freq3 = 1
            conf = Confirmations()
            if len(intersections[-1]) == 1:
                freq3 = freq_dict[list(intersections[-1])[0]]
            if len(inspected_set) == 1:
                freq1, freq2 = freq_dict[inspected_alleles[0]], freq_dict[list(set(all_alleles) - intersections[-1] - inspected_set)[0]]
                refutation = c.homo_refutation(freq_dict[inspected_alleles[0]])
                confirmation = conf.homo_confirmation(inspected_set, brothers_sets, freq1, freq2, inters_lens)
            else:
                freq1, freq2 = 1, 1
                for bro_set in brothers_sets:
                    if len(bro_set & inspected_set) != 0:
                        freq1 = freq_dict[list(bro_set & inspected_set)[0]]
                        freq2 = freq_dict[list(inspected_set - bro_set & inspected_set)[0]]
                freqs = [freq1, freq2, freq3]
                refutation = c.hetero_refutation(freq_dict[inspected_alleles[0]], freq_dict[inspected_alleles[1]])
                confirmation = conf.hetero_confirmation(inters_lens, intersections, sets, freqs)
            lr = confirmation / refutation
            return self.make_result(locus, lr, dict_make_result)


class Confirmations:
    @staticmethod
    def homo_confirmation(inspected_set, brothers_sets, freq1, freq2, inters_lens):
        print('called homo')
        if inspected_set in brothers_sets:  # aa aa any
            print('right')
            return 1

        homo_dict = {
            (1, 1, 1): 2 * freq1 / (2 - freq1 + 2 * freq2),  # aa ab ac
            (1, 0, 1): 2 * freq1 / (2 + freq2)  # aa ab bc
        }
        for key in homo_dict.keys():
            if key == inters_lens:
                return homo_dict[key]

    @staticmethod
    def hetero_confirmation(inters_lens, intersections, sets, freqs):
        c = Calculations()
        hetero_counter = c.hetero_counter(sets)
        freq1, freq2, freq3 = freqs
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
                    return hetero_dict[key]
            if inters_lens == key:
                target_dict = hetero_dict[key]
                if hetero_counter in target_dict.keys():
                    return target_dict[hetero_counter]
