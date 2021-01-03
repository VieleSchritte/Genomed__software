from __future__ import unicode_literals
from .base import Formula, Calculations


class BrotherFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 2)
        insp_set = sets[0]

        if self.is_gender_specific(locus):
            return self.result_gender_specific(locus, dict_make_result)

        c = Calculations()
        all_alleles, homo_counter = c.get_overall_alleles(alleles), c.homo_counter(sets)
        freq_dict = self.get_frequencies(locus, all_alleles)
        freq1, freq2, freq3, freq4 = self.freqs_order(sets, alleles, all_alleles, intersections[0], freq_dict)
        conditions_dict = {
            1: {
                (2, 1): (1, c.homo_refutation(freq1)),
                (1, 2): (2 * freq1 * freq2 * (2 * c.F(freq1) - freq1 * freq2) / (c.F(freq1) * c.F(freq2) - 2 * (freq1 * freq2) ** 2),
                         c.homo_refutation(freq1)),
                (2, 2): ((c.M(freq2, freq1)) ** 2, c.homo_refutation(freq1)),
                (1, 3): (4 * freq1 ** 2 * freq2 * freq3 / (c.F(freq2) * c.F(freq3) - 2 * (freq2 * freq3) ** 2),
                         c.homo_refutation(freq1))
            },
            2: {
                (1, 2): ((4 * freq1 * freq2 * c.F(freq1) - (2 * freq1 * freq2) ** 2) / (c.F(freq1)) ** 2,
                         c.hetero_refutation(freq1, freq2)),
                (0, 2): (1, c.hetero_refutation(freq1, freq2)),
                (0, 3): (2 * freq2 * freq3 * (c.F(freq1) - 2 * freq1 ** 2) / (c.F(freq1) * c.F(freq3) - 2 * (freq1 * freq3) ** 2),
                         c.hetero_refutation(freq1, freq2)),
                (1, 3): (2 * c.M(freq3, freq1) * c.M(freq3, freq2), c.hetero_refutation(freq1, freq2)),
                (0, 4): (2 * freq1 * freq2 * freq3 * freq4 / (c.F(freq3) * c.F(freq4) - 2 * (freq3 * freq4) ** 2),
                         c.hetero_refutation(freq1, freq2))
            }
        }
        if len(insp_set) in conditions_dict.keys():
            target_dict = conditions_dict[len(insp_set)]
            if (homo_counter, len(all_alleles)) in target_dict.keys():
                confirmation, refutation = target_dict[(homo_counter, len(all_alleles))]
                lr = confirmation / refutation
                return self.make_result(locus, lr, dict_make_result)

    @staticmethod
    def freqs_order(sets, alleles, all_alleles, intersection, freq_dict):
        insp_alleles, brother_alleles = alleles
        insp_set, brother_set = sets
        if len(all_alleles) == 1:
            return freq_dict[all_alleles[0]], 1, 1, 1
        if len(all_alleles) == 2:
            if len(intersection) == 0:
                return freq_dict[insp_alleles[0]], freq_dict[brother_alleles[0]], 1, 1
            if insp_set == brother_set:
                return freq_dict[insp_alleles[0]], freq_dict[insp_alleles[1]], 1, 1
            return freq_dict[list(intersection)[0]], freq_dict[list(set(all_alleles) - intersection)[0]], 1, 1
        if len(all_alleles) == 3:
            if len(insp_set) == 1:
                return freq_dict[insp_alleles[0]], freq_dict[brother_alleles[0]], freq_dict[brother_alleles[1]], 1
            if len(brother_set) == 1:
                return freq_dict[insp_alleles[0]], freq_dict[insp_alleles[1]], freq_dict[brother_alleles[0]], 1
            return freq_dict[list(intersection)[0]], freq_dict[list(insp_set - intersection)[0]], freq_dict[list(brother_set - intersection)[0]], 1
        return freq_dict[insp_alleles[0]], freq_dict[insp_alleles[1]], freq_dict[brother_alleles[0]], freq_dict[brother_alleles[1]]
