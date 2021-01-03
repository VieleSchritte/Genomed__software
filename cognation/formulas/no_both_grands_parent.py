from __future__ import unicode_literals
from .base import Formula, Calculations


# FORMULA_TYPE_GRANDPARENT
class NoBothGrandsParent(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 4)
        child_alleles, child_set = alleles[0], sets[0]

        if self.is_gender_specific(locus):
            return self.result_gender_specific(locus, dict_make_result)
        if intersections[0] == 0:
            return self.make_result(locus, 0, dict_make_result)

        freq_dict = self.get_frequencies(locus, child_alleles)
        freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
        calc, c = Calculations(), Confirmations()

        situations = {
            1: [calc.homo_refutation(freq1), c.homo_confirmation(sets, child_set)],
            2: [calc.hetero_refutation(freq1, freq2), c.hetero_confirmation(sets, child_alleles)]
        }
        refutation, confirmation = 0, 0
        for key in situations.keys():
            if len(child_set) == key:
                refutation, confirmation = situations[key]
        return self.make_result(locus, confirmation / refutation, dict_make_result)


class Confirmations:
    @staticmethod
    def homo_confirmation(sets, child_set):
        c, counter = Calculations(), 0
        grandparents_sets = [sets[2], sets[3]]
        child_allele = list(child_set)[0]
        for grandparent_set in grandparents_sets:
            if grandparent_set == child_set:  # aa an aa any
                return 1
            if child_allele in list(grandparent_set):
                counter += 1
        answers = {
            counter == 2: 0.75,  # aa an ab an, n != b
            counter != 2: 0.5  # aa an ab any != an
        }
        return c.get_lr_from_cond_dict_short(answers)

    @staticmethod
    def hetero_confirmation(sets, child_alleles):
        child_set, parent_set, grandparents_sets = sets[0], sets[1], [sets[2], sets[3]]
        if parent_set == child_set:
            occurrences = [0, 0]
            for grandparent_set in grandparents_sets:
                if grandparent_set == child_set or len(grandparent_set) == 1 and len(grandparent_set & child_set) == 1:
                    return 1
                for i in range(2):
                    if child_alleles[i] in grandparent_set:
                        occurrences[i] += 1
            occurrences_dict = {
                (2, 0): 0.75,  # ab ab an an n != b
                (0, 2): 0.75,  # ab ab bn bn n != b
                (1, 1): 0.75,  # ab ab an/bn an/bn n != b
                (1, 0): 0.5,  # ab ab an any != an, bn, n != b
                (0, 1): 0.5  # ab ab bn any != an, bn, n != b
            }
            for key in occurrences_dict.keys():
                if tuple(occurrences) == key:
                    return occurrences_dict[key]
        else:
            counter = 0
            unique_child_allele = child_set - parent_set
            for grandparent_set in grandparents_sets:
                if len(grandparent_set) == 1 and grandparent_set == unique_child_allele:
                    return 1
                if unique_child_allele in list(grandparent_set):
                    counter += 1
            dict_counters = {
                2: 0.75,  # ab an bn bn n != b
                1: 0.5,  # ab an bn any != bn
                0: 0  # ab an any != bn any != bn
            }
            for key in dict_counters.keys():
                if counter == key:
                    return dict_counters[key]
