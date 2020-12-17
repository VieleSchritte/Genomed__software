from __future__ import unicode_literals
from .base import Formula, Calculations
from .one_known_supposed import OneKnownSupposedFormula
from .two_known_supposed import TwoKnownSupposedFormula


class ThreeKnownSupposed(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 5)
        known_alleles, supposed_alleles, children_genotypes = alleles[0], alleles[1], alleles[2:]
        known_set, supposed_set, children_sets = sets[0], sets[1], sets[2:]

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        c = Calculations()
        children_alleles = c.get_overall_alleles(children_genotypes)
        freq_dict = self.get_frequencies(locus, children_alleles)
        len_counter, children_inters_lens = 0, []
        for child_set in children_sets:
            if len(child_set) == 2:
                len_counter += 1
        for i in range(7, 10):
            children_inters_lens.append(len(intersections[i]))
        if len_counter < 3 and children_inters_lens.count(0) != 1:  # intersections = 0 as usual
            for i in range(4, 7):
                if len(intersections[i]) == 0:
                    return self.make_result(locus, 0, dict_make_result)
        else:  # the last heterozygous case
            freq1, freq2 = freq_dict[supposed_alleles[0]], freq_dict[supposed_alleles[1]]
            lr = 2 * freq1 * freq2
            return self.make_result(locus, 1 / lr, dict_make_result)

        raw_values = [locus, '/'.join(known_alleles), '/'.join(supposed_alleles)]
        formulas = [OneKnownSupposedFormula(Formula), TwoKnownSupposedFormula(Formula)]
        lr = c.get_repeatable_lr(raw_values, children_genotypes, formulas)
        if lr:
            return self.make_result(locus, lr, dict_make_result)
        else:
            possible_parent_genotypes = c.get_possible_genotypes(children_alleles, children_genotypes, [known_set, 'known'])
            lr = c.get_lr_from_possible(possible_parent_genotypes, freq_dict)
            return self.make_result(locus, 1 / lr, dict_make_result)
