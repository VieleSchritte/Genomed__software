from __future__ import unicode_literals
from .base import Formula, Calculations
from .parent import ParentFormula


class TwoChildrenFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, part_alleles, part_sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        parent_alleles, children_genotypes = part_alleles[0], part_alleles[1:]
        parent_set, children_sets = part_sets[0], part_sets[1:]
        ch1ch2_inter = intersections[-1]

        if self.is_gender_specific(locus):
            return self.result_gender_specific(locus, dict_make_result)
        for i in range(2):
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)

        raw_values = [locus, '/'.join(parent_alleles)]
        c = Calculations()
        lr = c.get_repeatable_lr(raw_values, children_genotypes, [ParentFormula(Formula)])
        if lr:
            return self.make_result(locus, lr, dict_make_result)
        else:
            children_alleles = c.get_overall_alleles(children_genotypes)
            freq_dict = self.get_frequencies(locus, children_alleles)
            counter = 0
            if len(ch1ch2_inter) == 1:
                for child_set in children_sets:
                    if len(child_set) == 2:
                        counter += 1
            if counter == 2:  # ab ac an/bc
                lr = 2
                for child_set in children_sets:
                    lr *= freq_dict[list(child_set - ch1ch2_inter)[0]]
                lr += c.F(freq_dict[list(ch1ch2_inter)[0]])
                return self.make_result(locus, 1 / lr, dict_make_result)

            possible_parents_genotypes = c.get_possible_genotypes(children_alleles, children_genotypes, [parent_set, 'supposed'])
            lr = c.get_lr_from_possible(possible_parents_genotypes, freq_dict)
            return self.make_result(locus, 1 / lr, dict_make_result)
