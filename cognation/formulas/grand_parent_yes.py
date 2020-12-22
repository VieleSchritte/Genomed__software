from __future__ import unicode_literals
from cognation.formulas.base import Formula, Calculations


class GrandParentYes(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 3)
        child_alleles, grandparent_alleles, parent_alleles = alleles
        child_set, grandparent_set, parent_set, = sets

        if self.is_gender_specific(locus):
            return self.preparation_check(locus, dict_make_result)

        c = Calculations()
        all_alleles = c.all
        freq_dict = self.get_frequencies(locus, list(common_set))

        if len(child_set) == 1:  # Homozygous child
            if len(intersections[0]) != 0:
                freq = freq_dict[child_alleles[0]]
                lr = c.
