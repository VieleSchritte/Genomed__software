from __future__ import unicode_literals
from cognation.formulas.base import Formula, AllelesException, Calculations


class SiblingFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections = self.getting_alleles_locus(raw_values, 3)

        sibling_alleles, parent_alleles, child_alleles = alleles
        sibling_set, parent_set, child_set = sets
        sp_intersection, sc_intersection, cp_intersection = intersections

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr
        if self.is_gender_specific(locus):
            return self.make_result3(locus, '/'.join(child_alleles), '/'.join(parent_alleles), '/'.join(sibling_alleles), '-')

        # If it's not a gender specific locus there should be 2 alleles
        if len(child_alleles) != 2 or len(parent_alleles) != 2 or len(sibling_alleles) != 2:
            raise AllelesException()

        #  If there's no relation then return lr = 0 and start collecting mutations
        if sp_intersection == 0 or cp_intersection == 0:
            return self.make_result3(locus, '/'.join(child_alleles), '/'.join(parent_alleles), '/'.join(sibling_alleles), 0)

        # If everything's ok and locus isn't gender-specific, start calculating
        c = Calculations()
        freq_dict = self.get_frequencies(locus, child_alleles)
        parent2_alleles = [0, 0]

        if len(child_set) == 1:
            freq = freq_dict[child_alleles[0]]
            refutation = c.homo_refutation(freq)

        else:
            freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
            refutation = c.hetero_refutation(freq1, freq2)

        parent2_alleles = self.get_parent2_alleles(parent2_alleles, sibling_alleles, sp_intersection, 0)
        parent2_alleles = self.get_parent2_alleles(parent2_alleles, child_alleles, cp_intersection, 1)

        freq_dict = self.get_frequencies(locus, parent2_alleles)
        freq1, freq2 = freq_dict[parent2_alleles[0]], freq_dict[parent2_alleles[1]]
        confirmation = c.M(freq1, freq2)

        lr = confirmation / refutation

        return self.make_result3(locus, '/'.join(child_alleles), '/'.join(parent_alleles), '/'.join(sibling_alleles), lr)

    # A method to fill M(freq1, freq2) formula in base.Calculations. Returns list of alleles in required range
    @staticmethod
    def get_parent2_alleles(parent2_alleles, ch_sib_alleles, intersection, position):
        counter = 0

        for allele in ch_sib_alleles:
            if allele == intersection:
                counter += 1

        if counter == 2:
            parent2_alleles[position] = ch_sib_alleles[position]
        else:
            for allele in ch_sib_alleles:
                if allele == intersection:
                    continue
                else:
                    parent2_alleles[position] = allele

        return parent2_alleles
