from __future__ import unicode_literals
from .base import Formula, AllelesException


# FORMULA_TYPE_GRANDPARENT
class GrandParentFormula(Formula):
    def calculate_relation(self, raw_values):
        (gc_alleles, gp_alleles, locus, gc_set, gp_set, intersection) = self.getting_alleles_locus(raw_values)

        # Checking gender specificity of locus
        if self.is_gender_specific(locus):
            return self.make_result(locus, '/'.join(gc_alleles), '/'.join(gp_alleles), '-')

        if len(gc_alleles) != 2 or len(gp_alleles) != 2:
            raise AllelesException()

        freq_dict = self.get_frequencies(locus, gc_set)
        calc = Calculations()

        if len(gc_set) == 1:
            freq = freq_dict[next(iter(gc_set))]  # gets first
            refutation = calc.homo_gc_refutation(freq)
            confirmation = calc.homo_gc_confirmation(freq, intersection, gp_set)
        else:  # gc_set == 2
            gc1 = gc_alleles[0]
            gc2 = gc_alleles[1]

            if gc2 in gp_set:
                gc1, gc2 = gc2, gc1
            freq1 = freq_dict[gc1]
            freq2 = freq_dict[gc2]

            refutation = calc.hetero_gc_refutation(freq1, freq2)
            confirmation = calc.hetero_gc_confirmation(freq1, freq2, intersection, gp_set)
        lr = confirmation / refutation

        return self.make_result(locus, '/'.join(gc_alleles), '/'.join(gp_alleles), lr)


class Calculations:
    # A helper for the frequently used pattern F(Px) = Px * (2 - Px)
    @staticmethod
    def F(freq):
        return freq * (2 - freq)

    # A helper for the frequently used pattern Q(Px) = 0.5 - 0.5 * Px
    @staticmethod
    def Q(freq):
        return 0.5 + 0.5 * freq

    # Probability of relation theory refutation in case of grandchild's homozygosity
    def homo_gc_refutation(self, freq):
        return (self.F(freq)) ** 2

    # Probability of relation theory refutation in case of grandchild's heterozygosity
    def hetero_gc_refutation(self, freq1, freq2):
        return 2 * self.F(freq1) * self.F(freq2) - (2 * freq1 * freq2) ** 2

    # Probability of relation theory confirmation in case of grandchild's homozygosity
    def homo_gc_confirmation(self, freq, intersection, gp_set):
        if len(intersection) == 0:
            return freq * self.F(freq)

        if len(gp_set) == 1:
            return self.F(freq)

        return self.F(freq) * self.Q(freq)

    # Probability of relation theory confirmation in case of grandchild's heterozygosity
    def hetero_gc_confirmation(self, freq1, freq2, intersection, gp_set):
        if len(intersection) == 0:
            # case ab nn (nk)
            return freq1 * self.F(freq2) + freq2 * self.F(freq1)

        if len(intersection) == 2:
            # case ab ab
            return self.Q(freq1) * self.F(freq2) + self.Q(freq2) * self.F(freq1) - freq1 * freq2 * (freq1 + freq2)

        if len(gp_set) == 2:
            # case ab an
            return self.Q(freq1) * self.F(freq2) + freq2 * self.F(freq1) - freq1 * freq2 ** 2

        # default is ab aa case
        return self.F(freq2) + freq2 * (self.F(freq1) - 2 * freq1 * freq2)
