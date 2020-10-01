from __future__ import unicode_literals
from .base import Formula, AllelesException


# FORMULA_TYPE_GRANDPARENT
class GrandParentFormula(Formula):
    def calculate_relation(self, raw_values):
        (gc_alleles, gp_alleles, locus, gc_set, gp_set, intersection) = self.getting_alleles_locus(raw_values)

        # Checking gender specificity of locus
        lr = self.gender_specific(locus)
        
        if lr is None:
            if len(gc_alleles) != 2 or len(gp_alleles) != 2:
                raise AllelesException()

            freq_dict = self.get_frequencies(locus, gc_set)
            call_calc = Calculations()

            if len(gc_set) == 1:
                freq = freq_dict[gc_set[0]]
                refutation = call_calc.homo_gc_refutation(freq)
                confirmation = call_calc.homo_gc_confirmation(freq, intersection, gp_set)
            else:
                freq1 = freq_dict[gc_set[0]]
                freq2 = freq_dict[gc_set[1]]
                refutation = call_calc.hetero_gc_refutation(freq1, freq2)
                confirmation = call_calc.hetero_gc_confirmation(freq1, freq2, intersection, gp_set)
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
        return 0.5 - 0.5 * freq

    # Probability of relation theory refutation in case of grandchild's homozygosity
    def homo_gc_refutation(self, freq):
        return (self.F(freq)) ** 2

    # Probability of relation theory refutation in case of grandchild's heterozygosity
    def hetero_gc_refutation(self, freq1, freq2):
        return 2 * self.F(freq1) * self.F(freq2) - (2 * freq1 * freq2) ** 2

    # Probability of relation theory confirmation in case of grandchild's homozygosity
    def homo_gc_confirmation(self, freq, intersection, gp_set):
        if len(intersection) == 0:
            confirmation = freq * self.Q(freq)
        else:
            if len(gp_set) == 1:
                confirmation = self.F(freq)
            else:
                confirmation = self.F(freq) * self.Q(freq)

        return confirmation

    # Probability of relation theory confirmation in case of grandchild's heterozygosity
    def hetero_gc_confirmation(self, freq1, freq2, intersection, gp_set):
        if len(intersection) == 0:
            confirmation = freq1 * self.F(freq2) + freq2 * self.F(freq1)
        elif len(intersection) == 2:
            confirmation = self.Q(freq1) * self.F(freq2) + self.Q(freq2) * self.F(freq1) - freq1 * freq2 * (freq1 + freq2)
        else:
            if len(gp_set) == 2:
                confirmation = self.F(freq2) + freq2 * (self.F(freq1) - 2 * freq1 * freq2)
            else:
                confirmation = self.Q(freq1) * self.F(freq2) + freq2 * self.F(freq1) - freq1 * freq2 ** 2

        return confirmation
