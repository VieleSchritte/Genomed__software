from __future__ import unicode_literals
from .base import Formula, AllelesException
from .base import Calculations


# FORMULA_TYPE_GRANDPARENT
class GrandParentFormula(Formula):
    def calculate_relation(self, raw_values):
        (gc_alleles, gp_alleles, locus, gc_set, gp_set, intersection) = self.getting_alleles_locus(raw_values, 2)

        # Checking gender specificity of locus
        if self.is_gender_specific(locus):
            return self.make_result2(locus, '/'.join(gc_alleles), '/'.join(gp_alleles), '-')

        if len(gc_alleles) != 2 or len(gp_alleles) != 2:
            raise AllelesException()

        freq_dict = self.get_frequencies(locus, gc_set)
        calc = Calculations()
        conf = Confirmations()

        if len(gc_set) == 1:
            freq = freq_dict[next(iter(gc_set))]  # gets first
            refutation = calc.homo_refutation(freq)
            confirmation = conf.homo_gc_confirmation(freq, intersection, gp_set)
        else:  # gc_set == 2
            gc1 = gc_alleles[0]
            gc2 = gc_alleles[1]

            if gc2 in gp_set:
                gc1, gc2 = gc2, gc1
            freq1 = freq_dict[gc1]
            freq2 = freq_dict[gc2]

            refutation = calc.hetero_refutation(freq1, freq2)
            confirmation = conf.hetero_gc_confirmation(freq1, freq2, intersection, gp_set)
        lr = confirmation / refutation

        return self.make_result2(locus, '/'.join(gc_alleles), '/'.join(gp_alleles), lr)


class Confirmations:
    # Probability of relation theory confirmation in case of grandchild's homozygosity
    @staticmethod
    def homo_gc_confirmation(freq, intersection, gp_set):
        calc = Calculations()
        if len(intersection) == 0:
            return freq * calc.F(freq)

        if len(gp_set) == 1:
            return calc.F(freq)

        return calc.F(freq) * calc.Q(freq)

    # Probability of relation theory confirmation in case of grandchild's heterozygosity
    @staticmethod
    def hetero_gc_confirmation(freq1, freq2, intersection, gp_set):
        calc = Calculations()
        if len(intersection) == 0:
            # case ab nn (nk)
            return freq1 * calc.F(freq2) + freq2 * calc.F(freq1)

        if len(intersection) == 2:
            # case ab ab
            return calc.Q(freq1) * calc.F(freq2) + calc.Q(freq2) * calc.F(freq1) - freq1 * freq2 * (freq1 + freq2)

        if len(gp_set) == 2:
            # case ab an
            return calc.Q(freq1) * calc.F(freq2) + freq2 * calc.F(freq1) - freq1 * freq2 ** 2

        # default is ab aa case
        return calc.F(freq2) + freq2 * (calc.F(freq1) - 2 * freq1 * freq2)
