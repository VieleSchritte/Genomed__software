from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException


# FORMULA_TYPE_GRANDPARENT
class GrandParentFormula(Formula):
    def calculate_relation(self, raw_values):
        if len(raw_values) < 3:
            # skip line with warning
            raise LineFormatException()

        gc_alleles = self.split_sat(raw_values.pop())
        gp_alleles = self.split_sat(raw_values.pop())
        locus = ' '.join(raw_values)

        gc_set = list(set(gc_alleles))
        gp_set = list(set(gp_alleles))
        intersection = list(set(gc_alleles) & set(gp_alleles))
        freq_dict = self.get_frequencies(locus, gc_set)

        inst = Calculations()

        if len(gc_set) == 1:
            freq = freq_dict[gc_set[0]]
            refutation = inst.homo_gc_refutation(freq)
            confirmation = inst.homo_gc_confirmation(freq, intersection, gp_set)
        else:
            freq1 = freq_dict[gc_set[0]]
            freq2 = freq_dict[gc_set[1]]
            refutation = inst.hetero_gc_refutation(freq1, freq2)
            confirmation = inst.hetero_gc_confirmation(freq1, freq2, intersection, gp_set)
        lr = confirmation / refutation

        return self.make_result(locus, '/'.join(gp_alleles), '/'.join(gc_alleles), lr)


class Calculations:
    # A helper for the frequently used pattern F(Px) = Px * (2 - Px)
    @staticmethod
    def freq_pat1(freq):
        return freq * (2 - freq)

    # A helper for the frequently used pattern Q(Px) = 0.5 - 0.5 * Px
    @staticmethod
    def freq_pat2(freq):
        return 0.5 - 0.5 * freq

    # Probability of relation theory refutation in case of grandchild's homozygosity
    def homo_gc_refutation(self, freq):
        return (self.freq_pat1(freq)) ** 2

    # Probability of relation theory refutation in case of grandchild's heterozygosity
    def hetero_gc_refutation(self, freq1, freq2):
        return 2 * self.freq_pat1(freq1) * self.freq_pat1(freq2) - (2 * freq1 * freq2)**2

    # Probability of relation theory confirmation in case of grandchild's homozygosity
    def homo_gc_confirmation(self, freq, intersection, gp_set):
        if len(intersection) == 0:
            confirmation = freq * self.freq_pat2(freq)
        else:
            if len(gp_set) == 1:
                confirmation = self.freq_pat1(freq)
            else:
                confirmation = self.freq_pat1(freq) * self.freq_pat2(freq)

        return confirmation

    # Probability of relation theory confirmation in case of grandchild's heterozygosity
    def hetero_gc_confirmation(self, freq1, freq2, intersection, gp_set):
        if len(intersection) == 0:
            confirmation = freq1 * self.freq_pat1(freq2) + freq2 * self.freq_pat1(freq1)
        elif len(intersection) == 2:
            confirmation = self.freq_pat2(freq1) * self.freq_pat1(freq2) + self.freq_pat2(freq2) * self.freq_pat1(freq1) - freq1 * freq2 * (freq1 + freq2)
        else:
            if len(gp_set) == 2:
                confirmation = self.freq_pat1(freq2) + freq2 * (self.freq_pat1(freq1) - 2 * freq1 * freq2)
            else:
                confirmation = self.freq_pat2(freq1) * self.freq_pat1(freq2) + freq2 * self.freq_pat1(freq1) - freq1 * freq2 ** 2

        return confirmation
