from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException

# FORMULA_TYPE_GRANDPARENT
class GrandParentFormula(Formula):
    def calculate_relation(self, raw_values):
        if len(raw_values) < 3:
            # skip line with warning
            raise LineFormatException()

        gc_alleles = raw_values.pop()
        gp_alleles = raw_values.pop()
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

    # Probability of relation theory refutation in case of gc's homozygosity
    def homo_gc_refutation(self, freq):
        return (self.freq_pat1(freq)) ** 2

    # Probability of relation theory refutation in case of gc's heterozygosity
    def hetero_gc_refutation(self, freq1, freq2):
        return 2 * self.freq_pat1(freq1) * self.freq_pat1(freq2) - (2 * freq1 * freq2)**2

    # Probability of relation theory confirmation in case of gc's homozygosity
    def homo_gc_confirmation(self, freq, intersection, gp_set):
        if len(intersection) == 0:
            confirmation = freq * self.freq_pat2(freq)
        else:
            if len(gp_set) == 1:
                confirmation = self.freq_pat1(freq)
            else:
                confirmation = self.freq_pat1(freq) * self.freq_pat2(freq)

        return confirmation

    # Probability of relation theory confirmation in case of gc's heterozygosity
    def hetero_gc_confirmation(self, freq1, freq2, intersection, gp_set):
        if len(intersection) == 0:
            confirmation = freq1 * self.freq_pat1(freq2) + freq2 * self.freq_pat1(freq1)
        elif len(intersection) == 2:
            confirmation = self.freq_pat2(freq1) * self.freq_pat1(freq2) + self.freq_pat2(freq2) * self.freq_pat1(freq1) - freq1 * freq2 * (freq1 + freq2)
        else:
            if len(gp_set) == 2:
                confirmation = self.freq_pat1(freq2) + freq2 * (self.freq_pat1(freq1) - 2 * freq1 * freq2)
            else:
                confirmation = self.freq_pat2(freq1) * self.freq_pat1(freq2) + freq2 * self.freq_pat1(freq1)  - freq1 * freq2 **2

        return confirmation







"""
class GrandParentFormula(Formula):
    def calculate_relation(self, row_values):
        if len(row_values) < 3:
            # skip line with warning
            raise LineFormatException()

        raw_ab = row_values.pop()
        raw_cd = row_values.pop()
        locus = ' '.join(row_values)

        ab = self.split_sat(raw_ab)
        cd = self.split_sat(raw_cd)

        if len(ab) != 2 or len(cd) != 2:
            # skip line with warning
            raise AllelesException()

        b = ab.pop()
        a = ab.pop()
        d = cd.pop()
        c = cd.pop()

        lr = 0
        # Child AA
        if a == b:
            # Grand parent AA
            if c == a and d == a:
                lr = self.gen_aa_aa(locus, a)
            # Grand parent AB
            elif c == a or d == a:
                lr = self.gen_aa_ab(locus, a)
            # Grand parent BB or BC
            else:
                lr = self.gen_aa_bc(locus, a)
        # Child AB
        elif c == d:
            # GrandParent AA
            if a == c:
                lr = self.gen_ab_aa(locus, c, b)
            # GrandParent AA
            elif b == c:
                lr = self.gen_ab_aa(locus, c, a)
            # GrandParent CC
            else:
                lr = self.gen_ab_cd(locus, a, b)
        # GrandParent AB
        elif (a == c and b == d) or (a == d and b == c):
            lr = self.gen_ab_ab(locus, a, b)
        else:
            if a == c or a == d:
                lr = self.gen_ab_ac(locus, a, b)
            elif b == c or b == d:
                lr = self.gen_ab_ac(locus, b, a)
            else:
                lr = self.gen_ab_cd(locus, a, b)

        return self.make_result(locus, raw_ab, raw_cd, lr)

    # Child AA Formulas
    def gen_aa_aa(self, locus, a):
        p = self.get_frequencies(locus, {a: 0})
        return 0 if p[a] == 0 else self._2pa_sub_pa2(p, a) / self.prob_not_c_aa(p, a)

    def gen_aa_ab(self, locus, a):
        p = self.get_frequencies(locus, {a: 0})
        return 0 if p[a] == 0 else self._05_add_05pa(p, a) * self._2pa_sub_pa2(p, a) / self.prob_not_c_aa(p, a)

    def gen_aa_bc(self, locus, a):
        p = self.get_frequencies(locus, {a: 0})
        return 0 if p[a] == 0 else p[a] * self._2pa_sub_pa2(p, a) / self.prob_not_c_aa(p, a)

    # Child AB Formulas
    def gen_ab_aa(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            (self._2pa_sub_pa2(p, b) + p[b] * (self._2pa_sub_pa2(p, a) - self._2_pa_pb(p, a, b))) / divider

    def gen_ab_ab(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            (
                    self._05_add_05pa(p, a) * self._2pa_sub_pa2(p, b)
                    + self._05_add_05pa(p, b) * self._2pa_sub_pa2(p, a)
                    - .5 * (p[a] + p[b]) * self._2_pa_pb(p, a, b)
            ) / divider

    def gen_ab_ac(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            (
                    self._05_add_05pa(p, a) * self._2pa_sub_pa2(p, b)
                    + p[b] * self._2pa_sub_pa2(p, a)
                    - .5 * p[b] * self._2_pa_pb(p, a, b)
            ) / divider

    def gen_ab_cd(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            (p[a] * self._2pa_sub_pa2(p, b) + p[b] * self._2pa_sub_pa2(p, a)) / divider
    
"""
