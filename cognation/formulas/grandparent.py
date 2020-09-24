from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException


# FORMULA_TYPE_GRANDPARENT
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