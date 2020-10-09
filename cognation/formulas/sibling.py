from __future__ import unicode_literals
from cognation.formulas.base import Formula, AllelesException


class SiblingFormula(Formula):
    def calculate_relation(self, raw_values):

        (locus, alleles_list, sets_list, inter_list) = self.getting_alleles_locus(raw_values, 3)

        child_alleles, parent_alleles, sibling_alleles = alleles_list
        child_set, parent_set, sibling_set = sets_list

        #  cp = child and parent, cs = child and sibling, ps = parent and sibling
        intersection_cp, intersection_cs, intersection_ps = inter_list

        if self.is_gender_specific(locus):
            return self.make_result(locus, '/'.join(child_alleles), '/'.join(parent_alleles), '-', '/'.join(sibling_alleles))

        if len(child_alleles) != 2 or len(parent_alleles) != 2 or len(sibling_alleles) != 2:
            raise AllelesException()

        lr = 0
        if len(child_set) == 1:  # 1 - AA
            if len(intersection_cp) == 1:
                if len(parent_set) == 1:  # 2 - AA
                    if len(intersection_ps) == 1:
                        if len(sibling_set) == 1:  # 3 - AA
                            lr = self.gen_aa_aa_aa(locus, child_set.copy().pop())
                        else:  # 3 - AB
                            lr = self.gen_aa_aa_ab(locus, child_set.copy().pop(), (child_set ^ sibling_set).pop())
                    else:  # 3 - BC or CC
                        pass
                elif len(parent_set) == 2:  # 2 - AB
                    if len(intersection_cs) == 1:
                        if len(sibling_set) == 1:  # 3 - AA
                            lr = self.gen_aa_ab_aa(locus, child_set.copy().pop())
                        elif len(intersection_ps) == 2:  # 3 - AB
                            lr = self.gen_aa_ab_ab(locus, child_set.copy().pop(), (child_set ^ parent_set).pop())
                        elif len(intersection_ps) == 1:  # 3 - AC
                            lr = self.gen_aa_ab_ac(locus, child_set.copy().pop(), (child_set ^ sibling_set).pop())
                    else:  # 3 - CC or CD
                        pass
            else:  # 2 - BC
                pass
        # 1 - AB
        elif len(intersection_cp) == 1:
            if len(parent_set) == 1:  # 2 - AA
                a = parent_set.copy().pop()
                b = (parent_set ^ child_set).pop()
                if len(sibling_set) == 1 and len(intersection_ps) == 1:  # 3 - AA
                    lr = self.gen_ab_aa_aa(locus, a, b)
                elif len(sibling_set) == 2:
                    if len(intersection_cs) == 2:  # 3 - AB
                        lr = self.gen_ab_aa_ab(locus, a, b)
                    elif len(intersection_cs) == 1:  # 3 - AC
                        lr = self.gen_ab_aa_ac(locus, a, b, (parent_set ^ sibling_set).pop())
                    else:  # 3 - BC
                        pass
                else:  # 3 - CC
                    pass
            else:  # 2 - AC
                a = intersection_cp.copy().pop()
                b = (intersection_cp ^ child_set).pop()
                c = (intersection_cp ^ parent_set).pop()
                if len(sibling_set) == 1:
                    if a in sibling_set:  # 3 - AA
                        lr = self.gen_ab_ac_aa(locus, a, b)
                    elif c in sibling_set:  # 3 - CC
                        lr = self.gen_ab_ac_cc(locus, a, b, c)
                else:
                    if len(intersection_cs) == 2:  # 3 - AB
                        lr = self.gen_ab_ac_ab(locus, a, b)
                    elif len(intersection_ps) == 2:  # 3 - AC
                        lr = self.gen_ab_ac_ac(locus, a, b, c)
                    elif len(intersection_cs) == 1 and a in sibling_set:
                        lr = self.gen_ab_ac_ad(locus, a, b, (intersection_cs ^ sibling_set).pop())
                    elif b in sibling_set and c in sibling_set:  # 3 - BC
                        lr = self.gen_ab_ac_bc(locus, a, b)
                    elif len(intersection_cs) == 0 and len(intersection_ps) == 1:
                        lr = self.gen_ab_ac_cd(locus, a, b, (intersection_ps ^ sibling_set).pop())
        elif len(intersection_cp) == 2:  # 2 - AB
            if len(sibling_set) == 1 and len(intersection_cs) == 1:  # 3 - AA
                lr = self.gen_ab_ab_aa(locus, intersection_cs.pop(), (child_set ^ sibling_set).pop())
            elif len(sibling_set) == 2 and len(intersection_cs) == 2:  # 3 - AB
                lr = self.gen_ab_ab_ab(locus, child_set.pop(), child_set.pop())
            elif len(sibling_set) == 2 and len(intersection_cs) == 1:  # 3 - AC
                lr = self.gen_ab_ab_ac(locus, intersection_cs.copy().pop(), (intersection_cs ^ parent_set).pop(), (intersection_cs ^ sibling_set).pop())

        return self.make_result(locus, '/'.join(child_alleles), '/'.join(parent_alleles), lr, '/'.join(sibling_alleles))

    def gen_aa_aa_aa(self, locus, a):
        p = self.get_frequencies(locus, {a: 0})
        return 0 if p[a] == 0 else \
            1 / self.prob_not_c_aa(p, a)

    def gen_aa_aa_ab(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        return 0 if p[a] == 0 else \
            self._2_pa_pb(p, a, b) / self._2pa_sub_pa2(p, b) / self.prob_not_c_aa(p, a)

    def gen_aa_ab_aa(self, locus, a):
        p = self.get_frequencies(locus, {a: 0})
        return 0 if p[a] == 0 else \
            1 / self.prob_not_c_aa(p, a)

    def gen_aa_ab_ab(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        return 0 if p[a] == 0 else \
            (
                self._2pa_sub_pa2(p, a)
                / ((p[a] + p[b]) * (2 - (p[a] + p[b])))
            ) / self.prob_not_c_aa(p, a)

    def gen_aa_ab_ac(self, locus, a, c):
        p = self.get_frequencies(locus, {a: 0, c: 0})
        return 0 if p[a] == 0 else \
            self._2_pa_pb(p, a, c) / self._2pa_sub_pa2(p, c) / self.prob_not_c_aa(p, a)

    def gen_ab_aa_aa(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            self._2_pa_pb(p, a, b) / self._2pa_sub_pa2(p, a) / divider

    def gen_ab_aa_ab(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            1 / divider

    def gen_ab_aa_ac(self, locus, a, b, c):
        p = self.get_frequencies(locus, {a: 0, b: 0, c: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            self._2_pa_pb(p, a, c) / self._2pa_sub_pa2(p, c) / divider

    def gen_ab_ab_aa(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            1 / divider

    def gen_ab_ab_ab(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            1 / divider

    def gen_ab_ab_ac(self, locus, a, b, c):
        p = self.get_frequencies(locus, {a: 0, b: 0, c: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            (self._2_pa_pb(p, a, c) + self._2_pa_pb(p, b, c)) / self._2pa_sub_pa2(p, c) / divider

    def gen_ab_ac_aa(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            self._2_pa_pb(p, a, b) / self._2pa_sub_pa2(p, a) / divider

    def gen_ab_ac_ab(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            1 / divider

    def gen_ab_ac_ac(self, locus, a, b, c):
        p = self.get_frequencies(locus, {a: 0, b: 0, c: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            (self._2_pa_pb(p, a, b) + self._2_pa_pb(p, b, c)) / ((p[a] + p[c]) * (2 - (p[a] + p[c]))) / divider

    def gen_ab_ac_ad(self, locus, a, b, d):
        p = self.get_frequencies(locus, {a: 0, b: 0, d: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            self._2_pa_pb(p, b, d) / self._2pa_sub_pa2(p, d) / divider

    def gen_ab_ac_cc(self, locus, a, b, c):
        p = self.get_frequencies(locus, {a: 0, b: 0, c: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            self._2_pa_pb(p, b, c) / self._2pa_sub_pa2(p, c) / divider

    def gen_ab_ac_bc(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            1 / divider

    def gen_ab_ac_cd(self, locus, a, b, d):
        p = self.get_frequencies(locus, {a: 0, b: 0, d: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            self._2_pa_pb(p, b, d) / self._2pa_sub_pa2(p, d) / divider
