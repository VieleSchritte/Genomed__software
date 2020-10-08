from __future__ import unicode_literals
from cognation.formulas.for_3_participants.base import Formula, LineFormatException, AllelesException

class SiblingFormula(Formula):
    def calculate_relation(self, raw_values):
        if len(raw_values) < 4:
            raise LineFormatException()

        raw_ef = raw_values.pop()
        raw_cd = raw_values.pop()
        raw_ab = raw_values.pop()
        locus = ' '.join(raw_values)

        split_ab = self.split_sat(raw_ab)
        split_cd = self.split_sat(raw_cd)
        split_ef = self.split_sat(raw_ef)

        if len(split_ab) != 2 or len(split_cd) != 2 or len(split_ef) != 2:
            raise AllelesException()

        ab = set(split_ab)
        cd = set(split_cd)
        ef = set(split_ef)

        ab_cd = ab & cd
        ab_ef = ab & ef
        cd_ef = cd & ef

        lr = 0
        if len(ab) == 1:  # 1 - AA
            if len(ab_cd) == 1:
                if len(cd) == 1:  # 2 - AA
                    if len(cd_ef) == 1:
                        if len(ef) == 1:  # 3 - AA
                            lr = self.gen_aa_aa_aa(locus, ab.copy().pop())
                        else:  # 3 - AB
                            lr = self.gen_aa_aa_ab(locus, ab.copy().pop(), (ab ^ ef).pop())
                    else:  # 3 - BC or CC
                        pass
                elif len(cd) == 2:  # 2 - AB
                    if len(ab_ef) == 1:
                        if len(ef) == 1:  # 3 - AA
                            lr = self.gen_aa_ab_aa(locus, ab.copy().pop())
                        elif len(cd_ef) == 2:  # 3 - AB
                            lr = self.gen_aa_ab_ab(locus, ab.copy().pop(), (ab ^ cd).pop())
                        elif len(cd_ef) == 1:  # 3 - AC
                            lr = self.gen_aa_ab_ac(locus, ab.copy().pop(), (ab ^ ef).pop())
                    else:  # 3 - CC or CD
                        pass
            else:  # 2 - BC
                pass
        # 1 - AB
        elif len(ab_cd) == 1:
            if len(cd) == 1:  # 2 - AA
                a = cd.copy().pop()
                b = (cd ^ ab).pop()
                if len(ef) == 1 and len(cd_ef) == 1: # 3 - AA
                    lr = self.gen_ab_aa_aa(locus, a, b)
                elif len(ef) == 2:
                    if len(ab_ef) == 2:  # 3 - AB
                        lr = self.gen_ab_aa_ab(locus, a, b)
                    elif len(ab_ef) == 1:  # 3 - AC
                        lr = self.gen_ab_aa_ac(locus, a, b, (cd ^ ef).pop())
                    else:  # 3 - BC
                        pass
                else:  # 3 - CC
                    pass
            else:  # 2 - AC
                a = ab_cd.copy().pop()
                b = (ab_cd ^ ab).pop()
                c = (ab_cd ^ cd).pop()
                if len(ef) == 1:
                    if a in ef:  # 3 - AA
                        lr = self.gen_ab_ac_aa(locus, a, b)
                    elif c in ef:  # 3 - CC
                        lr = self.gen_ab_ac_cc(locus, a, b, c)
                else:
                    if len(ab_ef) == 2:  # 3 - AB
                        lr = self.gen_ab_ac_ab(locus, a, b)
                    elif len(cd_ef) == 2:  # 3 - AC
                        lr = self.gen_ab_ac_ac(locus, a, b, c)
                    elif len(ab_ef) == 1 and a in ef:
                        lr = self.gen_ab_ac_ad(locus, a, b, (ab_ef ^ ef).pop())
                    elif b in ef and c in ef:  # 3 - BC
                        lr = self.gen_ab_ac_bc(locus, a, b)
                    elif len(ab_ef) == 0 and len(cd_ef) == 1:
                        lr = self.gen_ab_ac_cd(locus, a, b, (cd_ef ^ ef).pop())
        elif len(ab_cd) == 2:  # 2 - AB
            if len(ef) == 1 and len(ab_ef) == 1:  # 3 - AA
                lr = self.gen_ab_ab_aa(locus, ab_ef.pop(), (ab ^ ef).pop())
            elif len(ef) == 2 and len(ab_ef) == 2:  # 3 - AB
                lr = self.gen_ab_ab_ab(locus, ab.pop(), ab.pop())
            elif len(ef) == 2 and len(ab_ef) == 1:  # 3 - AC
                lr = self.gen_ab_ab_ac(locus, ab_ef.copy().pop(), (ab_ef ^ cd).pop(), (ab_ef ^ ef).pop())

        result = self.make_result(locus, raw_ab, raw_cd, lr)
        result["ef"] = raw_ef
        return result

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

