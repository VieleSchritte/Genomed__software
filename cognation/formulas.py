from __future__ import unicode_literals

from .models import Locus

import abc
import re
from collections import Counter, OrderedDict

FORMULA_TYPE_PARENT = 1
FORMULA_TYPE_GRANDPARENT = 2
FORMULA_TYPE_SIBLING = 3

FORMULA_TYPE_UNCLE = 4
FORMULA_TYPE_COUSIN = 5
FORMULA_TYPE_STEPBROTHER = 6

#Tells pyhton what to do in every formula case (FORMULA_...)
def formula_builder(type, data):
    normalized_type = int(type)
    if normalized_type == FORMULA_TYPE_PARENT:
        return ParentFormula(data)
    elif normalized_type == FORMULA_TYPE_GRANDPARENT:
        return GrandParentFormula(data)
    elif normalized_type == FORMULA_TYPE_SIBLING:
        return SiblingFormula(data)
    #elif normalized_type == FORMULA_TYPE_UNCLE:
    #    return SiblingFormula(data)
    #elif normalized_type == FORMULA_TYPE_COUSIN:
    #    return SiblingFormula(data)
    #elif normalized_type == FORMULA_TYPE_STEPBROTHER:
    #    return SiblingFormula(data)
    #else:
    #    raise UnknownFormulaException(type)


# Exception - if something isn't right in data format
class UnknownFormulaException(Exception):
    def __init__(self, formula_type):
        self.formula_type = formula_type


class LineFormatException(Exception):
    def __str__(self):
        return 'Wrong line format'


class AllelesException(Exception):
    def __str__(self):
        return "Alleles count doesn't look right"


class UnknownAlleleException(Exception):
    def __init__(self, locus, sat):
        self.locus = locus
        self.sat = sat

    def __str__(self):
        return "Unknown allele found: " + str(self.sat)


# Abstract parent class
class Formula(abc.ABC):
    def __init__(self, user_data):
        self.user_data = str(user_data)

    def calculate(self):
        result = OrderedDict()
        lines = self.user_data.splitlines()
        for line in lines:
            line = line.strip(' \t\n\r')
            if len(line) == 0:
                continue

            try:
                relation = self.calculate_relation(re.split(r'[\s\t]+', line))
                result[relation['locus']] = relation
            except (LineFormatException, AllelesException, UnknownAlleleException) as exception:
                result[hash(line)] = {'exception': exception, 'line': line}

        return result

    # getting allele frequencies from DB
    def get_frequencies(self, locus, sat_list):
        result = {}
        for sat in sat_list:
            try:
                locus_object = Locus.objects.get(locus=locus, sat=self.normalize_sat(sat))
                result[sat] = locus_object.freq
            except Locus.DoesNotExist:
                raise UnknownAlleleException(locus, sat)

        return result

    def get_template(self):
        return 'cognation/formula/' + self.__class__.__name__.lower()[:-7] + '.html'

    @staticmethod
    def normalize_sat(value):
        if value == 'X':
            return 0.0
        elif value == 'Y':
            return 1.0
        else:
            return float(value)

    @staticmethod
    def split_sat(sat_string):
        return re.split(r'\/', sat_string)

    @staticmethod
    def get_sat_counter(sat_string):
        return Counter(Formula.split_sat(sat_string))

    @staticmethod
    def make_result(locus, ab, cd, lr):
        return {
            "locus": locus,
            "ab": ab,
            "cd": cd,
            "lr": lr
        }

    # Calculation Helpers
    @staticmethod
    def _2pa_sub_pa2(p, x):
        return p[x] * (2 - p[x])

    @staticmethod
    def _05_add_05pa(p, x):
        return .5 * (1 + p[x])

    @staticmethod
    def _2_pa_pb(p, a, b):
        return 2 * p[a] * p[b]

    def prob_not_c_aa(self, p, a):
        return self._2pa_sub_pa2(p, a) ** 2

    def prob_not_c_ab(self, p, a, b):
        return 2 * self._2pa_sub_pa2(p, a) * self._2pa_sub_pa2(p, b) \
               - self._2_pa_pb(p, a, b) * self._2_pa_pb(p, a, b)

    # Abstract methods
    @abc.abstractmethod
    def calculate_relation(self, row_values):
        """ Abstract method to calculate """


# FORMULA_TYPE_PARENT
class ParentFormula(Formula):
    def calculate_relation(self, row_values):
        if len(row_values) < 3:
            raise LineFormatException()

        child_alleles = self.split_sat(row_values.pop())
        parent_alleles = self.split_sat(row_values.pop())
        locus = ' '.join(row_values)  # for loci names contain space

        if len(parent_alleles) != 2 or len(child_alleles) != 2:
            raise AllelesException()

        parent_set = set(parent_alleles)  # unique parent alleles
        child_set = set(child_alleles)  # unique child alleles

        relation = 0
        freq_dict = self.get_frequencies(locus, list(parent_set & child_set))
        # i, j     i, j
        if len(freq_dict) == 2:
            divider = 4
            dividend = 0
            for key in freq_dict:
                divider = divider * freq_dict[key]
                dividend = dividend + freq_dict[key]
            # (f1 + f2) / 4 * f1 * f2

            if divider != 0:
                relation = dividend / divider
        elif len(freq_dict) == 1:
            # i, j    i, n
            # i, i    i, n
            # i, j    j, n
            # 1 / {1,2}*{1,2}*f
            freq = next(iter(freq_dict.values()))
            relation = 1 / (len(parent_set)*len(child_set) * freq)

        return self.make_result(locus, '/'.join(parent_alleles), '/'.join(child_alleles), relation)


# FORMULA_TYPE_GRANDPARENT
class GrandParentFormula(Formula):
    def calculate_relation(self, row_values):
        if len(row_values) < 3:
            # skip line with warning
            raise LineFormatException()

        raw_cd = row_values.pop()
        raw_ab = row_values.pop()
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


# FORMULA_TYPE_SIBLING
class SiblingFormula(Formula):
    def calculate_relation(self, row_values):
        if len(row_values) < 4:
            raise LineFormatException()

        raw_ef = row_values.pop()
        raw_cd = row_values.pop()
        raw_ab = row_values.pop()
        locus = ' '.join(row_values)

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

#FORMULA_TYPE_UNCLE
class UncleFormula(Formula):
    def calculate_relation(self, row_values):
        if len(row_values) < 3:
            raise LineFormatException()




#FORMULA_TYPE_COUSIN
#class CousinFormula(Formula):


#FORMULA_TYPE_STEPBROTHER
#class StepbrotherFormula(Formula):
