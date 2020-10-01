from __future__ import unicode_literals

from cognation.models import Locus

import abc
import re
from collections import Counter, OrderedDict


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

    # Checking out if the locus is gender-specific (so we don't need to add it to cpi calculation)
    @staticmethod
    def is_gender_specific(locus):
        gender_specific_loci = ['SRY', 'DYS391', 'Yindel']
        for i in range(len(gender_specific_loci)):
            if locus == gender_specific_loci[i]:
                return True
        return False

    def getting_alleles_locus(self, raw_values):
        if len(raw_values) < 3:
            # Skip line with warning
            raise LineFormatException()

        # child/grandchild for example
        pat1_alleles = self.split_sat(raw_values.pop())
        # parent/grandparent...
        pat2_alleles = self.split_sat(raw_values.pop())

        locus = ' '.join(raw_values)  # for loci names contain space
        pat1_set = list(set(pat1_alleles)) # unique alleles
        pat2_set = list(set(pat2_alleles))
        intersection = list(set(pat1_alleles) & set(pat2_alleles)) # common unique alleles

        return pat1_alleles, pat2_alleles, locus, pat1_set, pat2_set, intersection

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
