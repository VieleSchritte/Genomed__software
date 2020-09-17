from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException

class UncleFormula(Formula):
    def calculate_relation(self, raw_values):
        if len(raw_values) < 3:
            raise LineFormatException()

        nephew_alleles = self.split_sat(raw_values.pop())
        uncle_alleles = self.split_sat(raw_values.pop())
        locus = ''.join(raw_values)

        if len(nephew_alleles) != 2 or len(uncle_alleles) != 2:
            raise AllelesException()

