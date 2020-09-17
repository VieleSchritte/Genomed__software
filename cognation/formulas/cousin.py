from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException

class CousinFormula(Formula):
    def calculate_relation(self, raw_values):
        if len(raw_values) < 3:
            raise LineFormatException()

        cousin2_alleles = self.split_sat(raw_values.pop())
        cousin1_alleles = self.split_sat(raw_values.pop())

