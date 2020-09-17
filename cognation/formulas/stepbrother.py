from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException

class StepbrotherFormula(Formula):
    def calculate_relation(self, raw_values):
        if len(raw_values) < 3:
            raise LineFormatException()

        stbrother2_alleles = self.split_sat(raw_values.pop())
        stbrother1_alleles = self.split_sat(raw_values.pop())
