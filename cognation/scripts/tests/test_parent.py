# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.test import TestCase
from ...formulas.parent import ParentFormula
from ...formulas.base import Formula


class TestParentFormula(TestCase):

    @classmethod
    def setUpTestData(cls):
        print("setUpTestData: Run once to set up non-modified data for all class methods.")
        pass

    def setUp(self):
        print("setUp: Run once for every test method to setup clean data.")
        pass

    def testParentFormula(self):
        test_dict = {}
        self.cpi = 0
        self.p = 0.0
        with open('cognation/scripts/tests/test_cases/parent_cases/test_data_parent1.txt', 'r') as parent1_data:
            for line in parent1_data:
                line = line.strip().split('\t')
                if len(line) != 1:
                    locus = line.pop(0)
                    lr = line.pop()
                    test_dict[locus] = lr
                else:
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        self.cpi = int(line[1])
                    elif line[0] == 'P':
                        self.p = float(line[1])

        result_string = ParentFormula(Formula)
        print(result_string)