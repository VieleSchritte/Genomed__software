# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.test import TestCase
from ...formulas.parent import ParentFormula


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
        cpi = 0
        p = 0.0
        with open('./test_cases/parent_cases/test_data_parent1.txt', 'r') as parent1_data:
            for line in parent1_data:
                line_list = line.strip().split('\t')
                if len(line) != 0:
                    locus = line[0]
                    lr = line[3]
                    test_dict[locus] = lr
                else:
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    elif line[1] == 'P':
                        p = float(line[1])