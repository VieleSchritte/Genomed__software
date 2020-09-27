# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import django
import unittest
import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "genomed_software.settings")
from ...formulas.parent import ParentFormula
from django.test import TestCase
django.setup()

# all possible test cases
doc_names_list = ['parent1/reference_data_parent1.txt', 'parent2/reference_data_parent2.txt', 'parent3/reference_data_parent3_veri.txt']
overall_ref_dict = {}
overall_test_dict = {}

class GetParentsData(ParentFormula):
    # getting loci and lrs in the dictionary and also CPI and P from each patient's reference data
    @staticmethod
    def get_reference_data_parentx(doc_name):
        ref_dict = {}
        cpi = 0
        p = 0.0
        with open('cognation/scripts/tests/test_cases/parent_cases/' + doc_name, 'r') as parentx_data:
            for line in parentx_data:
                line = line.strip().split('\t')
                if len(line) != 1:
                    locus = line[0]
                    lr = line[3]
                    ref_dict[locus] = lr
                else:
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    elif line[0] == 'P':
                        p = float(line[1])

        return ref_dict, cpi, p

    # getting similar data from each patient's test data
    def get_test_data_parentx(self, doc_name):
        test_dict = {}
        test_cpi = 1
        with open('cognation/scripts/tests/test_cases/parent_cases/' + doc_name, 'r') as parentx:
            for line in parentx:
                parent_formula_dict = ParentFormula.calculate_relation(self, line)
                locus = parent_formula_dict['locus']
                lr = float(parent_formula_dict['lr'])
                test_cpi *= lr
                test_dict[locus] = lr

        test_cpi = int(test_cpi)
        test_p = (test_cpi / (1 + test_cpi)) * 100

        return test_dict, test_cpi, test_p

    for i in range(len(doc_names_list)):
        doc_path = doc_names_list[i]
        overall_ref_dict[doc_path] = get_reference_data_parentx(doc_path)
        overall_test_dict[doc_path] = self.get_test_data_parentx(doc_path)




if __name__ == '__main__':
    unittest.main()
