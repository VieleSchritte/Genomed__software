# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import django
import unittest
#import os
import re
#os.environ.setdefault("DJANGO_SETTINGS_MODULE", "genomed_software.settings")
from ...formulas.parent import ParentFormula
from django.test import TestCase
django.setup()

# all possible test cases
doc_refnames_list = ['parent1/reference_data_parent1.txt', 'parent2/reference_data_parent2.txt', 'parent3/reference_data_parent3_veri.txt']
doc_testnames_list = ['parent1/test_data_parent1.txt', 'parent2/test_data_parent2.txt', 'parent3/test_data_parent3_veri.txt']

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
            print(doc_name)
            for line in parentx:
                print(line)
                print(re.split(r'[\s\t]+', line))
                parent_formula_dict = self.calculate_relation(re.split(r'[\s\t]+', line))
                locus = parent_formula_dict['locus']
                lr = float(parent_formula_dict['lr'])
                test_cpi *= lr
                test_dict[locus] = lr

        test_cpi = int(test_cpi)
        test_p = (test_cpi / (1 + test_cpi)) * 100

        return test_dict, test_cpi, test_p

    def prep(self):
        for i in range(len(doc_refnames_list)):
            doc_ref_path = doc_refnames_list[i]
            overall_ref_dict[doc_ref_path] = GetParentsData.get_reference_data_parentx(doc_ref_path)
            doc_test_path = doc_testnames_list[i]
            overall_test_dict[doc_test_path] = self.get_test_data_parentx(doc_test_path)

instance = GetParentsData(ParentFormula)
instance.prep()

class TestParentFormula(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def FinalAssertion(self):
        for i in range(len(doc_refnames_list)):
            doc_ref_path = doc_refnames_list[i]
            doc_test_path = doc_testnames_list[i]

            parent_ref_tuple = overall_ref_dict[doc_ref_path]
            parent_test_tuple = overall_test_dict[doc_test_path]

            cpi_ref = parent_ref_tuple[1]
            cpi_test = parent_test_tuple[1]
            self.assertEqual(cpi_ref, cpi_test)

            p_ref = parent_ref_tuple[2]
            p_test = parent_test_tuple[2]
            self.assertEqual(p_ref, p_test)

            dict_loci_lrs_ref = parent_ref_tuple[0]
            dict_loci_lrs_test = parent_test_tuple[0]
            for key in dict_loci_lrs_ref.keys():
                lr_ref = dict_loci_lrs_ref[key]
                lr_test = dict_loci_lrs_test[key]
                self.assertEqual(lr_ref, lr_test)


if __name__ == '__main__':
    unittest.main()
    