# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import unittest
import re
from cognation.formulas.parent import ParentFormula
from django.test import TestCase
from cognation.scripts.tests import GetData

# all possible test cases
doc_refnames_list = ['parent1/reference_data_parent1.txt', 'parent2/reference_data_parent2.txt', 'parent3/reference_data_parent3_veri.txt']
doc_testnames_list = ['parent1/test_data_parent1.txt', 'parent2/test_data_parent2.txt', 'parent3/test_data_parent3_veri.txt']

short_path = 'cognation/scripts/tests/test_cases/parent_cases/'

overall_ref_dict = {}
overall_test_dict = {}


class GetParentsData(ParentFormula):
    # getting test data from each parent's file
    def get_test_data(self, doc_name):
        test_dict = {}
        test_cpi = 1
        with open(short_path + doc_name, 'r') as test_data:
            for line in test_data:
                line = line.strip()
                parent_formula_dict = self.calculate_relation(re.split(r'[\s\t]+', line))
                locus = parent_formula_dict['locus']
                lr = parent_formula_dict['lr']

                if lr != '-':
                    test_cpi *= lr
                    lr = float("{0:.2f}".format(lr))
                    test_dict[locus] = lr

                #  case of gender specific loci
                else:
                    test_dict[locus] = lr
                    continue

        test_cpi = float(test_cpi)
        test_p = (test_cpi / (1 + test_cpi)) * 100
        test_cpi = round(test_cpi)

        return test_dict, test_cpi, test_p

    #  preparing dictionaries for assertion
    def prep(self):
        get_ref = GetData()
        for i in range(len(doc_refnames_list)):
            doc_ref_path = doc_refnames_list[i]
            overall_ref_dict[doc_ref_path] = get_ref.get_reference_data(short_path, doc_ref_path)

            doc_test_path = doc_testnames_list[i]
            overall_test_dict[doc_test_path] = self.get_test_data(doc_test_path)


instance = GetParentsData(ParentFormula)
instance.prep()


class TestParentFormula(TestCase):
    def setUp(self):
        pass

    def test_final_assertion(self):
        for i in range(len(doc_refnames_list)):
            doc_ref_path = doc_refnames_list[i]
            doc_test_path = doc_testnames_list[i]

            parent_ref_tuple = overall_ref_dict[doc_ref_path]
            parent_test_tuple = overall_test_dict[doc_test_path]

            dict_loci_lrs_ref = parent_ref_tuple[0]
            dict_loci_lrs_test = parent_test_tuple[0]

            for key in dict_loci_lrs_ref.keys():
                lr_ref = dict_loci_lrs_ref[key]
                lr_test = dict_loci_lrs_test[key]
                self.assertEqual(lr_ref, lr_test)

            cpi_ref = parent_ref_tuple[1]
            cpi_test = parent_test_tuple[1]
            self.assertEqual(cpi_ref, cpi_test)

            p_ref = int(parent_ref_tuple[2] * 100) / 100
            p_test = int(parent_test_tuple[2] * 100) / 100
            self.assertEqual(p_ref, p_test)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
