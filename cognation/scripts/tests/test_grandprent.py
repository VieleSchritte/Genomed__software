# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import unittest
import re
from ...formulas.grandparent import GrandParentFormula
from django.test import TestCase

# all possible test cases
doc_refnames_list = ['grandparent1/reference_data_grandparent1.txt', 'grandparent2/reference_data_grandparent2.txt']
doc_testnames_list = ['grandparent1/test_data_grandparent1.txt', 'grandparent2/test_data_grandparent2.txt']

overall_ref_dict = {}
overall_test_dict = {}

class GetGrandParentData(GrandParentFormula):
    # getting loci and lrs in the dictionary and also CPI and P from each patient's reference data
    @staticmethod
    def get_reference_data_grandparentx(doc_name):
        ref_dict = {}
        cpi = 0
        p = 0.0
        with open('cognation/scripts/tests/test_cases/grandparent_cases/' + doc_name, 'r') as grandparentx_data:
            for line in grandparentx_data:
                line = line.strip().split('\t')

                # locus Yindel case - there's no lr
                if len(line[0]) == 3:
                    locus = line[0]
                    lr = ' '
                    ref_dict[locus] = lr

                # other loci - there is int meaning of lr
                elif len(line) == 4:
                    locus = line[0]
                    lr = line[3]
                    ref_dict[locus] = lr

                # case for getting cpi and p meanings
                elif len(line) == 1 and line[0] != '':
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    elif line[0] == 'P':
                        p = float(line[1])

        return ref_dict, cpi, p

    # getting similar data from each patient's test data
    def get_test_data_grandparentx(self, doc_name):
        test_dict = {}
        test_cpi = 1
        with open('cognation/scripts/tests/test_cases/grandparent_cases/' + doc_name, 'r') as parentx:
            for line in parentx:
                line = line.strip()
                grandparent_formula_dict = self.calculate_relation(re.split(r'[\s\t]+', line))
                locus = grandparent_formula_dict['locus']
                lr = grandparent_formula_dict['lr']
                test_dict[locus] = lr

                # if it's not Yindel locus
                if lr != '-':
                    lr = float(lr)
                    test_cpi *= lr
                else:
                    continue

        test_cpi = int(test_cpi)
        test_p = (test_cpi / (1 + test_cpi)) * 100

        return test_dict, test_cpi, test_p

    def prep(self):
        for i in range(len(doc_refnames_list)):
            doc_ref_path = doc_refnames_list[i]
            overall_ref_dict[doc_ref_path] = GetGrandParentData.get_reference_data_grandparentx(doc_ref_path)
            doc_test_path = doc_testnames_list[i]
            overall_test_dict[doc_test_path] = self.get_test_data_grandparentx(doc_test_path)

instance = GetGrandParentData(GrandParentFormula)
instance.prep()

class TestGrandParentFormula(TestCase):
    def setUp(self):
        pass

    def test_final_assertion(self):
        for i in range(len(doc_refnames_list)):
            doc_ref_path = doc_refnames_list[i]
            doc_test_path = doc_testnames_list[i]

            grandparent_ref_tuple = overall_ref_dict[doc_ref_path]
            grandparent_test_tuple = overall_test_dict[doc_test_path]

            cpi_ref = grandparent_ref_tuple[1]
            cpi_test = grandparent_test_tuple[1]
            self.assertEqual(cpi_ref, cpi_test)

            p_ref = grandparent_ref_tuple[2]
            p_test = grandparent_test_tuple[2]
            self.assertEqual(p_ref, p_test)

            dict_loci_lrs_ref = grandparent_ref_tuple[0]
            dict_loci_lrs_test = grandparent_test_tuple[0]
            for key in dict_loci_lrs_ref.keys():
                lr_ref = dict_loci_lrs_ref[key]
                lr_test = dict_loci_lrs_test[key]
                self.assertEqual(lr_ref, lr_test)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
