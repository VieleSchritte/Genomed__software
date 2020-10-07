# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import unittest
import re
from cognation.formulas.uncle import UncleFormula
from django.test import TestCase

# all possible test cases
doc_refnames_list = ['Aunt1/aunt_ref.txt']
doc_testnames_list = ['Aunt1/aunt_test.txt']

short_path = 'cognation/scripts/tests/test_cases/uncle_cases/'

overall_ref_dict = {}
overall_test_dict = {}


class GetUncleData(UncleFormula):
    # getting loci and lrs in the dictionary and also CPI and P from each patient's reference data
    @staticmethod
    def get_reference_data(doc_name):
        ref_dict = {}
        cpi = 0
        p = 0.0
        with open(short_path + doc_name, 'r') as ref_data:
            for line in ref_data:
                line = line.strip().split('\t')

                if len(line) == 4:
                    locus = line[0]
                    lr = line[3]

                    if lr == '-':
                        ref_dict[locus] = lr
                    else:
                        lr = float(lr)
                        ref_dict[locus] = lr

                # case for getting cpi and p meanings
                if len(line) == 1 and line[0] != '':
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    elif line[0] == 'P':
                        p = float(line[1])
        return ref_dict, cpi, p

    # getting similar data from each patient's test data
    def get_test_data(self, doc_name):
        test_dict = {}
        test_cpi = 1
        with open(short_path + doc_name, 'r') as test_data:
            for line in test_data:
                line = line.strip()
                uncle_formula_dict = self.calculate_relation(re.split(r'[\s\t]+', line))
                locus = uncle_formula_dict['locus']
                lr = uncle_formula_dict['lr']

                # if it's not gender-specific locus we can calculate cpi using lr meaning
                if lr == '-':
                    test_dict[locus] = lr
                    continue
                else:
                    test_cpi *= lr
                    lr = float("{0:.2f}".format(lr))
                    test_dict[locus] = lr

        test_p = (test_cpi / (1 + test_cpi)) * 100
        test_cpi = int(test_cpi)

        return test_dict, test_cpi, test_p

    def prep(self):
        for i in range(len(doc_refnames_list)):
            doc_ref_path = doc_refnames_list[i]
            overall_ref_dict[doc_ref_path] = GetUncleData.get_reference_data(doc_ref_path)
            doc_test_path = doc_testnames_list[i]
            overall_test_dict[doc_test_path] = self.get_test_data(doc_test_path)


instance = GetUncleData(UncleFormula)
instance.prep()


class TestUncleFormula(TestCase):
    def setUp(self):
        pass

    def test_final_assertion(self):
        for i in range(len(doc_refnames_list)):
            doc_ref_path = doc_refnames_list[i]
            doc_test_path = doc_testnames_list[i]

            uncle_ref_tuple = overall_ref_dict[doc_ref_path]
            uncle_test_tuple = overall_test_dict[doc_test_path]

            dict_loci_lrs_ref = uncle_ref_tuple[0]
            dict_loci_lrs_test = uncle_test_tuple[0]

            for key in dict_loci_lrs_ref.keys():
                lr_ref = dict_loci_lrs_ref[key]
                lr_test = dict_loci_lrs_test[key]
                self.assertEqual(lr_ref, lr_test)

            cpi_ref = uncle_ref_tuple[1]
            cpi_test = uncle_test_tuple[1]
            self.assertEqual(cpi_ref, cpi_test)

            p_ref = int(uncle_ref_tuple[2] * 100) / 100
            p_test = int(uncle_test_tuple[2] * 100) / 100
            self.assertEqual(p_ref, p_test)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()