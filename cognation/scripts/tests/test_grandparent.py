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

class GetGrandParentsData(GrandParentFormula):
    # getting loci and lrs in the dictionary and also CPI and P from each patient's reference data
    @staticmethod
    def get_reference_data_grandparentx(doc_name):
        ref_dict = {}
        cpi = 0
        p = 0.0
        with open('cognation/scripts/tests/test_cases/grandparent_cases/' + doc_name, 'r') as grandparentx_data:
            print('DOCNAME: ', doc_name)
            for line in grandparentx_data:
                line = line.strip().split('\t')
                locus = line[0]

                # There is no case to print the Yindel result, so now we skip it
                if locus == '' or locus[0] == 'Y':
                    continue

                # case for getting cpi and p meanings
                elif locus[0] == 'C' or locus[0] == 'P':
                    locus_split = locus.split('=')
                    if locus_split[0] == 'CPI':
                        cpi = int(locus_split[1])
                    elif locus_split[0] == 'P':
                        p = float(locus_split[1])

                # In other cases we just collect lrs
                else:
                    lr = float(line[3]) * 100 / 100
                    ref_dict[locus] = lr

        return ref_dict, cpi, p

    # getting similar data from each grandparent's test data
    def get_test_data_grandparentx(self, doc_name):
        test_dict = {}
        test_cpi = 1
        with open('cognation/scripts/tests/test_cases/grandparent_cases/' + doc_name, 'r') as grandparentx:
            for line in grandparentx:
                line = line.strip()

                # Yindel case
                if line[0] == 'Y':
                    continue

                grandparent_formula_dict = self.calculate_relation(re.split(r'[\s\t]+', line))
                locus = grandparent_formula_dict['locus']
                lr = grandparent_formula_dict['lr']

                # if it's not Yindel locus
                if lr != '-':
                    test_cpi *= lr
                    lr = float("{0:.2f}".format(lr))
                    test_dict[locus] = lr
                else:
                    test_dict[locus] = lr
                    continue

        test_cpi = float(test_cpi)
        test_p = (test_cpi / (1 + test_cpi)) * 100
        test_cpi = round(test_cpi)

        return test_dict, test_cpi, test_p

    def prep(self):
        for i in range(len(doc_refnames_list)):
            doc_ref_path = doc_refnames_list[i]
            overall_ref_dict[doc_ref_path] = GetGrandParentsData.get_reference_data_grandparentx(doc_ref_path)
            doc_test_path = doc_testnames_list[i]
            overall_test_dict[doc_test_path] = self.get_test_data_grandparentx(doc_test_path)

instance = GetGrandParentsData(GrandParentFormula)
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

            dict_loci_lrs_ref = grandparent_ref_tuple[0]
            dict_loci_lrs_test = grandparent_test_tuple[0]

            for key in dict_loci_lrs_ref.keys():
                lr_ref = dict_loci_lrs_ref[key]
                lr_test = dict_loci_lrs_test[key]
                self.assertEqual(lr_ref, lr_test)

            cpi_ref = grandparent_ref_tuple[1]
            cpi_test = grandparent_test_tuple[1]

            # There can be changes in the last digit of the large number, so this case would be submitted
            cond_exp = abs(cpi_test - cpi_ref) <=1
            self.assertTrue(cond_exp, True)

            p_ref = int(grandparent_ref_tuple[2] * 100) / 100
            p_test = int(grandparent_test_tuple[2] * 100) / 100
            self.assertEqual(p_ref, p_test)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()