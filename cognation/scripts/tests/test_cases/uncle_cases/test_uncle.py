# -*- coding: utf-8 -*-
import unittest
from django.test import TestCase
from cognation.scripts.tests import GetData

# all possible test cases
doc_refnames_list = ['Aunt1/aunt_ref.txt']
doc_testnames_list = ['Aunt1/aunt_test.txt']
short_path = 'cognation/scripts/tests/test_cases/uncle_cases/'

overall_ref_dict = {}
overall_test_dict = {}

UNCLE_TYPE = 3


class GetUncleData:
    #  preparing dictionaries for assertion
    def prep(self):
        get_ref = GetData()
        for i in range(len(doc_refnames_list)):
            doc_ref_path = doc_refnames_list[i]
            overall_ref_dict[doc_ref_path] = get_ref.get_reference_data(short_path, doc_ref_path)

            doc_test_path = doc_testnames_list[i]
            overall_test_dict[doc_test_path] = get_ref.get_test_data(short_path, doc_test_path, UNCLE_TYPE)


instance = GetUncleData()
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
