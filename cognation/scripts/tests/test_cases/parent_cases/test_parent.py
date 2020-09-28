# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import unittest
import re
from cognation.formulas.parent import ParentFormula
from django.test import TestCase

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

                # locus Yindel case - there's no lr
                if len(line) == 3:
                    locus = line[0]
                    lr = '-'
                    ref_dict[locus] = lr

                # other loci - there is int meaning of lr
                elif len(line) == 4:
                    locus = line[0]
                    lr = float(line[3]) * 100 / 100
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
    def get_test_data_parentx(self, doc_name):
        test_dict = {}
        test_cpi = 1
        with open('cognation/scripts/tests/test_cases/parent_cases/' + doc_name, 'r') as parentx:
            for line in parentx:
                line = line.strip()
                parent_formula_dict = self.calculate_relation(re.split(r'[\s\t]+', line))
                locus = parent_formula_dict['locus']
                lr = parent_formula_dict['lr']

                # if it's not Yindel locus
                if lr != '-':
                    test_cpi *= lr
                    lr = float("{0:.2f}".format(lr))
                    test_dict[locus] = lr
                else:
                    test_dict[locus] = lr
                    continue

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

            # There can be changes in the last digit of the large number, so this case would be submitted
            cond_exp = abs(cpi_test - cpi_ref) <=1
            self.assertTrue(cond_exp, True)

            p_ref = int(parent_ref_tuple[2] * 100) / 100
            p_test = int(parent_test_tuple[2] * 100) / 100
            self.assertEqual(p_ref, p_test)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
