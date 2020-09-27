# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import django
import unittest
import os
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "genomed_software.settings")
from django.test import TestCase
from ..tests import GetParentsData

# from ...models import Locus
django.setup()

# all possible test cases
doc_names_list = ['parent1/reference_data_parent1.txt', 'parent2/reference_data_parent2.txt', 'parent3/reference_data_parent3_veri.txt']
# places to collect data for all parents
overall_ref_dict = {}
overall_test_dict = {}

class AllTestsCases(GetParentsData):
    # getting reference and test data for all parents
    def all_data(self):
        for i in range(len(doc_names_list)):
            parent_path = doc_names_list[i]
            overall_test_dict[parent_path] = self.get_reference_data_parentx(parent_path)
            overall_ref_dict[parent_path] = self.get_test_data_parentx(parent_path)
        return overall_ref_dict, overall_test_dict

    over



# Here we test all that cases
class TestParentsData(TestCase):
    def setUp(self):
        pass

    def FinallyTest(self):
        for i in range(len(doc_names_list)):
            doc_ref_tuple = overall_ref_dict[doc_names_list[i]]
            cpi_ref = doc_ref_tuple[1]
            p_ref = doc_ref_tuple[2]
            doc_test_tuple = overall_test_dict[doc_names_list[i]]
            cpi_test = doc_test_tuple[1]
            p_test = doc_test_tuple[2]
            self.assertEqual(cpi_ref, cpi_test)
            self.assertEqual(p_ref, p_test)

            doc_ref_dict = doc_ref_tuple[0]
            doc_test_dict = doc_test_tuple[0]
            for key in doc_ref_dict.keys():
                lr_ref = doc_ref_dict[key]
                lr_test = doc_test_dict[key]
                self.assertEqual(lr_ref, lr_test)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
