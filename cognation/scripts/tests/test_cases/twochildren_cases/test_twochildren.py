# -*- coding: utf-8 -*-
import unittest
from django.test import TestCase
from cognation.scripts.tests import GetData
# import logging
"""
logger = logging.getLogger('django.db.backends')
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())
"""


class TestTwoChildrenFormula(TestCase):
    def setUp(self):

        self.reference_paths = ['1_calls_ParentHardy/1_calls_ParentHardy_ref']
        #self.reference_paths = ['1_calls_ParentHardy/1_calls_ParentHardy_ref', '2_aa_an_an/2_aa_an_an_ref',
                               # '3_aa_bb_ab/3_aa_bb_ab_ref', '4_aa_bc_abac/4_aa_bc_abac_ref', '5_ab_ac_anbc/5_ab_ac_anbc_ref',
                               # '6_ab_cc_acbc/6_ab_cc_acbc_ref', '7_ab_cd_acadbcbd/7_ab_cd_acadbcbd_ref']
        self.test_paths = ['1_calls_ParentHardy/1_calls_ParentHardy_test']
        #self.test_paths = ['1_calls_ParentHardy/1_calls_ParentHardy_test', '2_aa_an_an/2_aa_an_an_test',
                           #'3_aa_bb_ab/3_aa_bb_ab_test', '4_aa_bc_abac/4_aa_bc_abac_test', '5_ab_ac_anbc/5_ab_ac_anbc_test',
                          # '6_ab_cc_acbc/6_ab_cc_acbc_test', '7_ab_cd_acadbcbd/7_ab_cd_acadbcbd_test']
        short_path = 'cognation/scripts/tests/test_cases/twochildren_cases/'

        get_ref = GetData()
        self.overall_ref_dict = {}
        self.overall_test_dict = {}

        for i in range(len(self.reference_paths)):
            ref_path = self.reference_paths[i]
            self.overall_ref_dict[ref_path] = get_ref.get_reference_data(short_path, ref_path, 3)

            test_path = self.test_paths[i]
            self.overall_test_dict[test_path] = get_ref.get_test_data(short_path, test_path, 8)

    def test_formula(self):
        for i in range(len(self.reference_paths)):
            ref_path = self.reference_paths[i]
            test_path = self.test_paths[i]

            ref_tuple = self.overall_ref_dict[ref_path]
            test_tuple = self.overall_test_dict[test_path]

            dict_loci_lrs_ref = ref_tuple[0]
            dict_loci_lrs_test = test_tuple[0]

            for key in dict_loci_lrs_ref.keys():
                lr_ref = dict_loci_lrs_ref[key]
                lr_test = dict_loci_lrs_test[key]
                self.assertEqual(lr_ref, lr_test, key)

            cpi_ref = ref_tuple[1]
            cpi_test = test_tuple[1]
            self.assertEqual(cpi_ref, cpi_test)

            # p_ref = ref_tuple[2]
            # p_test = test_tuple[2]
            p_ref = int(ref_tuple[2] * 100) / 100
            p_test = int(test_tuple[2] * 100) / 100
            self.assertEqual(p_ref, p_test)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()