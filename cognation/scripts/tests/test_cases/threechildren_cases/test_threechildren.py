# -*- coding: utf-8 -*-
import unittest
from django.test import TestCase
from cognation.scripts.tests import GetData
# import logging

# logger = logging.getLogger('django.db.backends')
# logger.setLevel(logging.DEBUG)
# logger.addHandler(logging.StreamHandler())


class TestSiblingFormula(TestCase):
    def setUp(self):
        self.reference_paths = ['aa_ab_ac_an/aa_ab_ac_an_ref', 'aa_ab_bb_ab/aa_ab_bb_ab_ref', 'aa_ab_bc_abac/aa_ab_bc_abac_ref',
                                'aa_ab_cc_ac/aa_ab_cc_ac_ref', 'aa_bc_bn_ab/aa_bc_bn_ab_ref']
        self.test_paths = ['aa_ab_ac_an/aa_ab_ac_an_test', 'aa_ab_bb_ab/aa_ab_bb_ab_test', 'aa_ab_bc_abac/aa_ab_bc_abac_test',
                           'aa_ab_cc_ac/aa_ab_cc_ac_test', 'aa_bc_bn_ab/aa_bc_bn_ab_test']
        short_path = 'cognation/scripts/tests/test_cases/threechildren_cases/'

        get_ref = GetData()
        self.overall_ref_dict = {}
        self.overall_test_dict = {}

        for i in range(len(self.reference_paths)):
            ref_path = self.reference_paths[i]
            self.overall_ref_dict[ref_path] = get_ref.get_reference_data(short_path, ref_path, 4)

            test_path = self.test_paths[i]
            self.overall_test_dict[test_path] = get_ref.get_test_data(short_path, test_path, 9)

    def test_formula(self):
        for i in range(len(self.reference_paths)):
            ref_path = self.reference_paths[i]
            test_path = self.test_paths[i]

            sibling_ref_tuple = self.overall_ref_dict[ref_path]
            sibling_test_tuple = self.overall_test_dict[test_path]

            dict_loci_lrs_ref = sibling_ref_tuple[0]
            dict_loci_lrs_test = sibling_test_tuple[0]

            for key in dict_loci_lrs_ref.keys():
                lr_ref = dict_loci_lrs_ref[key]
                lr_test = dict_loci_lrs_test[key]
                self.assertEqual(lr_ref, lr_test, key)

            cpi_ref = sibling_ref_tuple[1]
            cpi_test = sibling_test_tuple[1]
            self.assertEqual(cpi_ref, cpi_test)

            p_ref = int(sibling_ref_tuple[2] * 100) / 100
            p_test = int(sibling_test_tuple[2] * 100) / 100
            self.assertEqual(p_ref, p_test)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()