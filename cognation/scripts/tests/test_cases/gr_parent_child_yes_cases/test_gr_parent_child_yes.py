# -*- coding: utf-8 -*-
import unittest
from django.test import TestCase
from cognation.scripts.tests import GetData
import logging

logger = logging.getLogger('django.db.backends')
logger.setLevel(logging.DEBUG)
logger.addHandler(logging.StreamHandler())


class TestFormula(TestCase):
    def setUp(self):

        self.reference_paths = ['aa_an_an/aa_an_an_ref', 'aa_bb_ab/aa_bb_ab_ref', 'aa_bc_abac/aa_bc_abac_ref',
                                'ab_aa_an/ab_aa_an_ref', 'ab_ab_anbn/ab_ab_anbn_ref', 'ab_ac_anbc/ab_ac_anbc_ref',
                                'ab_cc_acbc/ab_cc_acbc_ref', 'ab_cd_acadbcbd/ab_cd_acadbcbd_ref']
        self.test_paths = ['aa_an_an/aa_an_an_test', 'aa_bb_ab/aa_bb_ab_test', 'aa_bc_abac/aa_bc_abac_test',
                           'ab_aa_an/ab_aa_an_test', 'ab_ab_anbn/ab_ab_anbn_test', 'ab_ac_anbc/ab_ac_anbc_test',
                           'ab_cc_acbc/ab_cc_acbc_test', 'ab_cd_acadbcbd/ab_cd_acadbcbd_test']
        short_path ='cognation/scripts/tests/test_cases/gr_parent_child_yes_cases/'

        get_ref = GetData()
        self.overall_ref_dict = {}
        self.overall_test_dict = {}

        for i in range(len(self.reference_paths)):
            ref_path = self.reference_paths[i]
            self.overall_ref_dict[ref_path] = get_ref.get_reference_data(short_path, ref_path, 3)

            test_path = self.test_paths[i]
            self.overall_test_dict[test_path] = get_ref.get_test_data(short_path, test_path, 19)
        pass

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
                self.assertEqual(lr_ref, lr_test)

            cpi_ref = ref_tuple[1]
            cpi_test = test_tuple[1]
            self.assertEqual(cpi_ref, cpi_test)

            p_ref = int(ref_tuple[2] * 100) / 100
            p_test = int(test_tuple[2] * 100) / 100
            self.assertEqual(p_ref, p_test)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
