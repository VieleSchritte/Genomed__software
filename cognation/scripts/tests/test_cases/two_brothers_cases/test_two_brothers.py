import unittest
from django.test import TestCase
from cognation.scripts.tests import GetData
# import logging
# logger = logging.getLogger('django.db.backends')
# logger.setLevel(logging.DEBUG)
# logger.addHandler(logging.StreamHandler())


class TestFormula(TestCase):
    def setUp(self):
        self.reference_paths = [
            'aa_ab_cccd/aa_ab_cccd_ref',
            'ab_aa_cccd/ab_aa_cccd_ref',
            'ab_ac_ad/ab_ac_ad_ref',
            'ab_notab_notab/ab_notab_notab_ref',
            #'aa_aa_any/aa_aa_any_ref',
            #'ab_ab_any/ab_ab_any_ref',
            #'ab_aa_bn/ab_aa_bn_ref',
            #'aa_ab_ac/aa_ab_ac_ref',
            #'aa_ab_bc/aa_ab_bc_ref',
            #'ab_aa_ac/ab_aa_ac_ref',
            #'ab_ac_bc/ab_ac_bc_ref',
            #'ab_ac_bd/ab_ac_bd_ref',
            #'ab_ac_cc/ab_ac_cc_ref',
            #'ab_ac_cd/ab_ac_cd_ref',
            #'call_BrotherFormula/call_BrotherFormula_ref'
        ]

        self.test_paths = [
            'aa_ab_cccd/aa_ab_cccd_test',
            'ab_aa_cccd/ab_aa_cccd_test',
            'ab_ac_ad/ab_ac_ad_test',
            'ab_notab_notab/ab_notab_notab_test',
            #'aa_aa_any/aa_aa_any_test',
            #'ab_ab_any/ab_ab_any_test',
            #'ab_aa_bn/ab_aa_bn_test',
            #'aa_ab_ac/aa_ab_ac_test',
            #'aa_ab_bc/aa_ab_bc_test',
            #'ab_aa_ac/ab_aa_ac_test',
            #'ab_ac_bc/ab_ac_bc_test',
            #'ab_ac_bd/ab_ac_bd_test',
            #'ab_ac_cc/ab_ac_cc_test',
            #'ab_ac_cd/ab_ac_cd_test',
            #'call_BrotherFormula/call_BrotherFormula_test'
        ]

        short_path = 'cognation/scripts/tests/test_cases/two_brothers_cases/'
        get_ref = GetData()
        self.overall_ref_dict, self.overall_test_dict = {}, {}

        for i in range(len(self.reference_paths)):
            ref_path, test_path = self.reference_paths[i], self.test_paths[i]
            self.overall_ref_dict[ref_path] = get_ref.get_reference_data(short_path, ref_path, 3)
            self.overall_test_dict[test_path] = get_ref.get_test_data(short_path, test_path, 13)

    def test_formula(self):
        for i in range(len(self.reference_paths)):
            ref_path, test_path = self.reference_paths[i], self.test_paths[i]
            ref_tuple, test_tuple = self.overall_ref_dict[ref_path], self.overall_test_dict[test_path]
            dict_loci_lrs_ref, dict_loci_lrs_test = ref_tuple[0], test_tuple[0]

            for key in dict_loci_lrs_ref.keys():
                lr_ref, lr_test = dict_loci_lrs_ref[key], dict_loci_lrs_test[key]
                self.assertEqual(lr_ref, lr_test, key)

            cpi_ref, cpi_test = ref_tuple[1], test_tuple[1]
            p_ref, p_test = float(ref_tuple[2]), float(test_tuple[2])
            self.assertEqual(cpi_ref, cpi_test)
            self.assertEqual(p_ref, p_test)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()
