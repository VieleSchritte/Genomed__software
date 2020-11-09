import unittest
from django.test import TestCase
from cognation.scripts.tests import GetData
# import logging

# logger = logging.getLogger('django.db.backends')
# logger.setLevel(logging.DEBUG)
# logger.addHandler(logging.StreamHandler())


class TestFormula(TestCase):
    def setUp(self):
        self.reference_paths = ['no_intersections/no_intersections_ref', 'same_genotypes/same_genotypes_ref',
                                'aa_ab_ab_an/aa_ab_ab_an_ref', 'aa_ab_an_ab/aa_ab_an_ab_ref',
                                'aa_bb_ab_ab/aa_bb_ab_ab_ref', 'aa_bc_ab_ac/aa_bc_ab_ac_ref',
                                'ab_ac_ab_acbc/ab_ac_ab_acbc_ref', 'ab_ac_an_bc/ab_ac_an_bc_ref',
                                'ab_ac_bc_an/ab_ac_bc_an_ref', 'ab_cd_bd_ac/ab_cd_bd_ac_ref']

        self.test_paths = ['no_intersections/no_intersections_test', 'same_genotypes/same_genotypes_test',
                           'aa_ab_ab_an/aa_ab_ab_an_test', 'aa_ab_an_ab/aa_ab_an_ab_test',
                           'aa_bb_ab_ab/aa_bb_ab_ab_test', 'aa_bc_ab_ac/aa_bc_ab_ac_test',
                           'ab_ac_ab_acbc/ab_ac_ab_acbc_test', 'ab_ac_an_bc/ab_ac_an_bc_test',
                           'ab_ac_bc_an/ab_ac_bc_an_test', 'ab_cd_bd_ac/ab_cd_bd_ac_test']

        short_path = 'cognation/scripts/tests/test_cases/two_known_supposed_cases/'

        get_ref = GetData()
        self.overall_ref_dict = {}
        self.overall_test_dict = {}

        for i in range(len(self.reference_paths)):
            ref_path = self.reference_paths[i]
            self.overall_ref_dict[ref_path] = get_ref.get_reference_data(short_path, ref_path, 4)

            test_path = self.test_paths[i]
            self.overall_test_dict[test_path] = get_ref.get_test_data(short_path, test_path, 5)

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

            p_ref = float(ref_tuple[2])
            p_test = float(test_tuple[2])
            self.assertEqual(p_ref, p_test)

    def tearDown(self):
        pass


if __name__ == '__main__':
    unittest.main()