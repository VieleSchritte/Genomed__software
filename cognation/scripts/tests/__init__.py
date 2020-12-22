from cognation.formulas.parent import ParentFormula
from cognation.formulas.grandparent import GrandParentFormula
from cognation.formulas.uncle import UncleFormula
from cognation.formulas.cousin import CousinFormula
from cognation.formulas.sibling import SiblingFormula
from cognation.formulas.base import Formula
from cognation.formulas.brother import BrotherFormula
from cognation.formulas.stepbrother import StepbrotherFormula
from cognation.formulas.two_children import TwoChildrenFormula
from cognation.formulas.three_children import ThreeChildrenFormula
from cognation.formulas.one_known_supposed import OneKnownSupposedFormula
from cognation.formulas.two_known_supposed import TwoKnownSupposedFormula
from cognation.formulas.three_known_supposed import ThreeKnownSupposed
from cognation.formulas.couple import CoupleFormula
from cognation.formulas.two_couple import TwoCoupleFormula
from cognation.formulas.two_brothers import TwoBrothersFormula
from cognation.formulas.grand_parent_yes import YesGrandParent
from cognation.formulas.grand_parent_no import NoGrandParent
from cognation.formulas.three_couple import ThreeCoupleFormula
import re
from django.core.management import call_command


class GetData:
    @staticmethod
    def formula_usage(number):
        num_to_formula = {
            1: ParentFormula,
            2: TwoChildrenFormula,
            3: SiblingFormula,
            4: ThreeChildrenFormula,
            5: TwoKnownSupposedFormula,
            6: ThreeKnownSupposed,
            7: CoupleFormula,
            8: OneKnownSupposedFormula,
            9: TwoCoupleFormula,
            10: ThreeCoupleFormula,
            11: StepbrotherFormula,
            12: BrotherFormula,
            13: TwoBrothersFormula,
            14: CousinFormula,
            15: UncleFormula,
            16: [],
            17: GrandParentFormula,
            18: [],
            19: YesGrandParent,
            20: NoGrandParent,
            21: []
        }

        for key in num_to_formula.keys():
            if number == key:
                return num_to_formula[key](Formula)

    def get_test_data(self, short_path, doc_name, number):
        test_cpi, test_dict = 1, {}
        call_command("loaddata", "converted.json", verbosity=0)

        with open(short_path + doc_name, 'r') as test_data:
            for line in test_data:
                line = line.strip().replace(',', '.')
                case_formula = self.formula_usage(number)
                formula_dict = case_formula.calculate_relation(re.split(r'[\s\t]+', line))
                locus, lr = formula_dict['locus'], formula_dict['lr']

                if lr != '-':
                    lr = round(lr, 2)
                    test_cpi *= lr
                    test_dict[locus] = lr
                else:  # case of gender specific loci
                    test_dict[locus] = lr
                    continue

            test_cpi = float(test_cpi)
            test_p = (test_cpi / (1 + test_cpi)) * 100
            test_cpi = round(test_cpi)
            return test_dict, test_cpi, test_p

    @staticmethod
    def get_reference_data(short_path, doc_name, part_number):
        ref_dict, cpi, p = {}, 0, 0.0
        with open(short_path + doc_name, 'r') as ref_data:
            for line in ref_data:
                line = line.strip().split('\t')
                if len(line) != 1:  # not gender specific loci - there is int meaning of lr
                    locus, lr = line[0], line[-1]
                    if part_number == 3:
                        lr = line[4]
                    if lr == '-':  # case of gender specific loci
                        ref_dict[locus] = lr
                        continue
                    else:  # case of int lr meanings of loci
                        lr = float(lr) * 100 / 100
                        ref_dict[locus] = lr
                elif len(line) == 1 and line[0] != '':  # case for getting cpi and p meanings
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    p = float(line[1])
        return ref_dict, cpi, p
