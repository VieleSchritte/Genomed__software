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
            7: [],
            8: OneKnownSupposedFormula,
            9: [],
            10: [],
            11: StepbrotherFormula,
            12: BrotherFormula,
            13: [],
            14: CousinFormula,
            15: UncleFormula,
            16: [],
            17: GrandParentFormula,
            18: [],
            19: [],
            20: [],
            21: []
        }

        for key in num_to_formula.keys():
            if number == key:
                return num_to_formula[key](Formula)

    def get_test_data(self, short_path, doc_name, number):
        test_cpi = 1
        test_dict = {}
        call_command("loaddata", "converted.json", verbosity=0)

        with open(short_path + doc_name, 'r') as test_data:
            for line in test_data:
                line = line.strip()
                case_formula = self.formula_usage(number)
                formula_dict = case_formula.calculate_relation(re.split(r'[\s\t]+', line))
                locus = formula_dict['locus']
                lr = formula_dict['lr']

                if lr != '-':
                    test_cpi *= lr
                    lr = float("{0:.2f}".format(lr))
                    test_dict[locus] = lr

                #  case of gender specific loci
                else:
                    test_dict[locus] = lr
                    continue

            test_cpi = float(test_cpi)
            test_p = (test_cpi / (1 + test_cpi)) * 100
            test_cpi = round(test_cpi)

            return test_dict, test_cpi, test_p

    @staticmethod
    def get_reference_data(short_path, doc_name, part_number):
        ref_dict = {}
        cpi = 0
        p = 0.0
        with open(short_path + doc_name, 'r') as ref_data:
            for line in ref_data:
                line = line.strip().split('\t')

                # not gender specific loci - there is int meaning of lr
                if len(line) != 1:
                    locus = line[0]
                    lr = line[-1]
                    if part_number == 3:
                        lr = line[4]

                    #  case of gender specific loci
                    if lr == '-':
                        ref_dict[locus] = lr
                        continue
                    # case of int lr meanings of loci
                    else:
                        lr = float(lr) * 100 / 100
                        ref_dict[locus] = lr

                # case for getting cpi and p meanings
                elif len(line) == 1 and line[0] != '':
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    elif line[0] == 'P':
                        p = float(line[1])

        return ref_dict, cpi, p

