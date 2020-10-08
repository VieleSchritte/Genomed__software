from cognation.formulas.for_2_participants.parent import ParentFormula
from cognation.formulas.for_2_participants.grandparent import GrandParentFormula
from cognation.formulas.for_2_participants.uncle import UncleFormula
from cognation.formulas.for_2_participants.cousin import CousinFormula
from cognation.formulas.for_3_participants.sibling import SiblingFormula
from cognation.formulas.for_2_participants.base import Formula
import re

PARENT_TYPE = 1
GRANDPARENT_TYPE = 2
UNCLE_TYPE = 3
COUSIN_TYPE = 4
SIBLING_TYPE = 5


class GetData:
    @staticmethod
    def get_reference_data(short_path, doc_name):
        ref_dict = {}
        cpi = 0
        p = 0.0
        with open(short_path + doc_name, 'r') as ref_data:
            for line in ref_data:
                line = line.strip().split('\t')

                # other loci - there is int meaning of lr
                if len(line) == 4:
                    locus = line[0]
                    lr = line[3]

                    #  case of gender specific loci
                    if lr == '-':
                        ref_dict[locus] = lr
                        continue
                    # case of int lr meanings of loci
                    else:
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

    @staticmethod
    def formula_usage(number):
        if number == PARENT_TYPE:
            return ParentFormula(Formula)
        if number == GRANDPARENT_TYPE:
            return GrandParentFormula(Formula)
        if number == UNCLE_TYPE:
            return UncleFormula(Formula)
        if number == COUSIN_TYPE:
            return CousinFormula(Formula)
        if number == SIBLING_TYPE:
            return SiblingFormula(Formula)

    def get_test_data(self, short_path, doc_name, number):
        test_cpi = 1
        test_dict = {}

        with open(short_path + doc_name, 'r') as test_data:
            for line in test_data:
                line = line.strip()
                case_formula = self.formula_usage(number)
                parent_formula_dict = case_formula.calculate_relation(re.split(r'[\s\t]+', line))
                locus = parent_formula_dict['locus']
                lr = parent_formula_dict['lr']

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
