from __future__ import unicode_literals
import re
from cognation.formulas.parent import ParentFormula
from cognation.formulas.grandparent import GrandParentFormula
from cognation.formulas.uncle import UncleFormula
from cognation.formulas.sibling import SiblingFormula
from cognation.formulas.cousin import CousinFormula

FORMULA_TYPE_PARENT = 1
FORMULA_TYPE_GRANDPARENT = 2
FORMULA_TYPE_UNCLE = 3
FORMULA_TYPE_SIBLING = 4
FORMULA_TYPE_COUSIN = 5


def formula_usage(number, data):
    if number == FORMULA_TYPE_PARENT:
        return ParentFormula(data)
    if number == FORMULA_TYPE_GRANDPARENT:
        return GrandParentFormula(data)
    if number == FORMULA_TYPE_UNCLE:
        return UncleFormula(data)
    if number == FORMULA_TYPE_SIBLING:
        return SiblingFormula(data)
    if number == FORMULA_TYPE_COUSIN:
        return CousinFormula(data)

# all possible test cases
# doc_refnames_list = ['1', '2', '3']
# doc_testnames_list = ['1', '2', '3']


class GetRefTestData():
    # getting loci and lrs in the dictionary and also CPI and P from each person's reference data
    @staticmethod
    def get_reference_data(short_path, doc_name):
        #  short_path is a path to the reference data file until case folder
        #  doc_name is a path directly to the reference data file from case folder - in parents case we get it from doc_refnames_list

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

                    if lr == '-':
                        ref_dict[locus] = lr
                    else:
                        lr = float(line[3]) * 100 / 100
                        ref_dict[locus] = lr

                # case for getting cpi and p meanings
                if len(line) == 1 and line[0] != '':
                    line = line[0].split('=')

                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    else:
                        p = float(line[1])

        return ref_dict, cpi, p

    # getting similar data from each patient's test data
    def get_test_data_parentx(self, short_path, doc_name, formula_number):
        locus_lr_dict = {}
        with open(short_path + doc_name, 'r') as test_data:
            for line in test_data:
                line = line.strip()

                formula_name = formula_usage(formula_number)
                parent_formula_dict = formula_name.calculate_relation(re.split(r'[\s\t]+', line))

                # parent_formula_dict = self.calculate_relation(re.split(r'[\s\t]+', line))
                locus = parent_formula_dict['locus']
                lr = parent_formula_dict['lr']
                locus_lr_dict[locus] = lr

        return locus, lr

