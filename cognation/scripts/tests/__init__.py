from ...formulas.parent import ParentFormula
from ...formulas.cousin import CousinFormula
from ...formulas.grandparent import GrandParentFormula
from ...formulas.sibling import SiblingFormula
from ...formulas.uncle import UncleFormula


doc_names_list = []
overall_test_dict = {}
overall_ref_dict = {}
# Here we get all data for all test cases for paternity
class GetParentsData(ParentFormula):
    # getting loci and lrs in the dictionary and also CPI and P from each patient's reference data
    @staticmethod
    def get_reference_data_parentx(doc_name):
        ref_dict = {}
        cpi = 0
        p = 0.0
        with open('cognation/scripts/tests/test_cases/parent_cases/' + doc_name, 'r') as parentx_data:
            for line in parentx_data:
                line = line.strip().split('\t')
                if len(line) != 1:
                    locus = line[0]
                    lr = line[3]
                    ref_dict[locus] = lr
                else:
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    elif line[0] == 'P':
                        p = float(line[1])

        return ref_dict, cpi, p

    # getting similar data from each patient's test data
    def get_test_data_parentx(self, doc_name):
        test_dict = {}
        test_cpi = 1
        with open('cognation/scripts/tests/test_cases/parent_cases/' + doc_name, 'r') as parentx:
            for line in parentx:
                parent_formula_dict = self.calculate_relation(line)
                locus = parent_formula_dict['locus']
                lr = float(parent_formula_dict['lr'])
                test_cpi *= lr
                test_dict[locus] = lr

        test_cpi = int(test_cpi)
        test_p = (test_cpi / (1 + test_cpi)) * 100

        return test_dict, test_cpi, test_p

    # getting reference and test data for all parents
    for i in range(len(doc_names_list)):
        parent_path = doc_names_list[i]
        overall_test_dict[parent_path] = get_reference_data_parentx(parent_path)
        overall_ref_dict[parent_path] = get_test_data_parentx(parent_path)

# Here we get all data for all test cases for cousins
class GetCousinData(CousinFormula):
    # getting loci and lrs in the dictionary and also CPI and P from each patient's reference data
    @staticmethod
    def get_reference_data_cousinx(doc_name):
        ref_dict = {}
        cpi = 0
        p = 0.0
        with open('cognation/scripts/tests/test_cases/cousin_cases/' + doc_name, 'r') as cousinx_data:
            for line in cousinx_data:
                line = line.strip().split('\t')
                if len(line) != 1:
                    locus = line[0]
                    lr = line[3]
                    ref_dict[locus] = lr
                else:
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    elif line[0] == 'P':
                        p = float(line[1])

        return ref_dict, cpi, p

    # getting similar data from each patient's test data
    def get_test_data_cousinx(self, doc_name):
        test_dict = {}
        test_cpi = 1
        with open('cognation/scripts/tests/test_cases/cousin_cases/' + doc_name, 'r') as cousinx:
            for line in cousinx:
                cousin_formula_dict = self.calculate_relation(line)
                locus = cousin_formula_dict['locus']
                lr = float(cousin_formula_dict['lr'])
                test_cpi *= lr
                test_dict[locus] = lr

        test_cpi = int(test_cpi)
        test_p = (test_cpi / (1 + test_cpi)) * 100

        return test_dict, test_cpi, test_p

    # getting reference and test data for all parents
    for i in range(len(doc_names_list)):
        cousin_path = doc_names_list[i]
        overall_test_dict[cousin_path] = get_reference_data_cousinx(cousin_path)
        overall_ref_dict[cousin_path] = get_test_data_cousinx(cousin_path)

# Here we get all data for all test cases for grandparents
class GetGrandparentsData(GrandParentFormula):
    # getting loci and lrs in the dictionary and also CPI and P from each patient's reference data
    @staticmethod
    def get_reference_data_grandparentx(doc_name):
        ref_dict = {}
        cpi = 0
        p = 0.0
        with open('cognation/scripts/tests/test_cases/grandparent_cases/' + doc_name, 'r') as grandparentx_data:
            for line in grandparentx_data:
                line = line.strip().split('\t')
                if len(line) != 1:
                    locus = line[0]
                    lr = line[3]
                    ref_dict[locus] = lr
                else:
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    elif line[0] == 'P':
                        p = float(line[1])

        return ref_dict, cpi, p

    # getting similar data from each patient's test data
    def get_test_data_grandparentx(self, doc_name):
        test_dict = {}
        test_cpi = 1
        with open('cognation/scripts/tests/test_cases/grandparent_cases/' + doc_name, 'r') as grandparentx:
            for line in grandparentx:
                grandparent_formula_dict = self.calculate_relation(line)
                locus = grandparent_formula_dict['locus']
                lr = float(grandparent_formula_dict['lr'])
                test_cpi *= lr
                test_dict[locus] = lr

        test_cpi = int(test_cpi)
        test_p = (test_cpi / (1 + test_cpi)) * 100

        return test_dict, test_cpi, test_p

    # getting reference and test data for all parents
    for i in range(len(doc_names_list)):
        grandparent_path = doc_names_list[i]
        overall_test_dict[grandparent_path] = get_reference_data_grandparentx(grandparent_path)
        overall_ref_dict[grandparent_path] = get_test_data_grandparentx(grandparent_path)

# Here we get all data for all test cases for siblings
class GetSiblingsData(SiblingFormula):
    # getting loci and lrs in the dictionary and also CPI and P from each patient's reference data
    @staticmethod
    def get_reference_data_siblingx(doc_name):
        ref_dict = {}
        cpi = 0
        p = 0.0
        with open('cognation/scripts/tests/test_cases/siblings_cases/' + doc_name, 'r') as siblingsx:
            for line in siblingsx:
                line = line.strip().split('\t')
                if len(line) != 1:
                    locus = line[0]
                    lr = line[3]
                    ref_dict[locus] = lr
                else:
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    elif line[0] == 'P':
                        p = float(line[1])

        return ref_dict, cpi, p

    # getting similar data from each patient's test data
    def get_test_data_siblingx(self, doc_name):
        test_dict = {}
        test_cpi = 1
        with open('cognation/scripts/tests/test_cases/siblings_cases/' + doc_name, 'r') as siblingx:
            for line in siblingx:
                siblings_formula_dict = self.calculate_relation(line)
                locus = siblings_formula_dict['locus']
                lr = float(siblings_formula_dict['lr'])
                test_cpi *= lr
                test_dict[locus] = lr

        test_cpi = int(test_cpi)
        test_p = (test_cpi / (1 + test_cpi)) * 100

        return test_dict, test_cpi, test_p

    # getting reference and test data for all parents
    for i in range(len(doc_names_list)):
        siblings_path = doc_names_list[i]
        overall_test_dict[siblings_path] = get_reference_data_siblingx(siblings_path)
        overall_ref_dict[siblings_path] = get_test_data_siblingx(siblings_path)

# Here we get all data for all test cases for uncles
class GetUncleData(UncleFormula):
    # getting loci and lrs in the dictionary and also CPI and P from each patient's reference data
    @staticmethod
    def get_reference_data_unclex(doc_name):
        ref_dict = {}
        cpi = 0
        p = 0.0
        with open('cognation/scripts/tests/test_cases/uncle_cases/' + doc_name, 'r') as unclex:
            for line in unclex:
                line = line.strip().split('\t')
                if len(line) != 1:
                    locus = line[0]
                    lr = line[3]
                    ref_dict[locus] = lr
                else:
                    line = line[0].split('=')
                    if line[0] == 'CPI':
                        cpi = int(line[1])
                    elif line[0] == 'P':
                        p = float(line[1])

        return ref_dict, cpi, p

    # getting similar data from each patient's test data
    def get_test_data_unclex(self, doc_name):
        test_dict = {}
        test_cpi = 1
        with open('cognation/scripts/tests/test_cases/uncle_cases/' + doc_name, 'r') as unclex:
            for line in unclex:
                uncle_formula_dict = self.calculate_relation(line)
                locus = uncle_formula_dict['locus']
                lr = float(uncle_formula_dict['lr'])
                test_cpi *= lr
                test_dict[locus] = lr

        test_cpi = int(test_cpi)
        test_p = (test_cpi / (1 + test_cpi)) * 100

        return test_dict, test_cpi, test_p

    # getting reference and test data for all parents
    for i in range(len(doc_names_list)):
        uncle_path = doc_names_list[i]
        overall_test_dict[uncle_path] = get_reference_data_unclex(uncle_path)
        overall_ref_dict[uncle_path] = get_test_data_unclex(uncle_path)