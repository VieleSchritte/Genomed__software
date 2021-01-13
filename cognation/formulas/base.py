from cognation.models import Locus
import abc
import re
from collections import OrderedDict


class UnknownFormulaException(Exception):
    def __init__(self, formula_type):
        self.formula_type = formula_type


class LineFormatException(Exception):
    def __init__(self, line):
        self.line = line

    def __str__(self):
        return_string = 'Неверный формат ввода: '
        if type(self.line) == list:
            for item in self.line:
                return_string += str(item) + ' '
            return return_string
        if type(self.line) == str:
            return return_string + ' ' + self.line


class AllelesException(Exception):
    def __init__(self, locus, part):
        self.locus = locus
        self.part = part

    def __str__(self):
        return "Неверное число аллелей: " + str(self.part) + " в локусе " + str(self.locus)


class UnknownAlleleException(Exception):
    def __init__(self, locus, sat):
        self.locus = locus
        self.sat = sat

    def __str__(self):
        return "В локусе " + str(self.locus) + " найден неизвестный аллель: " + str(self.sat)


class UnknownSymbolInAlleles(Exception):
    def __init__(self, locus, alleles, symbol):
        self.locus = locus
        self.alleles = alleles
        self.symbol = symbol

    def __str__(self):
        return "Неизвестный символ '" + str(self.symbol) + "' найден среди аллелей " + str(
            '/'.join(self.alleles)) + " в локусе " + str(self.locus)


class TooManyDelimitingSymbols(Exception):
    def __init__(self, locus, alleles):
        self.locus = locus
        self.alleles = alleles

    def __str__(self):
        return "Слишком много разделяющих символов в списке аллелей " + str(self.alleles) + " в локусе " + str(
            self.locus) + ". Используйте только одну точку или запятую"


class DelimitingLast(Exception):
    def __init__(self, alleles, locus):
        self.alleles = alleles
        self.locus = locus

    def __str__(self):
        return "Разделяющий символ в конце числа в локусе " + str(self.locus) + ": " + str(self.alleles)


class DelimitingFirst(Exception):
    def __init__(self, alleles, locus):
        self.alleles = alleles
        self.locus = locus

    def __str__(self):
        return "Разделяющий символ в начале числа в локусе : " + str(self.locus) + ": " + str(self.alleles)


class LociSetDoesNotEqual(Exception):
    def __str__(self):
        return "Числа локусов участников не совпадают. Введите одинаковое количество локусов для каждого участника"


class UnknownLocus(Exception):
    def __init__(self, locus):
        self.locus = locus

    def __str__(self):
        return "Введен неизвестный локус: " + str(self.locus) + ". Проверьте правильность вводимых данных."


class EmptyAlleles(Exception):
    def __init__(self, locus):
        self.locus = locus

    def __str__(self):
        return "Отсутствуют аллели в локусе " + str(self.locus) + ". Проверьте правильность вводимых данных."


# Abstract parent class
class Formula(abc.ABC):
    def __init__(self, user_data):
        self.user_data = user_data

    # Checking out if the locus is gender-specific (so we don't need to add it to cpi calculation)
    @staticmethod
    def is_gender_specific(locus):
        gender_specific_loci = ['SRY', 'DYS391', 'Yindel', 'AMEL']
        for g_locus in gender_specific_loci:
            if locus == g_locus:
                return True
        return False

    def result_gender_specific(self, locus, dict_make_result):
        if locus == 'AMEL':
            return self.make_result(locus, 1, dict_make_result)
        return self.make_result(locus, '-', dict_make_result)

    def getting_alleles_locus(self, raw_values, part_number):
        if len(raw_values) < part_number + 1:  # Skip line with warning
            raise LineFormatException(raw_values)

        part_alleles, dict_make_result, locus = [], {}, raw_values[0]
        if len(raw_values) > part_number + 1:
            locus += ' ' + raw_values[1]
            for i in range(2, part_number + 2):
                part_alleles.append(self.split_alleles(raw_values[i]))
                dict_make_result['part' + str(i - 1)] = '/'.join(part_alleles[i - 2])
        else:
            for i in range(1, part_number + 1):
                part_alleles.append(self.split_alleles(raw_values[i]))
                dict_make_result['part' + str(i)] = '/'.join(part_alleles[i - 1])

        part_sets, intersections = [], []
        for part in part_alleles:
            if not self.is_gender_specific(locus) and len(part) != 2:
                raise AllelesException(locus, part)
            part_sets.append(set(part))

        if not self.is_gender_specific(locus):
            self.get_frequencies(locus, part_alleles[0] + part_alleles[1])  # to see exception if needed
        for i in range(len(part_alleles)):
            for j in range(len(part_alleles)):
                if j > i:
                    intersections.append(part_sets[i] & part_sets[j])

        return locus, part_alleles, part_sets, intersections, dict_make_result

    # Checks if the string is integer or float
    @staticmethod
    def is_digit(string):
        if string.isdigit():
            return True
        else:
            try:
                float(string)
                return True
            except ValueError:
                return False

    def alleles_check(self, alleles, locus):  # Checks if everything's right in alleles, if not - raises correct exception
        for allele in alleles:
            if locus == 'AMEL':
                allowed_alleles = ['X', 'Х', 'Y']  # Both English and Russian X variants
                if allele not in allowed_alleles:
                    counter = 0
                    for symbol in allele:
                        if symbol not in allowed_alleles:
                            raise UnknownSymbolInAlleles(locus, alleles, symbol)
                        else:
                            counter += 1
                    if counter > 1:
                        raise LineFormatException(alleles)
            elif not self.is_gender_specific(locus):
                exception_counter = 0
                if not self.is_digit(allele):
                    for symbol in allele:
                        if not symbol.isdigit() and symbol != '.':
                            raise UnknownSymbolInAlleles(locus, alleles, symbol)
                        if symbol == '.':
                            exception_counter += 1
                    if exception_counter > 1:
                        raise TooManyDelimitingSymbols(locus, alleles)
            elif self.is_gender_specific(locus):
                return

            if allele[-1] == '.':
                raise DelimitingLast(allele, locus)
            if allele[0] == '.':
                raise DelimitingFirst(allele, locus)
        if len(alleles) > 2:
            raise AllelesException(locus, alleles)

    @staticmethod
    def locus_check(locus):
        possible_loci = ['AMEL', 'D3S1358', 'vWA', 'D16S539', 'CSF1PO', 'TPOX', 'D8S1179', 'D21S11', 'SE33', 'D18S51',
                         'Penta E', 'D2S441', 'D19S433', 'TH01', 'FGA', 'D22S1045', 'D5S818', 'D13S317', 'D7S820',
                         'D6S1043', 'D10S1248', 'D1S1656', 'D12S391', 'D2S1338', 'Penta D', 'Yindel', 'DYS391', 'SRY']
        if locus not in possible_loci:
            raise UnknownLocus(locus)

    @staticmethod
    def is_locus(string):
        possible_loci = ['AMEL', 'D3S1358', 'vWA', 'D16S539', 'CSF1PO', 'TPOX', 'D8S1179', 'D21S11', 'SE33', 'D18S51',
                         'Penta E', 'D2S441', 'D19S433', 'TH01', 'FGA', 'D22S1045', 'D5S818', 'D13S317', 'D7S820',
                         'D6S1043', 'D10S1248', 'D1S1656', 'D12S391', 'D2S1338', 'Penta D', 'Yindel', 'DYS391', 'SRY']
        if string in possible_loci:
            return True
        return False

    def calculate(self):
        result = OrderedDict()
        processed_user_data = []
        for participant in self.user_data:
            participant_lines = sorted(participant.splitlines())
            participant = []
            for line in participant_lines:
                if len(line) == 0:  # Skip empty line
                    continue
                line = re.split(r'[\s\t]+', line.strip())
                print(line)
                if len(line) == 1:
                    if self.is_locus(line[0]):
                        print('should raise empty alleles!')
                        raise EmptyAlleles(line[0])
                    print('went to line format')
                    raise LineFormatException(line)

                locus, alleles = line[0], []
                self.locus_check(locus)
                if locus == 'AMEL':
                    alleles = line[1:]
                else:
                    for i in range(1, len(line)):
                        if line[i].isalpha():
                            locus += ' ' + line[i]
                        else:
                            if ',' in line[i]:
                                alleles.append(line[i].replace(',', '.'))
                            else:
                                alleles.append(line[i])
                self.alleles_check(alleles, locus)
                alleles = '/'.join(alleles)
                homozygous_cases = [
                    not self.is_gender_specific(locus) and '/' not in alleles,
                    locus == 'AMEL' and '/' not in alleles
                ]
                for condition in homozygous_cases:
                    if condition:
                        alleles += '/' + alleles
                participant.append([locus, alleles])
            processed_user_data.append(participant)

        loci_numbers = []
        for participant in processed_user_data:
            loci_numbers.append(len(participant))
        for i in range(len(loci_numbers)):
            for j in range(len(loci_numbers)):
                if j > i:
                    if loci_numbers[i] != loci_numbers[j]:
                        raise LociSetDoesNotEqual

        overall_participants = []
        for i in range(len(processed_user_data[0])):
            pair = []
            for j in range(len(processed_user_data)):
                target = processed_user_data[j][i]
                pair.append([target[0], target[1]])
            base = pair[0]
            for k in range(1, len(pair)):
                base.append(pair[k][1])
            overall_participants.append(base)

        for line in overall_participants:
            if len(line) == 0:
                continue
            try:
                relation = self.calculate_relation(line)
                result[relation['locus']] = relation
            except (LineFormatException, AllelesException, UnknownAlleleException, UnknownSymbolInAlleles,
                    TooManyDelimitingSymbols, DelimitingLast, DelimitingFirst, LociSetDoesNotEqual) as exception:
                result[hash(str(line))] = {'exception': exception, 'line': str(line)}
        return result

    # getting allele frequencies from DB
    @staticmethod
    def get_frequencies(locus, alleles_list):
        result = {}
        for allele in alleles_list:
            try:
                locus_object = Locus.objects.get(locus=locus, sat=float(allele))
                result[allele] = locus_object.freq
            except Locus.DoesNotExist:
                raise UnknownAlleleException(locus, allele)
        return result

    @staticmethod
    def split_alleles(alleles_string):
        return re.split(r'/', alleles_string)

    def get_template(self):
        return 'cognation/formula/' + self.__class__.__name__.lower()[:-7] + '.html'

    @staticmethod
    def make_result(locus, lr, dict_make_result):
        dict_make_result["locus"], dict_make_result["lr"] = locus, lr
        return dict_make_result

    def get_division_lr(self, locus, key_set, alleles_list, confirmation):
        c = Calculations()
        freq_dict = self.get_frequencies(locus, alleles_list)
        freq1, freq2 = freq_dict[alleles_list[0]], freq_dict[alleles_list[1]]
        answers = {
            len(key_set) == 1: c.homo_refutation(freq1),
            len(key_set) != 1: c.hetero_refutation(freq1, freq2)
        }
        refutation = c.get_lr_from_cond_dict_short(answers)
        return confirmation / refutation

    # Abstract methods
    @abc.abstractmethod
    def calculate_relation(self, row_values):
        """ Abstract method to calculate """


# Calculation Helpers
class Calculations:
    #  A helper for the frequently used pattern F(Px) = Px * (2 - Px)
    @staticmethod
    def F(freq):
        return freq * (2 - freq)

    #  A helper for the frequently used pattern Q(Px) = 0.5 - 0.5 * Px (for GrandParentFormula)
    @staticmethod
    def Q(freq):
        return 0.5 + 0.5 * freq

    #  A helper for the frequently used pattern M(Px, Py) = 2 * Px * Py / F(Px) (for SiblingFormula)
    def M(self, freq1, freq2):
        return 2 * freq1 * freq2 / self.F(freq1)

    #  Probability of relation theory refutation in case of inspected person's homozygosity
    def homo_refutation(self, freq):
        return (self.F(freq)) ** 2

    #  Probability of relation theory refutation in case of inspected person's heterozygosity
    def hetero_refutation(self, freq1, freq2):
        return 2 * self.F(freq1) * self.F(freq2) - (2 * freq1 * freq2) ** 2

    # Returns unique and repeatable children's genotypes in case of two and three children
    @staticmethod
    def get_repeat_unique(children_genotypes):
        repeats_dict = {}
        for genotype in children_genotypes:
            repeats_dict[tuple(genotype)] = children_genotypes.count(genotype)
        return repeats_dict

    def get_repeatable_lr(self, raw_values, children_genotypes, formulas):
        repeats_dict = self.get_repeat_unique(children_genotypes)
        for key in repeats_dict:
            raw_values.append('/'.join(key))
        pos_dict = {
            2: False,
            1: formulas[0].calculate_relation
        }
        if len(children_genotypes) == 3:
            pos_dict = {
                3: False,
                2: formulas[1].calculate_relation,
                1: formulas[0].calculate_relation
            }
        for key in pos_dict.keys():
            if len(repeats_dict.keys()) == key:
                if not pos_dict[key]:
                    return pos_dict[key]
                else:
                    return pos_dict[key](raw_values)['lr']

    def get_possible_genotypes(self, children_alleles, children_genotypes, parents_data):
        single_alleles_combinations = []
        parent_set, key_word = parents_data[0], parents_data[1]
        for children_allele in children_alleles:
            single_alleles_combinations.append(self.get_combinations(list(parent_set), [children_allele]))
        if self.is_get_F_case(single_alleles_combinations, children_genotypes, parents_data, children_alleles):
            return set(children_genotypes[0]) & set(children_genotypes[1])  # case lr = F(Pa)

        possible_parent_genotypes = self.get_possible_parent_genotypes(children_alleles)
        answer = []
        if key_word == 'known':
            for parent_genotype in possible_parent_genotypes:
                known_supposed_alleles = self.get_overall_alleles([list(parent_genotype), list(parent_set)])
                if len(known_supposed_alleles) < len(children_alleles):
                    continue
                possible_children_genotypes = self.get_combinations(list(parent_genotype), list(parent_set))
                answer = self.answer_genotypes_selection(key_word, children_genotypes, possible_children_genotypes, answer, parent_genotype)
            return answer
        if key_word == 'supposed':
            for possible_parent in possible_parent_genotypes:
                answer = self.answer_genotypes_selection(key_word, children_genotypes, possible_parent_genotypes, answer, possible_parent)
            return answer

    @staticmethod
    def answer_genotypes_selection(key_word, children_genotypes, possible_genotypes, answer, parent_genotype):
        if key_word == 'known':
            counter = 0
            for child_genotype in children_genotypes:
                counter += possible_genotypes.count(set(child_genotype))
            if counter >= len(children_genotypes):
                answer.append(parent_genotype)
            return answer
        if key_word == 'supposed':
            counter = 0
            for child_genotype in children_genotypes:
                if len(set(child_genotype) & parent_genotype) != 0:
                    counter += 1
            if counter == len(children_genotypes):
                answer.append(parent_genotype)
            return answer

    @staticmethod
    def F_check(single_alleles_combinations, children_genotypes):
        for item in single_alleles_combinations:
            counter = 0
            for child_genotype in children_genotypes:
                counter += item.count(set(child_genotype))
            if counter == len(children_genotypes):
                return True

    def is_get_F_case(self, single_alleles_combinations, children_genotypes, parents_data, children_alleles):
        if self.F_check(single_alleles_combinations, children_genotypes):
            return True
        elif parents_data[1] == 'supposed':
            another_combinations = []
            for parent_allele in list(parents_data[0]):
                another_combinations.append(self.get_combinations(children_alleles, [parent_allele]))
            if self.F_check(another_combinations, children_genotypes):
                return True

    @staticmethod
    def get_possible_parent_genotypes(children_alleles):
        possible_parent_genotypes = []
        for i in range(len(children_alleles)):
            for j in range(len(children_alleles)):
                if j > i:
                    possible_parent_genotypes.append({children_alleles[i], children_alleles[j]})
        return possible_parent_genotypes

    @staticmethod
    def get_combinations(known_alleles, children_allele):
        combinations = []
        for other_allele in children_allele:
            for children_allele in known_alleles:
                combinations.append({children_allele, other_allele})
        return combinations

    @staticmethod
    def combination_processing(combinations):
        new_combinations = []
        for combination in combinations:
            if combinations.count(combination) == len(combinations):
                return [combination]
            if combination not in new_combinations:
                new_combinations.append(combination)
        return new_combinations

    # function for getting supposed parent's possible genotypes using known genotype and grandparent's one
    def get_supposed_one_child(self, child_alleles, grandparent_alleles, known_set):
        combinations = self.combination_processing(self.get_combinations(grandparent_alleles, child_alleles))
        if len(set(child_alleles)) == 1:  # Homozygous child
            len_one_combination = self.is_len_one(combinations)
            if len_one_combination:  # case of lr = F(freq)
                return set(len_one_combination)
            else:
                return combinations  # other cases can be processed through function for getting lr from possible genotypes

        if known_set != set(child_alleles):
            parent_ch_allele, new_combinations = list(set(child_alleles) - known_set)[0], []
            for combination in combinations:
                if parent_ch_allele in combination:
                    new_combinations.append(combination)
            return new_combinations
        return combinations

    @staticmethod
    def is_len_one(combinations):
        for combination in combinations:
            if len(combination) == 1:
                return combination

    @staticmethod
    def get_overall_alleles(genotypes):
        overall_alleles = []
        for genotype in genotypes:
            overall_alleles += genotype
        overall_alleles = list(set(overall_alleles))
        return overall_alleles

    def get_lr_from_possible(self, possible_parent_genotypes, freq_dict):
        lr = 0
        if type(possible_parent_genotypes) == set:  # cases where lr = c.F(Pa)
            lr = self.F(freq_dict[list(possible_parent_genotypes)[0]])
            return lr
        for genotype in possible_parent_genotypes:
            genotype = list(genotype)
            freq1, freq2 = freq_dict[genotype[0]], freq_dict[genotype[1]]
            lr += 2 * freq1 * freq2
        return lr

    @staticmethod
    def multiply_lr_on_children_allele(lr, children_alleles, freq_dict):
        for allele in children_alleles:
            lr *= freq_dict[allele]
        return lr

    @staticmethod
    def homo_counter(target_sets):
        homo_counter = 0
        for single_set in target_sets:
            if len(single_set) == 1:
                homo_counter += 1
        return homo_counter

    @staticmethod
    def hetero_counter(target_sets):
        hetero_counter = 0
        for child_set in target_sets:
            if len(child_set) == 2:
                hetero_counter += 1
        return hetero_counter

    @staticmethod
    def get_correct_frequency_order_couple(children_sets, freq_dict, children_alleles):
        freq1, other_alleles = 1, []
        for child_set in children_sets:
            if len(child_set) == 1:
                freq1, other_alleles = freq_dict[list(child_set)[0]], list(set(children_alleles) - child_set)
        if len(other_alleles) == 0:
            target_inter_list = list(children_sets[0] & children_sets[1])
            freq1, other_alleles = freq_dict[target_inter_list[0]],  list(set(children_alleles) - set(target_inter_list))
        if len(other_alleles) == 1:
            return freq1, freq_dict[other_alleles[0]], 0
        return freq1, freq_dict[other_alleles[0]], freq_dict[other_alleles[1]]

    @staticmethod
    def get_lr_from_dict_couple(overall_dict, hetero_counter, children_alleles_len):
        for key in overall_dict.keys():
            if key == hetero_counter:
                types = [int, float]
                if type(overall_dict[hetero_counter]) in types:
                    return overall_dict[hetero_counter]
                target_dict = overall_dict[hetero_counter]
                for length in target_dict.keys():
                    if length == children_alleles_len:
                        return target_dict[length]

    @staticmethod
    def get_lr_from_cond_dict_long(answers):
        counter = 0
        for condition in answers.keys():
            if not condition:
                counter += 1
            if condition:
                return 1 / answers[condition]
        if counter == len(answers.keys()):
            return 0

    @staticmethod
    def get_lr_from_cond_dict_short(answers):
        for key in answers.keys():
            if key:
                return answers[key]
