from cognation.models import Locus
import abc
import re
from collections import OrderedDict


class UnknownFormulaException(Exception):
    def __init__(self, formula_type):
        self.formula_type = formula_type


class LineFormatException(Exception):
    def __str__(self):
        return 'Wrong line format'


class AllelesException(Exception):
    def __init__(self, locus, part):
        self.locus = locus
        self.part = part

    def __str__(self):
        return "Alleles count doesn't look right: " + str(self.part) + " in locus " + str(self.locus)


class UnknownAlleleException(Exception):
    def __init__(self, locus, sat):
        self.locus = locus
        self.sat = sat

    def __str__(self):
        return "Unknown allele found in locus " + str(self.locus) + ": " + str(self.sat)


class UnknownSymbolInAlleles(Exception):
    def __init__(self, locus, alleles, symbol):
        self.locus = locus
        self.alleles = alleles
        self.symbol = symbol

    def __str__(self):
        return "Unknown symbol '" + str(self.symbol) + "' found in alleles " + str('/'.join(self.alleles)) + " of locus " + str(self.locus)


class TooManyDelimitingSymbols(Exception):
    def __init__(self, locus, alleles):
        self.locus = locus
        self.alleles = alleles

    def __str__(self):
        return "Too many delimiting symbols found in alleles " + str(self.alleles) + " of locus " + str(self.locus) + ". Use only one '.' symbol in case of float number"


class DelimitingLast(Exception):
    def __init__(self, alleles):
        self.alleles = alleles

    def __str__(self):
        return "Delimiting character in the end of alleles line: " + str(self.alleles)


class LociSetDoesNotEqual(Exception):
    def __str__(self):
        return "participants' loci sets are not match"


# Abstract parent class
class Formula(abc.ABC):
    def __init__(self, user_data):
        self.user_data = user_data

    # Checking out if the locus is gender-specific (so we don't need to add it to cpi calculation)
    @staticmethod
    def is_gender_specific(locus):
        gender_specific_loci = ['SRY', 'DYS391', 'Yindel', 'AMEL']
        for i in range(len(gender_specific_loci)):
            if locus == gender_specific_loci[i]:
                return True
        return False

    def preparation_check(self, locus, dict_make_result):
        if locus == 'AMEL':
            return self.make_result(locus, 1, dict_make_result)
        return self.make_result(locus, '-', dict_make_result)

    def getting_alleles_locus(self, raw_values, part_number):
        if len(raw_values) < part_number + 1:
            # Skip line with warning
            raise LineFormatException()

        part_alleles = []

        dict_make_result = {}
        locus = raw_values[0]

        if len(raw_values) > part_number + 1:
            locus += ' ' + raw_values[1]
            for i in range(2, part_number + 2):
                part_alleles.append(self.split_sat(raw_values[i]))
                key = 'part' + str(i-1)
                dict_make_result[key] = '/'.join(part_alleles[i - 2])
        else:
            for i in range(1, part_number + 1):
                part_alleles.append(self.split_sat(raw_values[i]))
                key = 'part' + str(i)
                dict_make_result[key] = '/'.join(part_alleles[i - 1])

        part_sets = []
        for part in part_alleles:
            if not self.is_gender_specific(locus) and len(part) != 2:
                raise AllelesException(locus, part)

            part_sets.append(set(part))

        if not self.is_gender_specific(locus):
            self.get_frequencies(locus, part_alleles[0] + part_alleles[1])

        intersections = []
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

    def alleles_check(self, alleles, locus):
        for allele in alleles:
            if locus == 'AMEL':
                allowed_alleles = ['X', 'Y']
                if allele not in allowed_alleles:
                    counter = 0
                    for symbol in allele:
                        if symbol not in allowed_alleles:
                            raise UnknownSymbolInAlleles(locus, alleles, symbol)
                        else:
                            counter += 1
                    if counter > 1:
                        raise LineFormatException

            else:
                exception_counter = 0
                if not self.is_digit(allele):
                    for symbol in allele:
                        if not symbol.isdigit() and symbol != '.':
                            raise UnknownSymbolInAlleles(locus, alleles, symbol)
                        if symbol == '.':
                            exception_counter += 1
                    if exception_counter > 1:
                        raise TooManyDelimitingSymbols(locus, alleles)
            if allele[-1] == '.':
                raise DelimitingLast(alleles)

    def calculate(self):
        result = OrderedDict()
        processed_user_data = []
        for participant in self.user_data:
            participant_lines = sorted(participant.splitlines())
            participant = []
            for line in participant_lines:
                # Skip empty line
                if len(line) == 0:
                    continue
                line = re.split(r'[\s\t]+', line.strip())
                if len(line) == 1:
                    raise LineFormatException

                locus, alleles = line[0], []
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
                if not self.is_gender_specific(locus) and '/' not in alleles:
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
            except (LineFormatException, AllelesException, UnknownAlleleException) as exception:
                result[hash(str(line))] = {'exception': exception, 'line': str(line)}
        return result

    # getting allele frequencies from DB
    @staticmethod
    def get_frequencies(locus, sat_set):
        result = {}
        for sat in sat_set:
            try:
                locus_object = Locus.objects.get(locus=locus, sat=float(sat))
                result[sat] = locus_object.freq
            except Locus.DoesNotExist:
                raise UnknownAlleleException(locus, sat)
        return result

    @staticmethod
    def split_sat(sat_string):
        return re.split(r'/', sat_string)

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
        if len(key_set) == 1:
            refutation = c.homo_refutation(freq1)
        else:
            refutation = c.hetero_refutation(freq1, freq2)
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
    def get_repeat_unique(child1_alleles, child2_alleles, child3_alleles):
        children_genotypes = [child1_alleles, child2_alleles, child3_alleles]
        unique_genotype, repeatable_genotype = [], []

        for i in range(len(children_genotypes)):
            for j in range(len(children_genotypes)):
                for k in range(len(children_genotypes)):

                    if j > i and children_genotypes[i] == children_genotypes[j]:
                        if children_genotypes[k] != children_genotypes[i]:
                            unique_genotype = children_genotypes[k]
                            repeatable_genotype = children_genotypes[i]

        if children_genotypes[0] == children_genotypes[1] == children_genotypes[2]:
            repeatable_genotype = children_genotypes[0]

        return unique_genotype, repeatable_genotype
