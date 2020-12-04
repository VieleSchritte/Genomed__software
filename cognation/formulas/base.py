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


# Abstract parent class
class Formula(abc.ABC):
    def __init__(self, user_data):
        self.user_data = user_data

    # Checking out if the locus is gender-specific (so we don't need to add it to cpi calculation)
    @staticmethod
    def is_gender_specific(locus):
        gender_specific_loci = ['SRY', 'DYS391', 'Yindel']
        for i in range(len(gender_specific_loci)):
            if locus == gender_specific_loci[i]:
                return True
        return False

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

    def calculate(self):
        result = OrderedDict()

        processed_user_data = []
        for participant in self.user_data:
            participant_lines = participant.splitlines()
            processed_participant = []
            for line in participant_lines:
                line = line.strip()
                processed_participant.append(re.split(r'[\s\t]+', line))
            processed_user_data.append(processed_participant)

        overall_participants = []
        for i in range(len(processed_user_data[0])):
            pair = []
            for j in range(len(processed_user_data)):
                target = processed_user_data[j][i]
                locus = target[0]
                if len(target) > 3:
                    locus += ' ' + target[1]
                    alleles = '/'.join(target[2:])
                else:
                    alleles = '/'.join(target[1:])
                pair.append([locus, alleles])
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
    def get_frequencies(self, locus, sat_set):
        result = {}
        for sat in sat_set:
            try:
                locus_object = Locus.objects.get(locus=locus, sat=self.normalize_sat(sat))
                result[sat] = locus_object.freq
            except Locus.DoesNotExist:
                raise UnknownAlleleException(locus, sat)

        return result

    @staticmethod
    def normalize_sat(value):
        # English and Russian variants
        if value == 'X' or value == 'Х':
            return 0.0
        elif value == 'Y':
            return 1.0
        else:
            return float(value)

    @staticmethod
    def split_sat(sat_string):
        return re.split(r'/', sat_string)

    def get_template(self):
        return 'cognation/formula/' + self.__class__.__name__.lower()[:-7] + '.html'

    @staticmethod
    def make_result(locus, lr, dict_alleles):
        if len(dict_alleles.keys()) == 2:
            return {
                "locus": locus,
                "part1": dict_alleles['part1'],
                "part2": dict_alleles['part2'],
                "lr": lr
            }
        elif len(dict_alleles.keys()) == 3:
            return {
                "locus": locus,
                "part1": dict_alleles['part1'],
                "part2": dict_alleles['part2'],
                "part3": dict_alleles['part3'],
                "lr": lr
            }
        elif len(dict_alleles.keys()) == 4:
            return {
                "locus": locus,
                "part1": dict_alleles['part1'],
                "part2": dict_alleles['part2'],
                "part3": dict_alleles['part3'],
                "part4": dict_alleles['part4'],
                "lr": lr
            }
        elif len(dict_alleles.keys()) == 5:
            return {
                "locus": locus,
                "part1": dict_alleles['part1'],
                "part2": dict_alleles['part2'],
                "part3": dict_alleles['part3'],
                "part4": dict_alleles['part4'],
                "part5": dict_alleles['part5'],
                "lr": lr
            }

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
