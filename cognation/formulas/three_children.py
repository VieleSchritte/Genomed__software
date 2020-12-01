from __future__ import unicode_literals
from .base import Formula, Calculations
from .two_children import TwoChildrenFormula
from .parent import ParentFormula


class ThreeChildrenFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 4)
        parent_alleles, child1_alleles, child2_alleles, child3_alleles = alleles

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr = '-'
        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        if locus == 'AMEL':
            return self.make_result(locus, 1, dict_make_result)

        common_set = set(child1_alleles + child2_alleles + child3_alleles + parent_alleles)
        freq_dict = self.get_frequencies(locus, list(common_set))
        c = Calculations()
        lr = 0

        # If there are no intersections between children and parent, return lr = 0
        for i in range(0, 3):
            if intersections[i] == 0:
                return self.make_result(locus, lr, dict_make_result)

        print(locus)
        print('freqs: ', freq_dict)

        unique_genotype, repeat_genotype = c.get_repeat_unique(child1_alleles, child2_alleles, child3_alleles)
        # There are repeatable children genotypes
        if len(repeat_genotype) != 0:
            raw_values = [locus, '/'.join(parent_alleles), '/'.join(repeat_genotype)]

            # all children genotypes are same, use ParentFormula
            if len(unique_genotype) == 0:
                result = ParentFormula(Formula).calculate_relation(raw_values)
                lr = result['lr']
                return self.make_result(locus, lr, dict_make_result)

            # Two children have same genotypes
            else:
                raw_values.append('/'.join(unique_genotype))
                result = TwoChildrenFormula(Formula).calculate_relation(raw_values)
                lr = result['lr']
                return self.make_result(locus, lr, dict_make_result)

        # special cases (aa ab ac an) and (ab ac ad an)
        alleles_list = child1_alleles + child2_alleles + child3_alleles
        children_genotypes = alleles[1:]

        for allele in alleles_list:
            homo_counter = 0
            repeats_number = alleles_list.count(allele)
            for genotype in children_genotypes:
                homo_counter = genotype.count(allele)
                if homo_counter == 2:
                    break

            if repeats_number == 3 and homo_counter == 1 or repeats_number == 4 and homo_counter == 2:
                print('special case')
                freq = freq_dict[allele]
                lr = c.F(freq)
                return self.make_result(locus, lr, dict_make_result)

        children_alleles = child1_alleles + child2_alleles + child3_alleles
        children_set = set(children_alleles)

        print('default case through function')
        lr = self.lr_from_possible_genotypes(children_set, children_genotypes, freq_dict)
        print(lr)
        print()
        return self.make_result(locus, lr, dict_make_result)

    @staticmethod
    def lr_from_possible_genotypes(children_set, children_genotypes, freq_dict):
        all_children_alleles = list(children_set)
        possible_parents = []
        for i in range(len(all_children_alleles)):
            for j in range(len(all_children_alleles)):
                if j > i:
                    combination = [all_children_alleles[i], all_children_alleles[j]]
                    counter = 0
                    for allele in combination:
                        print('combination: ', combination)
                        for genotype in children_genotypes:
                            print('allele, genotype: ', allele, genotype)
                            print('allele not in genotype is ', allele not in genotype)
                            if allele not in genotype:
                                counter += 1
                        if counter < 2:
                            possible_parents.append(combination)
        print('possible parents: ', possible_parents)
        comb_sum = 0
        for combination in possible_parents:
            comb_sum += freq_dict[combination[0]] * freq_dict[combination[1]]
        lr = 2 * comb_sum
        return lr
