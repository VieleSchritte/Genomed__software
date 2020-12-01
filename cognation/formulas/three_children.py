from __future__ import unicode_literals
from .base import Formula, Calculations
from .two_children import TwoChildrenFormula
from .parent import ParentFormula


class ThreeChildrenFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 4)
        child3_alleles, child2_alleles, child1_alleles, parent_alleles, = alleles
        child3_set, child2_set, child1_set, parent_set = sets
        ch2ch3_inter, ch1ch3_inter, ch3p_inter, ch1ch2_inter = intersections[0:4]

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
        if len(intersections[2]) == 0:
            for i in range(4, 6):
                if len(intersections[i]) == 0:
                    return self.make_result(locus, lr, dict_make_result)

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
        if ch1ch2_inter == ch2ch3_inter == ch3p_inter:
            freq = freq_dict[list(ch1ch2_inter)[0]]
            lr = c.F(freq)
            return self.make_result(locus, lr, dict_make_result)

        #  Homozygous 1st child
        if len(child1_set) == 1:
            # aa ab bc ab/ac
            if len(ch1ch2_inter) != 0 and len(child2_set) == len(child3_set) == len(parent_set) == 2:
                freq1 = freq_dict[child1_alleles[0]]
                freq2, freq3 = freq_dict[child3_alleles[0]], freq_dict[child3_alleles[1]]
                lr = 2 * freq1 * (freq2 + freq3)
                return self.make_result(locus, lr, dict_make_result)

            # default is lr = 2 * P(parent_alleles[0]) * P(parent_alleles[1])
            freq1, freq2 = freq_dict[parent_alleles[0]], freq_dict[parent_alleles[1]]
            lr = 2 * freq1 * freq2
            return self.make_result(locus, lr, dict_make_result)

        # Heterozygous 1st child
        if len(child1_set) == 2:
            # ab ac bc ab/ac/bc
            if len(common_set) == 3:
                freq1, freq2 = freq_dict[child1_alleles[0]], freq_dict[child1_alleles[1]]
                freq3 = freq_dict[list(child2_set - child1_set)[0]]
                lr = 2 * (freq1 * freq2 + freq1 * freq3 + freq2 * freq3)
                return self.make_result(locus, lr, dict_make_result)

            if len(ch1ch2_inter) == 1:
                # ab ac bd ad/bc
                if len(common_set) == 4 and len(ch2ch3_inter) == 0:
                    freq1, freq2 = freq_dict[list(ch1ch2_inter)[0]], freq_dict[list(ch1ch3_inter)[0]]
                    freq3, freq4 = freq_dict[list(child2_set - child1_set)[0]], freq_dict[list(child3_set - child1_set)[0]]
                    lr = 2 * (freq1 * freq4 + freq2 * freq3)
                    return self.make_result(locus, lr, dict_make_result)

                # ab ac cd ac/ad/bc
                if len(common_set) == 4 and len(ch1ch2_inter) == len(ch2ch3_inter) == 1:
                    freq1, freq2 = freq_dict[list(ch1ch2_inter)[0]], freq_dict[list(child1_set - child2_set)[0]]
                    freq3, freq4 = freq_dict[list(ch2ch3_inter)[0]], freq_dict[list(child3_set - child2_set)[0]]
                    lr = 2 * (freq1 * freq3 + freq1 * freq4 + freq2 * freq3)
                    return self.make_result(locus, lr, dict_make_result)

                # ab ac de ad/ae
                if len(common_set) == 5 and len(ch2ch3_inter) == 0:
                    freq1 = freq_dict[list(ch1ch2_inter)[0]]
                    freq2, freq3 = freq_dict[child3_alleles[0]], freq_dict[child3_alleles[1]]
                    lr = 2 * freq1 * (freq2 + freq3)
                    return self.make_result(locus, lr, dict_make_result)

            # ab cd cn (n != a, b, d), ac/bc
            if len(ch1ch2_inter) == 0:
                freq1, freq2 = freq_dict[child1_alleles[0]], freq_dict[child1_alleles[1]]
                freq3 = freq_dict[list(ch2ch3_inter)[0]]
                lr = 2 * freq3 * (freq1 + freq2)
                return self.make_result(locus, lr, dict_make_result)
