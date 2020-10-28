from __future__ import unicode_literals
from .base import Formula, Calculations


class ParentFormula(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 2)
        child1_alleles, child2_alleles, child3_alleles, parent_alleles = alleles
        child1_set, child2_set, child3_set, parent_set = sets
        ch1ch2_inter, ch1ch3_inter, ch1p_inter, ch2ch3_inter, ch2p_inter, ch2ch3_inter, ch3p_inter = intersections

        c = Calculations()

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr = '-'
        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        if ch1p_inter == 0 or ch2p_inter == 0 or ch3p_inter == 0:
            return self.make_result(locus, 0, dict_make_result)

        common_alleles = child1_alleles + child2_alleles + child3_alleles + parent_alleles
        common_set = set(common_alleles)
        freq_dict = self.get_frequencies(locus, list(common_set))

        #  Homozygous 1st child
        if len(child1_set) == 1:
            if len(ch1ch2_inter) == len(ch2ch3_inter) == len(ch3p_inter) == 1:
                freq = freq_dict[list(ch1ch2_inter)[0]]
                lr = c.F(freq)
                return self.make_result(locus, lr, dict_make_result)

            # Heterozygous 2nd child with intersection with 1st child
            if len(child2_set) == 2 and len(ch1ch2_inter) != 0:

                # aa ab bb ab
                if len(child3_set) == 1 and len(common_set) == 2:
                    freq1, freq2 = freq_dict[list(common_set)[0]], freq_dict[list(common_set)[1]]
                    lr = 2 * freq1 * freq2
                    return self.make_result(locus, lr, dict_make_result)

                # aa ab bc ab/ac
                if len(child3_set) == 2 and len(common_set) == 3:
                    freq1 = freq_dict[list(ch1ch2_inter)[0]]
                    freq2, freq3 = freq_dict[child3_alleles[0]], freq_dict[child3_alleles[1]]
                    lr = 2 * freq1 * (freq2 + freq3)
                    return self.make_result(locus, lr, dict_make_result)

                # aa ab cc ac
                if len(child3_set) == 1 and len(common_set) == 3:
                    freq1, freq2 = freq_dict[child1_alleles[0]], freq_dict[child3_alleles[0]]
                    lr = 2 * freq1 * freq2
                    return self.make_result(locus, lr, dict_make_result)
