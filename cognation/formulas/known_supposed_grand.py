from __future__ import unicode_literals
from .base import Formula, Calculations


class GrandKnownSupposed(Formula):
    def calculate_relation(self, raw_values):
        locus, alleles, sets, intersections, dict_make_result = self.getting_alleles_locus(raw_values, 4)
        child_alleles, grandparent_alleles = alleles[0], alleles[2]
        child_set, known_set, grandparent_set, supposed_set = sets
        target_alleles = [alleles[0], alleles[2]]
        ch_grand_inter = intersections[1]

        zeros_list = [0, 2, 5]
        for item in zeros_list:
            if len(intersections[item]) == 0:
                return self.make_result(locus, 0, dict_make_result)
        if self.is_gender_specific(locus):
            return self.result_gender_specific(locus, dict_make_result)
        c = Calculations()
        all_alleles = c.get_overall_alleles(target_alleles)
        freq_dict = self.get_frequencies(locus, all_alleles)

        supposed_genotypes = c.get_supposed_one_child(child_alleles, grandparent_alleles, known_set)
        if len(child_set) == 1:  # Homozygous child
            if type(supposed_genotypes) == list and supposed_set not in supposed_genotypes:
                return self.make_result(locus, 0, dict_make_result)  # checked if parent's combination is proper for counting lr, if not => lr=0
            lr = c.get_lr_from_possible(supposed_genotypes, freq_dict)
            return self.make_result(locus, 1 / lr, dict_make_result)

        if child_set != known_set:  # and heterozygous child
            for genotype in supposed_genotypes:
                if len(genotype) == 1:
                    lr = c.get_lr_from_possible(genotype, freq_dict)
                    return self.make_result(locus, 1 / lr, dict_make_result)
            if type(supposed_genotypes) == list and supposed_set not in supposed_genotypes:
                return self.make_result(locus, 0, dict_make_result)  # checked if parent's combination is proper for counting lr, if not => lr=0
            lr = c.get_lr_from_possible(supposed_genotypes, freq_dict)
            return self.make_result(locus, 1 / lr, dict_make_result)

        if len(ch_grand_inter) == 1 and len(grandparent_set) == 1:  # ab ab aa an
            lr = c.F(freq_dict[grandparent_alleles[0]])
            return self.make_result(locus, 1 / lr, dict_make_result)
        if grandparent_set == child_set:  # ab ab ab an/bn
            freq1, freq2 = freq_dict[child_alleles[0]], freq_dict[child_alleles[1]]
            lr = (freq1 + freq2) * (2 - (freq1 + freq2))
            return self.make_result(locus, 1 / lr, dict_make_result)
        if len(ch_grand_inter) == 1 and len(grandparent_set) == 2:  # ab ab ac an/bc
            freq1 = freq_dict[list(ch_grand_inter)[0]]
            freq2, freq3 = freq_dict[list(child_set - grandparent_set)[0]], freq_dict[list(grandparent_set - child_set)[0]]
            lr = c.F(freq1) + 2 * freq2 * freq3
            return self.make_result(locus, 1 / lr, dict_make_result)

        if supposed_set not in supposed_genotypes:
            return self.make_result(locus, 0, dict_make_result)
        lr = c.get_lr_from_possible(supposed_genotypes, freq_dict)
        return self.make_result(locus, 1 / lr, dict_make_result)
