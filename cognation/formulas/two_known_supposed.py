from __future__ import unicode_literals
from .base import Formula, Calculations
from .one_known_supposed import OneKnownSupposedFormula


class TwoKnownSupposedFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, alleles, sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 4)
        # known child1 child2 supposed => supposed child2 child1 known
        supposed_alleles, child2_alleles, child1_alleles, known_alleles = alleles
        supposed_set, child2_set, child1_set, known_set = sets
        sch2_inter, sch1_inter, sk_inter, ch1ch2_inter, kch2_inter, kch1_inter = intersections

        print(locus)

        if self.is_gender_specific(locus):
            return self.make_result(locus, '-', dict_make_result)

        # If there are no intersections between children and parents, return lr = 0 and start counting mutations
        for i in range(1, 5):
            if i == 2 or i == 3:
                continue
            if len(intersections[i]) == 0:
                return self.make_result(locus, 0, dict_make_result)

        common_set = set(child1_alleles + child2_alleles + known_alleles + supposed_alleles)
        freq_dict = self.get_frequencies(locus, list(common_set))
        c = Calculations()

        # If children's genotypes are same, use OneKnownSupposedFormula
        if child1_set == child2_set:
            raw_values = [locus, '/'.join(known_alleles), '/'.join(child1_alleles), '/'.join(supposed_alleles)]
            result = OneKnownSupposedFormula(Formula).calculate_relation(raw_values)
            lr = result['lr']
            return self.make_result(locus, lr, dict_make_result)

        # homozygous first child
        if len(child1_set) == 1:
            # aa ab ab an
            if (len(common_set) == 2 or len(common_set) == 3) and child2_set == known_set:
                freq = freq_dict[child1_alleles[0]]
                lr = c.F(freq)
                return self.make_result(locus, lr, dict_make_result)

            freq1, freq2 = freq_dict[supposed_alleles[0]], freq_dict[supposed_alleles[1]]
            lr = 2 * freq1 * freq2
            return self.make_result(locus, lr, dict_make_result)

        # Heterozygous first child
        else:
            # ab ac ab ac/bc
            unique_ch2_allele = list(child2_set - child1_set)[0]
            unique_ch1_allele = list(child1_set - child2_set)[0]
            if child1_set == known_set and unique_ch2_allele in supposed_alleles and unique_ch2_allele not in child1_set:
                freq1, freq2, freq3 = freq_dict[child1_alleles[0]], freq_dict[child1_alleles[1]], freq_dict[list(child2_set - child1_set)[0]]
                lr = 2 * freq3 * (freq1 + freq2)
                return self.make_result(locus, lr, dict_make_result)

            # ab ac an (n != b) bc
            if ch1ch2_inter == kch2_inter and unique_ch1_allele not in known_alleles:
                freq1, freq2 = freq_dict[supposed_alleles[0]], freq_dict[supposed_alleles[1]]
                lr = 2 * freq1 * freq2
                return self.make_result(locus, lr, dict_make_result)

            # ab ac bc an
            if list(ch1ch2_inter)[0] in supposed_alleles and list(ch1ch2_inter)[0] not in known_alleles:
                lr = c.F(list(ch1ch2_inter)[0])
                return self.make_result(locus, lr, dict_make_result)

            # ab cd bd ac
            if len(ch1ch2_inter) == 0 and len(sk_inter) == 0:
                freq1, freq2 = freq_dict[supposed_alleles[0]], freq_dict[supposed_alleles[1]]
                lr = 2 * freq1 * freq2
                return self.make_result(locus, lr, dict_make_result)
