from __future__ import unicode_literals
from .base import Formula, Calculations
from .parent import ParentFormula


class TwoChildrenFormula(Formula):
    def calculate_relation(self, raw_values):
        (locus, part_alleles, part_sets, intersections, dict_make_result) = self.getting_alleles_locus(raw_values, 3)
        child1_alleles, child2_alleles, parent_alleles = part_alleles
        child1_set, child2_set, parent_set = part_sets
        ch1p_intersection, ch2p_intersection, ch1ch2_intersection = intersections

        print(locus)
        print('parent_alleles: ', parent_alleles)
        print('child1_alleles: ', child1_alleles)
        print('child2_alleles: ', child2_alleles)

        if self.is_gender_specific(locus):
            print('gender specific')
            print()
            return self.make_result(locus, '-', dict_make_result)

        if child1_set == child2_set:
            print('calls ParentFormula')
            print()
            raw_values = [locus, '/'.join(parent_alleles), '/'.join(child2_alleles)]
            result = ParentFormula(Formula).calculate_relation(raw_values)
            result['part3'] = '/'.join(child1_alleles)
            return result

        c = Calculations()
        freq_dict = self.get_frequencies(locus, child1_alleles + child2_alleles + parent_alleles)
        lr = 0

        if len(ch1p_intersection) >= 1 and len(ch2p_intersection) >= 1:
            # Homozygous 1st child
            if len(child1_set) == 1:
                print('homozygous child')
                freq1, freq2, freq3 = freq_dict[child1_alleles[0]], freq_dict[child2_alleles[0]], freq_dict[child2_alleles[1]]

                # case aa an an
                if len(ch1ch2_intersection) != 0:
                    print('aa an an')
                    print()
                    lr = c.F(freq1)
                    return self.make_result(locus, lr, dict_make_result)

                else:
                    # case aa bb ab
                    if len(child2_set) == 1:
                        print('aa bb ab')
                        print()
                        lr = 2 * freq1 * freq2
                        return self.make_result(locus, lr, dict_make_result)

                    # case aa bc ab/ac
                    else:
                        print('aa bc ab/ac')
                        print()
                        lr = 2 * freq1 * (freq2 + freq3)
                        return self.make_result(locus, lr, dict_make_result)

            # Heterozygous 1st child
            else:
                print('heterozygous child')
                # case ab cc ac/bc
                if len(child2_set) == 1:
                    print('ab cc ac/bc')
                    print()
                    freq1, freq2, freq3 = freq_dict[child2_alleles[0]], freq_dict[child1_alleles[0]], freq_dict[child1_alleles[1]]
                    lr = 2 * freq1 * (freq2 + freq3)
                    return self.make_result(locus, lr, dict_make_result)

                # case ab ac an/bc
                if len(child2_set) == 2 and len(ch1ch2_intersection) == 1:
                    print('ab ac an/bc')
                    print()
                    freq2, freq3 = freq_dict[child1_alleles[0]], freq_dict[child1_alleles[1]]
                    lr = c.F(freq_dict[list(ch1ch2_intersection)[0]]) + 2 * freq2 * freq3
                    return self.make_result(locus, lr, dict_make_result)

                # case ab cd ac/ad/bc/bd
                if len(child2_set) == 2 and len(ch1ch2_intersection) == 0:
                    print('ab cd ac/ad/bc/bd')
                    print()
                    freq1, freq2 = freq_dict[child1_alleles[0]], freq_dict[child1_alleles[1]]
                    freq3, freq4 = freq_dict[child2_alleles[0]], freq_dict[child2_alleles[1]]
                    lr = 2 * (freq1 + freq2) * (freq3 + freq4)
                    return self.make_result(locus, lr, dict_make_result)
        print('no intersections, return lr = 0')
        print()
        return self.make_result(locus, lr, dict_make_result)
