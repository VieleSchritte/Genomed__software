from __future__ import unicode_literals
from .base import Formula
from .base import AllelesException
from .base import Calculations


class BrotherFormula(Formula):
    def calculate_relation(self, raw_values):
        (brother_alleles, insp_alleles, locus, brother_set, insp_set, intersection) = self.getting_alleles_locus(raw_values, 2)

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr
        if self.is_gender_specific(locus):
            return self.make_result2(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), '-')

        if len(brother_alleles) != 2 or len(insp_alleles) != 2:
            raise AllelesException()

        c = Calculations()
        # In cases aa aa or ab ab lr = 1
        if len(intersection) == 2 or len(intersection) == len(brother_set) == len(insp_set) == 1:
            confirmation = 1
            freq_dict = self.get_frequencies(locus, list(insp_set))

            # Homozygous inspected person
            if len(insp_set) == 1:
                freq = freq_dict[insp_alleles[0]]
                refutation = c.homo_refutation(freq)

            #  Heterozygous inspected person
            else:
                freq1, freq2 = freq_dict[insp_alleles[0]], freq_dict[insp_alleles[1]]
                refutation = c.hetero_refutation(freq1, freq2)

            lr = confirmation / refutation

            return self.make_result2(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), lr)

        freq_dict = self.get_frequencies(locus, insp_alleles + brother_alleles)
        conf = Confirmations()
        calc = Calculations()

        #  Homozygous inspected person
        if len(insp_set) == 1:
            freq = freq_dict[insp_alleles[0]]
            confirmation = conf.homo_insp_conf(freq_dict, intersection, brother_set, brother_alleles, insp_alleles)

            refutation = calc.homo_refutation(freq)

        #  Heterozygous inspected person
        else:
            freq1, freq2 = freq_dict[insp_alleles[0]], freq_dict[insp_alleles[1]]
            refutation = calc.hetero_refutation(freq1, freq2)
            confirmation = conf.hetero_insp_conf(freq_dict, intersection, brother_set, brother_alleles, insp_alleles)

        lr = confirmation / refutation
        return self.make_result2(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), lr)


class Confirmations:
    def homo_insp_conf(self, freq_dict, intersection, brother_set, brother_alleles, insp_alleles):
        calc = Calculations()
        freq1 = freq_dict[insp_alleles[0]]

        #  No common alleles
        if len(intersection) == 0:
            #  Homozygous brother
            if len(brother_set) == 1:
                freq2 = freq_dict[brother_alleles[0]]
                return (calc.M(freq2, freq1)) ** 2

            #  Heterozygous brother
            freq2, freq3 = freq_dict[brother_alleles[0]], freq_dict[brother_alleles[1]]
            return 4 * freq1 ** 2 * freq2 * freq3 / (calc.F(freq2) * calc.F(freq3) - 2 * (freq2 * freq3) ** 2)

        #  One common allele
        allele2, allele1 = self.get_alleles_order(brother_alleles, list(intersection))
        freq1, freq2 = freq_dict[allele1], freq_dict[allele2]
        return 2 * freq1 * freq2 * (2 * calc.F(freq1) - freq1 * freq2) / (calc.F(freq1) * calc.F(freq2) - 2 * (freq1 * freq2) ** 2)

    def hetero_insp_conf(self, freq_dict, intersection, brother_set, brother_alleles, insp_alleles):
        calc = Calculations()

        #  1 common allele
        if len(intersection) == 1:

            #  Homozygous brother
            if len(brother_set) == 1:
                (allele2, allele1) = self.get_alleles_order(insp_alleles, brother_alleles)
                freq1, freq2 = freq_dict[allele1], freq_dict[allele2]
                return (4 * freq1 * freq2 * calc.F(freq1) - (2 * freq1 * freq2) ** 2) / (calc.F(freq1)) ** 2

            #  Heterozygous brother a (1) = inter, b (2) = unique insp, c (3) = unique bro
            (allele2, allele1) = self.get_alleles_order(insp_alleles, list(intersection))
            (allele3, allele1) = self.get_alleles_order(brother_alleles, list(intersection))
            freq1, freq2, freq3 = freq_dict[allele1], freq_dict[allele2], freq_dict[allele3]
            return 2 * freq2 * freq3 * (calc.F(freq1) - 2 * freq1 ** 2) / (calc.F(freq1) * calc.F(freq3) - 2 * (freq1 * freq3) ** 2)

        #  No common alleles
        #  Homozygous brother
        if len(brother_set) == 1:
            freq1, freq2, freq3 = freq_dict[insp_alleles[0]], freq_dict[insp_alleles[1]], freq_dict[brother_alleles[0]]
            return 2 * calc.M(freq3, freq1) * calc.M(freq3, freq2)

        #  Heterozygous brother
        freq1, freq2 = freq_dict[insp_alleles[0]], freq_dict[insp_alleles[1]]
        freq3, freq4 = freq_dict[brother_alleles[0]], freq_dict[brother_alleles[1]]
        return 2 * freq1 * freq2 * freq3 * freq4 / (calc.F(freq3) * calc.F(freq4) - 2 * (freq3 * freq4) ** 2)

    @staticmethod
    def get_alleles_order(run_list, check_list):
        allele1 = ''
        allele2 = ''

        for allele in run_list:
            if allele not in check_list:
                allele1 = allele
            else:
                allele2 = allele
        #  unique allele returns first
        return allele1, allele2
