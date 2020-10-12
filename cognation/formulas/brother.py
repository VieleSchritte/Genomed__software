from __future__ import unicode_literals
from .base import Formula
from .base import AllelesException
from .base import Calculations


class BrotherFormula(Formula):
    def calculate_relation(self, raw_values):
        (insp_alleles, brother_alleles, locus, insp_set, brother_set, intersection) = self.getting_alleles_locus(raw_values, 2)

        # Function in base.py for checking out if the locus is gender-specific; if yes return lr
        if self.is_gender_specific(locus):
            return self.make_result(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), '-')

        if len(brother_alleles) != 2 or len(insp_alleles) != 2:
            raise AllelesException()

        calc = Calculations()
        lr = 1

        #  Homozygous inspected person
        if len(insp_set) == 1:
            freq_dict = self.get_frequencies(locus, list(insp_set))
            freq1 = freq_dict[insp_alleles[0]]
            refutation = calc.F(freq1)

            #  One common allele
            if len(intersection) == 1:

                #  Homozygous brother
                if len(brother_set) == 1:
                    return self.make_result(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), lr)

                #  Heterozygous brother
                freq_dict = self.get_frequencies(locus, list(brother_set))
                freq2 = 0
                for i in range(len(brother_alleles)):
                    if brother_alleles[i] not in insp_alleles:
                        brother_allele = brother_alleles[i]
                        freq2 = freq_dict[brother_allele]

                confirmation = 2 * freq1 * freq2 * (2 * calc.F(freq1) - freq1 * freq2) / (calc.F(freq1) * calc.F(freq2) - 2 * (freq1 * freq2) ** 2)
                lr = confirmation / refutation
                return self.make_result(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), lr)

            #  No common alleles
            freq_dict = self.get_frequencies(locus, list(brother_set))

            #  homozygous brother
            if len(brother_set) == 1:
                freq2 = freq_dict[0]
                confirmation = (calc.M(freq2, freq1)) ** 2
                lr = confirmation / refutation
                return self.make_result(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), lr)

            #  heterozygous brother
            freq2, freq3 = freq_dict[brother_alleles[0]], freq_dict[brother_alleles[1]]
            confirmation = 4 * freq1 ** 2 * freq2 * freq3 / (calc.F(freq2) * calc.F(freq3) - 2 * (freq2 * freq3) ** 2)
            lr = confirmation / refutation
            return self.make_result(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), lr)

        # Heterozygous inspected person
        if len(insp_set) == 2:

            #  2 common alleles
            if len(intersection) == 2:
                return self.make_result(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), lr)

            freq_dict = self.get_frequencies(locus, insp_alleles)
            freq1, freq2 = freq_dict[insp_alleles[0]], freq_dict[insp_alleles[1]]
            refutation = 2 * calc.F(freq1) * calc.F(freq2) - (2 * freq1 * freq2) ** 2

            #  1 common allele
            if len(intersection) == 1:

                if len(brother_set) == 1:
                    #  case ab aa: freq2 is the frequency of allele b - unique inspected person's allele
                    for allele in insp_alleles:
                        if allele == brother_alleles[0]:
                            freq1 = freq_dict[allele]
                        else:
                            freq2 = freq_dict[allele]

                    confirmation = (4 * freq1 * freq2 * calc.F(freq1) - (2 * freq1 * freq2) ** 2) / (calc.F(freq1)) ** 2
                    lr = confirmation / refutation
                    return self.make_result(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), lr)

                #  Heterozygous brother
                freq_dict = self.get_frequencies(locus, insp_alleles + brother_alleles)
                freq3 = 0

                for allele in brother_alleles:
                    if allele == list(intersection)[0]:
                        freq1 = freq_dict[allele]
                    else:
                        freq3 = freq_dict[allele]

                for allele in insp_alleles:
                    if allele not in brother_alleles:
                        freq2 = freq_dict[allele]

                confirmation = 2 * freq2 * freq3 * (calc.F(freq1) - 2 * freq1 ** 2) / (calc.F(freq1) * calc.F(freq3) - 2 * (freq1 * freq3) ** 2)
                lr = confirmation / refutation
                return self.make_result(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), lr)

            #  No common alleles
            freq_dict = self.get_frequencies(locus, insp_alleles + brother_alleles)

            #  homozygous brother
            if len(brother_set) == 1:
                freq3 = freq_dict[brother_alleles[0]]
                confirmation = 2 * calc.M(freq3, freq1) * calc.M(freq3, freq2)
                lr = confirmation / refutation
                return self.make_result(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), lr)

            #  heterozygous brother
            freq3, freq4 = freq_dict[brother_alleles[0]], freq_dict[brother_alleles[1]]
            confirmation = 2 * freq1 * freq2 * freq3 * freq4 / (calc.F(freq3) * calc.F(freq4) - 2 * (freq3 * freq4) ** 2)
            lr = confirmation / refutation
            return self.make_result(locus, '/'.join(brother_alleles), '/'.join(insp_alleles), lr)
