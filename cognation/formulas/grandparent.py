from __future__ import unicode_literals
from .base import Formula, LineFormatException, AllelesException


class GrandParentFormula(Formula):
    def calculate_relation(self, raw_values):
        if len(raw_values) < 3:
            raise LineFormatException()

        # gc - grandchild, gp - grandparent
        gc_alleles = self.split_sat(raw_values.pop())
        gp_alleles = self.split_sat(raw_values.pop())
        locus = ' '.join(raw_values)  # for loci names contain space

        if len(gp_alleles) != 2 or len(gc_alleles) != 2:
            raise AllelesException()

        gc_set = set(gc_alleles)  # unique grandparent alleles
        gp_set = set(gp_alleles)  # unique grandchild alleles
        intersection = list(gc_set & gp_set)

        len_gc_set = len(gc_set)
        len_gp_set = len(gp_set)
        len_inter = len(intersection)

        freq_dict = self.get_frequencies(locus, gc_alleles)

        if len_gc_set == 1:
            # Homozygous grandchild
            freq = freq_dict[gc_alleles[0]]
            refutation = self._F_(freq)**2
            if len_inter == 0:
                # no common alleles
                confirmation = freq * self._Q_(freq)
            else:
                if len_gp_set == 1:
                    # homozygous grandparent
                    confirmation = self._F_(freq)
                else:
                    # heterozygous grandparent
                    confirmation = self._F_(freq) * self._Q_(freq)
        else:
            # heterozygous grandchild
            freq1 = freq_dict[gc_alleles[0]]
            freq2 = freq_dict[gc_alleles[1]]
            refutation = 2 * self._F_(freq1) * self._F_(freq2) - (2 * freq1 * freq2)**2

            if len_inter == 0:
                # no common alleles
                confirmation = freq1 * self._F_(freq2) + freq2 * self._F_(freq1)
            elif len_inter == 2:
                # both heterozygous, 2 common alleles
                confirmation = self._Q_(freq1) * self._F_(freq2) + self._Q_(freq2) * self._F_(freq1) - freq1 * freq2 * (freq1 + freq2)
            else:
                if len_gp_set == 1:
                    # homozygous grandparent
                    confirmation = self._F_(freq2) + freq2 * (self._F_(freq1) - 2 * freq1 * freq2)
                else:
                    # heterozygous grandparent
                    confirmation = self._Q_(freq1) * self._F_(freq2) + freq2 * self._F_(freq1) - freq1 * (freq2) **2

        lr = confirmation / refutation

        return self.make_result(locus, '/'.join(gc_alleles), '/'.join(gp_alleles), lr)

    # Helpers for frequently used patterns
    @staticmethod
    def _F_(freq):
        return freq * (2 - freq)

    @staticmethod
    def _Q_(freq):
        return 0.5 + 0.5 * freq


"""
class GrandParentFormula(Formula):
    def calculate_relation(self, raw_values):
        if len(raw_values) < 3:
            # skip line with warning
            raise LineFormatException()

        grandchild_alleles = self.split_sat(raw_values.pop())
        grandparent_alleles = self.split_sat(raw_values.pop())
        locus = ' '.join(raw_values)

        if len(grandchild_alleles) != 2 or len(grandparent_alleles) != 2:
            # skip line with warning
            raise AllelesException()
        # GC - grandchild, GP - grandparent
        GC_allele1 = grandchild_alleles[0]
        GC_allele2 = grandchild_alleles[1]
        GP_allele2 = grandparent_alleles[1]
        GP_allele1 = grandparent_alleles[0]

        GC_genotype = GC_allele1+'/'+GC_allele2
        GP_genotype = GP_allele1 + '/' + GP_allele2

        lr = 0
        # Child AA
        if GC_allele2 == GC_allele1:
            # Grand parent AA
            if GP_allele1 == GC_allele2 and GP_allele2 == GC_allele2:
                lr = self.gen_aa_aa(locus, GC_allele2)
            # Grand parent AB
            elif GP_allele1 == GC_allele2 or GP_allele2 == GC_allele2:
                lr = self.gen_aa_ab(locus, GC_allele2)
            # Grand parent BB or BC
            else:
                lr = self.gen_aa_bc(locus, GC_allele2)
        # Child AB
        elif GP_allele1 == GP_allele2:
            # GrandParent AA
            if GC_allele2 == GP_allele1:
                lr = self.gen_ab_aa(locus, GP_allele1, GC_allele1)
            # GrandParent AA
            elif GC_allele1 == GP_allele1:
                lr = self.gen_ab_aa(locus, GP_allele1, GC_allele2)
            # GrandParent CC
            else:
                lr = self.gen_ab_cd(locus, GC_allele2, GC_allele1)
        # GrandParent AB
        elif (GC_allele2 == GP_allele1 and GC_allele1 == GP_allele2) or (GC_allele2 == GP_allele2 and GC_allele1 == GP_allele1):
            lr = self.gen_ab_ab(locus, GC_allele2, GC_allele1)
        else:
            if GC_allele2 == GP_allele1 or GC_allele2 == GP_allele2:
                lr = self.gen_ab_ac(locus, GC_allele2, GC_allele1)
            elif GC_allele1 == GP_allele1 or GC_allele1 == GP_allele2:
                lr = self.gen_ab_ac(locus, GC_allele1, GC_allele2)
            else:
                lr = self.gen_ab_cd(locus, GC_allele2, GC_allele1)

        return self.make_result(locus, GP_genotype, GC_genotype, lr)

    # Child AA Formulas
    def gen_aa_aa(self, locus, a):
        p = self.get_frequencies(locus, {a: 0})
        return 0 if p[a] == 0 else self._2pa_sub_pa2(p, a) / self.prob_not_c_aa(p, a)

    def gen_aa_ab(self, locus, a):
        p = self.get_frequencies(locus, {a: 0})
        return 0 if p[a] == 0 else self._05_add_05pa(p, a) * self._2pa_sub_pa2(p, a) / self.prob_not_c_aa(p, a)

    def gen_aa_bc(self, locus, a):
        p = self.get_frequencies(locus, {a: 0})
        return 0 if p[a] == 0 else p[a] * self._2pa_sub_pa2(p, a) / self.prob_not_c_aa(p, a)

    # Child AB Formulas
    def gen_ab_aa(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            (self._2pa_sub_pa2(p, b) + p[b] * (self._2pa_sub_pa2(p, a) - self._2_pa_pb(p, a, b))) / divider

    def gen_ab_ab(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            (
                    self._05_add_05pa(p, a) * self._2pa_sub_pa2(p, b)
                    + self._05_add_05pa(p, b) * self._2pa_sub_pa2(p, a)
                    - .5 * (p[a] + p[b]) * self._2_pa_pb(p, a, b)
            ) / divider

    def gen_ab_ac(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            (
                    self._05_add_05pa(p, a) * self._2pa_sub_pa2(p, b)
                    + p[b] * self._2pa_sub_pa2(p, a)
                    - .5 * p[b] * self._2_pa_pb(p, a, b)
            ) / divider

    def gen_ab_cd(self, locus, a, b):
        p = self.get_frequencies(locus, {a: 0, b: 0})
        divider = self.prob_not_c_ab(p, a, b)
        return 0 if divider == 0 else \
            (p[a] * self._2pa_sub_pa2(p, b) + p[b] * self._2pa_sub_pa2(p, a)) / divider
"""