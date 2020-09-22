# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import unittest
from django.test import TestCase
from .models import Locus
from .formulas import parent

"""
class TestLocus(unittest.TestCase):
    def checkingLocus(self):
        self.assertEqual(Locus.locus.isalpha(), 'locus object in class Locus is a string')
        self.assertEqual(self.is_number(Locus.sat), 'sat object in class Locus is float type')
        self.assertEqual(self.is_number(Locus.freq), 'freq object in class Locus is float type')

    # checking float number
    def is_number(self, sat_freq):
        try:
            float(sat_freq)
            int(sat_freq)
            return True
        except ValueError:
            return False
"""

class