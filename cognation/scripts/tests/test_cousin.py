# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from django.test import TestCase
from ...formulas.cousin import CousinFormula

class TestCousinFormula(TestCase):

    @classmethod
    def setUpTestData(cls):
        print("setUpTestData: Run once to set up non-modified data for all class methods.")
        pass

    def setUp(self):
        print("setUp: Run once for every test method to setup clean data.")
        pass