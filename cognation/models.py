# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models


class Locus(models.Model):
    locus = models.CharField(max_length=32)
    sat = models.FloatField()
    freq = models.FloatField()

    class Meta:
        indexes = [
            models.Index(fields=['locus', 'sat'])
        ]

    def __str__(self):
        return self.locus + '  |  ' + str(self.sat) + ' | ' + str(self.freq)
