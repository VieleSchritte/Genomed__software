# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import render
from django.http import JsonResponse
from django.template.defaultfilters import register

from .formulas import formula_builder
from .models import Locus

def index(request):
    return render(request, 'cognation/index.html')


def calculate(request):
    data_type = request.POST.get('type', 0)
    data = request.POST.get('data', '')

    formula = formula_builder(data_type, data)
    result = formula.calculate()

    cpi = 1
    participants = 0
    mutations = 0

    for key in result:
        if participants == 0:
            participants = 3 if "ef" in result[key] else 2

        if "lr" in result[key]:
            if result[key]["lr"] == '-':
                continue
            elif result[key]["lr"] > 0:
                cpi = cpi * result[key]["lr"]
            else:
                mutations = mutations + 1

    if mutations > 2:
        cpi = 0
        mutations = 0

    ctx = {
        'result': result,
        'cpi': cpi,
        'prob': ((cpi / (1. + cpi)) * 100.),
        'participants': participants,
        'mutations': mutations,
        'order': make_order(result),
        'original': {
            'type': data_type,
            'data': data,
        }
    }

    return render(request, formula.get_template(), ctx)


def save_allele(request):
    locus = request.POST.get('locus')
    sat = float(request.POST.get('sat'))
    freq = float(request.POST.get('freq'))

    try:
        locus_object = Locus.objects.get(locus=locus, sat=sat)
        locus_object.freq = freq
        locus_object.save()
    except Locus.DoesNotExist:
        locus_object = Locus(locus=locus, sat=sat, freq=freq)
        locus_object.save()

    return JsonResponse({'result': 'ok'})


def make_order(results):
    if len(results) <= 20:
        return make_order_cordis()
    else:
        return make_order_verifiler()


def make_order_cordis():
    return ['AMEL', 'D3S1358', 'TH01',  'D12S391', 'D1S1656',  'D10S1248', 'D22S1045',  'D2S441',  'D7S820',  'D13S317',
            'FGA',  'TPOX', 'D18S51', 'D16S539', 'D8S1179',  'CSF1PO',  'D5S818',  'vWA',  'D21S11', 'SE33']


def make_order_verifiler():
    return ['AMEL', 'Yindel', 'D3S1358', 'vWA', 'D16S539', 'CSF1PO', 'TPOX', 'D8S1179', 'D21S11', 'D18S51', 'Penta E', 'D2S441',
            'D19S433', 'TH01', 'FGA', 'D22S1045', 'D5S818', 'D13S317', 'D7S820', 'D6S1043', 'D10S1248', 'D1S1656',
            'D12S391', 'D2S1338', 'Penta D']

def make_order_yfiler():
    return ['DYS576', 'DYS389I', 'DYS635', 'DYS389II', 'DYS627', 'DYS460', 'DYS458', 'DYS19', 'YGATAH4',
            'DYS448', 'DYS391', 'DYS456', 'DYS390', 'DYS438', 'DYS392', 'DYS518', 'DYS570', 'DYS437',
            'DYS385', 'DYS449', 'DYS533', 'DYS393', 'DYS439', 'DYS481', 'DYF387S1']

def make_order_cordis_exp():
    return ['AMEL', 'D3S1358', 'TH01', 'D12S391', 'D5S818', 'TPOX', 'D2S441', 'D7S820', 'D13S317', 'FGA',
            'D22S1045', 'D18S51', 'D16S539', 'D8S1179', 'CSF1PO', 'D6S1043', 'vWA', 'D21S11', 'SE33',
            'D10S1248', 'D1S1656', 'D19S433', 'D2S1338', 'SRY', 'DYS391', 'Yindel']


@register.filter
def get_value(value, key):
    return value.get(key)
