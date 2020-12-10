# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.shortcuts import render
from django.http import JsonResponse
from django.template.defaultfilters import register
from .formulas.base import LineFormatException, AllelesException, DelimitingFirst, TooManyDelimitingSymbols, UnknownSymbolInAlleles, DelimitingLast, LociSetDoesNotEqual
import json
# from django.core.serializers.json import DjangoJSONEncoder

from .formulas import formula_builder
from .models import Locus

formulas = [
    {
        'name': 'Отцовство/материнство для двух участников',
        'desc': """
        Исследованию подвергаются биологические образцы проверяемого ребенка и одного
        предполагаемого родителя. Биологический образец второго родителя для исследования
        недоступен""",
        'image': 'ParentFormula.png',
        'value': 1,
        'participants_number': 2,
        'headers': ['Введите гаплотип предполагаемого родителя', 'Введите гаплотип проверяемого лица'],
    },
    {
        'name': 'Двое детей, предполагаемый родитель',
        'desc': """
        Исследованию подвергают биологически образцы проверяемых детей и одного предполагаемого
        родителя. Биологическое родство между проверяемыми детьми достоверно известно. Биологический
        образец второго родителя для исследования недоступен""",
        'image': 'TwoChildrenFormula.png',
        'value': 2,
        'participants_number': 3,
        'headers': ['Введите гаплотип предполагаемого родителя', 'Введите гаплотип первого проверяемого лица',
                    'Введите гаплотип второго проверяемого лица']
    },
    {
        'name': 'Проверяемое лицо, предполагаемые родитель и брат/сестра',
        'desc': """
    Исследованию подвергают биологические образцы проверяемого лица, одного предполагаемого
    родителя и одного предполагаемого брата или сестры. Предполагается, что проверяемое лицо
    и предполагаемый брат (сестра) имеют общих родителей, биологический образец одного из
    которых для исследования недоступен""",
        'image': 'SiblingFormula.png',
        'value': 3,
        'participants_number': 3,
        'headers': ['Введите гаплотип предполагаемого родителя', 'Введите гаплотип предполагаемого брата/сестры',
                    'Введите гаплотип проверяемого лица']
    },
    {
        'name': 'Трое детей и предполагаемый родитель',
        'desc': """
    Исследованию подвергают биологические образцы проверяемых детей и одного предполагаемого
    родителя. Биологическое родство между проверяемыми детьми достоверно известно. Биологический
    образец второго родителя для исследования недоступен""",
        'image': 'ThreeChildrenFormula.png',
        'value': 4,
        'participants_number': 4,
        'headers': ['Введите гаплотип предполагаемого родителя', 'Введите гаплотип первого проверяемого лица',
                    'Введите гаплотип второго проверяемого лица', 'Введите гаплотип третьего проверяемого лица']
    },
    {
        'name': 'Двое детей, известный и предполагаемый родители',
        'desc': """
Исследованию подвергают биологические образцы проверяемых детей и одного предполагаемого
родителя. Биологическое родство между проверяемыми детьми достоверно известно. Биологический
образец второго родителя для исследования недоступен""",
        'image': 'ThreeChildrenFormula.png',
        'value': 5,
        'participants_number': 4,
        'headers': ['Введите гаплотип известного родителя', 'Введите гаплотип первого проверяемого лица',
                    'Введите гаплотип второго проверяемого лица', 'Введите гаплотип предполагаемого родителя'],
    },
    {
        'name': 'Трое детей, известный и предполагаемый родители',
        'desc': """
Исследованию подвергают биологические образцы проверяемых детей, одного известного и одного предполагаемого 
родителя. Биологическое родство между проверяемыми детьми и известным родителем достоверно известно""",
        'image': 'ThreeKnownSupposed.png',
        'value': 6,
        'participants_number': 5,
        'headers': ['Введите гаплотип известного родителя', 'Введите гаплотип первого проверяемого лица',
                    'Введите гаплотип второго проверяемого лица', 'Введите гаплотип третьего проверяемого лица',
                    'Введите гаплотип предполагаемого родителя'],
    },
    {
        'name': 'Ребенок и cупружеская пара',
        'desc': """
Исследованию подвергают биологические образцы проверяемого ребенка и двух предполагаемых родителей. 
Достоверно известно, что предполагаемые родители являются родительской парой, а не отдельными лицами, 
с которыми проверяется родство""",
        'image': 'CoupleFormula.png',
        'value': 7,
        'participants_number': 3,
        'headers': ['Введите гаплотип предполагаемого отца', 'Введите гаплотип предполагаемой матери',
                    'Введите гаплотип проверяемого лица'],
    },
    {
        'name': 'Ребенок, известный и предполагаемый родители',
        'desc': """
Исследованию подвергают биологические образцы проверяемого ребенка, одного известного 
родителя и одного предполагаемого родителя. Биологическое родство между проверяемым ребенком и 
известным родителем достоверно известна""",
        'image': 'OneKnownSupposed.png',
        'value': 8,
        'participants_number': 3,
        'headers': ['Введите гаплотип известного родителя', 'Введите гаплотип проверяемого лица',
                    'Введите гаплотип предполагаемого родителя']
    },
    {
        'name': 'Двое детей и cупружеская пара',
        'desc': """
Исследованию повергают биологические образцы проверяемых детей и двух предполагаемых родителей. 
Достоверно известно, что предполагаемые родители являются родительской парой, а не отдельными 
лицами, с которыми проверяется родство. Также достоверно известно биологическое родство между 
проверяемыми детьми""",
        'image': 'TwoCoupleFormula.png',
        'value': 9,
        'participants_number': 4,
        'headers': ['Введите гаплотип предполагаемого отца', 'Введите гаплотип предполагаемой матери',
                    'Введите гаплотип первого проверяемого лица', 'Введите гаплотип второго проверяемого лица'],
    },
    {
        'name': 'Трое детей и cупружеская пара',
        'desc': """
Исследованию подвергают биологические образцы проверяемых детей и двух предполагаемых родителей. 
Достоверно известно, что предполагаемые родители являются родительской парой, а не отдельными лицами, 
с которыми проверяется родство. Также достоверно известно биологическое родство между проверяемыми детьми""",
        'image': 'ThreeCoupleFormula.png',
        'value': 10,
        'participants_number': 5,
        'headers': ['Введите гаплотип предполагаемого отца', 'Введите гаплотип предполагаемой матери',
                    'Введите гаплотип первого проверяемого лица', 'Введите гаплотип второго проверяемого лица',
                    'Введите гаплотип третьего проверяемого лица'],
    },
    {
        'name': 'Сводные/единокровные братья/сестры',
        'desc': """
Исследованию подвергают биологические образцы лиц, имеющих одного общего родителя. Это могут быть единокровные 
братья/сестры, при условии, что у них общий отец или единоутробные - если у них общая мать""",
        'image': '',
        'value': 11,
        'participants_number': 2,
        'headers': ['Введите гаплотип проверяемого лица', 'Введите гаплотип единокровного/сводного брата/сестры']
    },
    {
        'name': 'Родные братья/сестры',
        'desc': """
Исследованию подвергают биологические образцы проверяемого лица и одного предполагаемого брата или сестры. 
Предполагается, что исследуемые лица имеют общих родителей, биологические образцы которых для исследования 
недоступны""",
        'image': 'BrotherFormula.png',
        'value': 12,
        'participants_number': 2,
        'headers': ['Введите гаплотип проверяемого лица', 'Введите гаплотип предполагаемого брата/сестры'],
    },
    {
        'name': 'Проверяемое лицо и двое братьев/сестер',
        'desc': """
Исследованию подвергают биологические образцы проверяемого лица и двух предполагаемых братьев или сестер. 
Предполагается, что исследуемые лица имеют общих родителей, биологические образцы которых для исследования 
недоступны. Экспертная оценка более достоверна, чем при сравнении гаплотипов исследуемого лица и одного 
брата или сестры""",
        'image': 'TwoBrothersFormula.png',
        'value': 13,
        'participants_number': 3,
        'headers': ['Введите гаплотип проверяемого лица', 'Введите гаплотип первого предполагаемого брата/сестры',
                    'Введите гаплотип второго предполагаемого брата/сестры'],
    },
    {
        'name': 'Двоюродные братья/сестры',
        'desc': """
Исследованию подвергают биологические образцы проверяемого лица и предполагаемого двоюродного брата или сестры. 
Предполагается, что исследуемые лица имеют общих бабушку и дедушку, а один из родителей проверяемого лица является
родным братом/сестрой одному из родителей второго участника. Биологические образцы родственников участников для 
исследования недоступны""",
        'image': '',
        'value': 14,
        'participants_number': 2,
        'headers': ['Введите гаплотип проверяемого лица', 'Введите гаплотип двоюродного брата/сестры'],
    },
    {
        'name': 'Дядя/тетя и племянник/племянница',
        'desc': """
Исследованию подвергают биологические образцы проверяемого лица и дяди/тети. Предполагается, что один из родителей
проверяемого лица является родным братом или сестрой второго учасника. Биологический образец этого родителя для
исследования недоступен""",
        'image': '',
        'value': 15,
        'participants_number': 2,
        'headers': ['Введите гаплотип проверяемого лица', 'Введите гаплотип предполагаемого дяди/тети'],
    },
    {
        'name': 'Проверяемое лицо и предполагаемая бабушка/дедушка',
        'desc': """
Исследованию подвергают биологические образцы проверяемого лица и предполагаемой бабушки(дедушки). Предполагается, 
что один из родителей проверяемого лица является дочерью или сыном предполагаемой бабушки (дедушки). Биологические 
образцы предполагаемых родителей проверяемого лица для исследования недоступны""",
        'image': 'GrandParentFormula.png',
        'value': 17,
        'participants_number': 2,
        'headers': ['Введите гаплотип проверяемого лица', 'Введите гаплотип предполагаемой бабушки/дедушки'],
    },
    {
        'name': 'Проверяемое лицо и предполагаемая бабушка/дедушка',
        'desc': """
Исследованию подвергают биологические образцы проверяемого лица и предполагаемой бабушки(дедушки). Предполагается, 
что один из родителей проверяемого лица является дочерью или сыном предполагаемой бабушки (дедушки). Биологические 
образцы предполагаемых родителей проверяемого лица для исследования недоступны""",
        'image': '',
        'value': 19,
        'participants_number': 3,
        'headers': ['Введите гаплотип предполагаемой бабушки/дедушки', 'Введите гаплотип предполагаемого родителя',
                    'Введите гаплотип проверяемого лица'],
    }
]


def index(request):
    args = {'formulas': formulas, 'formulas_dump': json.dumps(formulas)}
    return render(request, 'cognation/index.html', args)


def calculate(request):
    data_type = request.POST.get('type', 0)
    participants_number = 0
    participants_data = []

    for formula_data in formulas:
        if formula_data['value'] == int(data_type):
            headers = formula_data['headers']
            participants_number = len(headers)

            for i in range(len(headers)):
                each_participant_data = request.POST.get('part'+str(i+1), '')
                participants_data.append(each_participant_data)
            break

    formula = formula_builder(data_type, participants_data)
    try:
        result = formula.calculate()
        cpi = 1
        mutations = 0
        for key in result:
            if "lr" in result[key]:
                if result[key]["lr"] == '-':
                    continue
                elif result[key]["lr"] > 0:
                    cpi *= round(result[key]["lr"], 2)
                else:
                    mutations += 1
        if mutations > 2:
            cpi = 0
            mutations = 0
        ctx = {
            'result': result,
            'cpi': cpi,
            'prob': ((cpi / (1. + cpi)) * 100.),
            'participants': participants_number,
            'mutations': mutations,
            'order': make_order(result),
            'original': {
                'type': data_type,
                'data': participants_data,
            }
        }

        return render(request, formula.get_template(), ctx)
    except (LineFormatException, AllelesException, UnknownSymbolInAlleles,
            TooManyDelimitingSymbols, DelimitingLast, DelimitingFirst, LociSetDoesNotEqual) as exception:
        return render(request, 'cognation/exception.html', {"exception": exception})


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
    return ['AMEL', 'D3S1358', 'vWA', 'D16S539', 'CSF1PO', 'TPOX', 'D8S1179', 'D21S11', 'SE33', 'D18S51', 'Penta E',
            'D2S441', 'D19S433', 'TH01', 'FGA', 'D22S1045', 'D5S818', 'D13S317', 'D7S820', 'D6S1043', 'D10S1248',
            'D1S1656', 'D12S391', 'D2S1338', 'Penta D', 'Yindel', 'DYS391', 'SRY']


@register.filter
def get_value(value, key):
    return value.get(key)
