import openpyxl
import sys
import os
import django
sys.path.append('')
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "paternity_app.settings")
django.setup()

from cognation.models import Locus

if len(sys.argv) != 2:
    print("Wrong arguments count")
    exit(0)

file = sys.argv[1]
if not os.path.isfile(file):
    print("File not exists")
    exit(0)

locus_name = file[file.rfind('/')+1:-5]

work_book = openpyxl.load_workbook(file)
work_sheet = work_book.worksheets[0]
print(work_sheet)
if work_sheet.title[-3:] != "rus":
    print("First sheet not rus")
    exit(0)

for row in work_sheet.iter_rows():
    try:
        sat = row[0].value
        freq = row[1].value
        if sat is None or freq is None:
            continue

        locus = Locus(locus=locus_name, sat=float(sat), freq=float(freq))
        locus.save()
    except ValueError:
        print("Skip " + str(row[0].value) + ", " + str(row[1].value))