#!/usr/bin/env python3

import h5py
import csv
import argparse

parser = argparse.ArgumentParser('dpcheck', description='Data pipeline table checker')

parser.add_argument('dp_path', help='data_pipeline hdf5 file')
parser.add_argument('component', help='component name of table in hdf5 file')
parser.add_argument('csv_file', help='csv file to compare against')

args = parser.parse_args()

# CSV file

csv.register_dialect('eera', skipinitialspace=True)

with open(args.csv_file, newline='') as csvfile:
    csvdata = csv.reader(csvfile, dialect='eera')
    csv_headings = csvdata.__next__()
    csv_data = []

    for row in csvdata:
        csv_data.append([float(i) for i in row])

# HDF5 file

table = h5py.File(args.dp_path, 'r')

if (args.component != ''):
    for name in args.component.split('/'):
        table = table[name]

table = table['table']
table_headings = [name for name in table.dtype.names]

# Compare data
# This actually seems really slow, I assume there are better ways...

header_error = False

if (csv_headings != table_headings):
    print('Headings do not match:')
    print(f'    table: {table_headings}')
    print(f'    csv  : {csv_headings}')
    header_error = True

if (len(csv_data) != len(table)):
    print(f'Number of rows mismatch: table {len(table)}, CSV {len(csv_data)}')
    header_error = True

rows_to_check = min(len(csv_data), table.len())
error = False

x = 0

for csv_row in csv_data:
    if (len(csv_row) != len(table[x])):
        print(f'Row {x+1} size mismatch: table {len(table[x])}, csv {len(csv_row)}')
        error = True
    else:
        y = 0

        for csv_el in csv_row:
            if (csv_el != table[x][y]):
                print(f'Row {x+1} el {y+1} mismatch: table {table[x][y]}, csv {csv_el}')
                error = True

            y = y + 1

    x = x + 1

    if (error):
        print('...')
        break

# for x in range(rows_to_check):
#     if (len(csv_data[x]) != len(table[x])):
#         print(f'Row {x+1} size mismatch: table {len(table[x])}, csv {len(csv_data[x])}')
#         error = True

#     cols_to_check = min(len(csv_data[x]), len(table[x]))

#     for y in range(cols_to_check):
#         if (csv_data[x][y] != table[x][y]):
#             print(f'Row {x+1} el {y+1} mismatch: table {table[x][y]}, csv {csv_data[x][y]}')
#             error = True

#     if (error):
#         break

if (error | header_error):
    print(f'DATA UNMATCHED')
    exit(1)
else:
    exit(0)
