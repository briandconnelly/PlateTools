#!/usr/bin/env python
# -*- coding: utf-8 -*-

# TODO: format time to be HH:MM:SS.  What if there are days?? anything in manual?
# TODO: add support for groups?

__author__ = "Brian Connelly <bdc@msu.edu>"
__version__ = 0.1

from PlateTools.formats import SoftMaxPro

import argparse
import csv
import numpy as np
import string
import sys

def main():
    parser = argparse.ArgumentParser(description='Convert plate data from a SoftMax® Pro experiment file to CSV',
                                     epilog='SoftMax is a registered trademark of Molecular Devices, LLC.')

    parser.add_argument('-c', '--comments', dest='comment_char', default='#',
                        metavar='C',
                        help='character to be used for comments (default #)')
    parser.add_argument('-d', '--delimiter', dest='delimiter', default=',',
                        metavar='D', help='delimiter (default ,)')
    parser.add_argument('-H', '--noheader', dest='header', action='store_false',
                        default=True, help='do not write header information')
    parser.add_argument('infile', type=argparse.FileType('rb'),
                        help='read data from given file')
    parser.add_argument('-i', '--info', action='store_true', default=False,
                        help='display information about the experiment or selected plate')
    parser.add_argument('-N', '--nonotes', dest='notes', action='store_false',
                        default=True, help='do not write notes as comments')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout, help='write data to given file')
    parser.add_argument('-p', '--plate', action='append', default=None,
                        help='use specific plate(s)')
    parser.add_argument('-T', '--transpose', action='store_true', default=False,
                        help='transpose resulting plate read data')
    parser.add_argument('-v', '--version', action='version',
                        version='%(prog)s {0}'.format(__version__))
    args = parser.parse_args()

    # Read in the experiment data
    try:
        experiment = SoftMaxPro.Experiment(args.infile)
    except:
        print("ERROR: Could not read input file")
        sys.exit(1)
    else:
        target_plates = []
        if args.plate and len(args.plate) > 0:
            plates_specified = True

            for p in args.plate:
                try:
                    target_plates.append(experiment.plates[p])
                except KeyError as err:
                    if p.lower() == "all":
                        target_plates = list(experiment.plates.values())
                    else:
                        print("Error: Plate '{0}' does not exist in experiment".format(p))
            target_plates = set(target_plates)
        else:
            target_plates = list(experiment.plates.values())

        # TODO: if info specified, show it and quit
        #num_plates = len(experiment.plates)
        #print('{0}: {1}'.format('Number of Plates', num_plates))
        #for k,v in experiment.plates.iteritems():
        #    print("\t* {0}".format(v))

        if args.info:
            for p in target_plates:
                p.print_information()
            sys.exit(0)


        writer = csv.writer(args.outfile, delimiter=args.delimiter)

        plate_num = 0
        for p in target_plates:
            plate_name = string.replace(p.info['name'], '#', '_')

            if plate_num == 0:
                if p.info['data_mode'] == 'Absorbance':
                    plate_type = '{m}'.format(m='Absorbance')
                elif p.info['data_mode'] == 'Luminescence':
                    plate_type = '{m}'.format(m='Luminescence')
                else:
                    plate_type = 'Value'

                if args.header:
                    writer.writerow(['Plate','Time','Temperature','Row','Column',plate_type])

                if args.notes:
                    for n in experiment.notes:
                        if len(n) > 0:
                            writer.writerow(['{c} {n}'.format(c=args.comment_char, n=n)])

            for r in p.reads:
                try:
                    i = np.nditer(r.data)
                except ValueError:
                    # Sometimes bad reads happen.  Skip these.
                    continue
                else:
                    if r.info['time'] == None:
                        read_time = "00:00:00"
                    else:
                        # TODO: format this for HH:MM:SS (force hours)
                        read_time = r.info['time']


                    for row in range(r.data.shape[0]):
                        for col in range(r.data.shape[1]):
                            val = float(r.data[row,col])

                            adj_row  = p.info['first_row_read'] + row - 1
                            adj_col = p.info['first_column_read'] + col - 1

                            writer.writerow([plate_name, read_time, r.info['temperature'], adj_row, adj_col, val])

            plate_num += 1

if __name__ == "__main__":
    main()
