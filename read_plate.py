#!/usr/bin/env arch -i386 python
# TODO: change this line back to what it should be
# -*- coding: utf-8 -*-

# TODO: use action=append for the lists (plate,group,cuvette)

__author__ = "Brian Connelly <bdc@msu.edu>"
__version__ = 0.1

try:
    import numpy
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    pass

import PlateTools
from PlateTools.formats import SoftMaxPro

import argparse
import sys


def main():
    parser = argparse.ArgumentParser(description='Interact with data from a SoftMax(R) Pro experiment file',
                                     epilog='SoftMax is a registered trademark of Molecular Devices, LLC.')
    parser.add_argument('infile', type=argparse.FileType('rb'), help='read data from given file')
    parser.add_argument('-c', '--cuvette', nargs='+', help='use specific cuvette(s)', default=None)
    parser.add_argument('-g', '--group', nargs='+', help='use specific group(s)', default=None)
    parser.add_argument('-i', '--info', action='store_true', default=False, help='display information about the experiment or selected plate/cuvette')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout, help='write data to given file')
    parser.add_argument('-p', '--plate', nargs='+', help='use specific plate(s)', default=None)
    #parser.add_argument('--cuvette_info', metavar='C', help='display information about the given cuvette')
    if HAS_NUMPY:
        parser.add_argument('-T', '--transpose', action='store_true', default=False, help='transpose resulting data')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
    args = parser.parse_args()
    
    # Read in the experiment data
    # TODO: wrap with a try/except.  On except, print an error and quit
    experiment = SoftMaxPro.Experiment(args.infile)

    if args.info:
        if args.plate and len(args.plate) > 0:
            for p in args.plate:
                try:
                    p = experiment.plates[p]
                    p.print_information()
                except KeyError as err:
                    print("Error: Plate '{0}' does not exist in experiment".format(p))
        elif args.cuvette and len(args.cuvette) > 0:
            for c in args.cuvette:
                # TODO: implment
                pass
        elif args.group and len(args.group) > 0:
            for g in args.group:
                # TODO: implment
                pass
        else:
            for n in experiment.notes:
                print('{0}: {1}'.format('Note', n))

            num_plates = len(experiment.plates)
            print('{0}: {1}'.format('Number of Plates', num_plates))
            for k,v in experiment.plates.iteritems():
                print("\t* {0}".format(v))

            num_cuvettes = len(experiment.cuvettes)
            print('{0}: {1}'.format('Number of Cuvettes', num_cuvettes))
            for k,v in experiment.cuvettes.iteritems():
                print("\t* {0}".format(v))

            # TODO: display groups info
    else:
        target_groups = []
        target_plates = []

        if args.plate and len(args.plate) > 0:
            for p in args.plate:
                try:
                    target_plates.append(experiment.plates[p])
                except KeyError as err:
                    print("Error: Plate '{0}' does not exist in experiment".format(p))
        else:
            for k, v in experiment.plates.iteritems():
                target_plates.append(v)

        for p in target_plates:
            if args.transpose:
                pass
            # TODO - print the data.  default to CSV.

        # TODO: do the same for cuvettes

if __name__ == "__main__":
    main()