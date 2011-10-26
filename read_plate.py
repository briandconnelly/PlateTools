#!/usr/bin/env arch -i386 python
# TODO: change this line back to what it should be
# -*- coding: utf-8 -*-

__author__ = "Brian Connelly <bdc@msu.edu>"
__version__ = 0.1

from PlateTools.formats import SoftMaxPro

import argparse
import numpy
import sys

def main():
    parser = argparse.ArgumentParser(description='Interact with data from a SoftMax(R) Pro experiment file',
                                     epilog='SoftMax is a registered trademark of Molecular Devices, LLC.')
    parser.add_argument('infile', type=argparse.FileType('rb'), help='read data from given file')
    parser.add_argument('-c', '--cuvette', action='append', help='use specific cuvette(s)', default=None)
    parser.add_argument('-g', '--group', action='append', help='use specific group(s)', default=None)
    parser.add_argument('-i', '--info', action='store_true', default=False, help='display information about the experiment or selected plate/cuvette')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout, help='write data to given file')
    parser.add_argument('-p', '--plate', action='append', help='use specific plate(s)', default=None)
    parser.add_argument('-P', '--pretty', action='store_true', default=False, help='display plate read data as a formatted table')
    parser.add_argument('-T', '--transpose', action='store_true', default=False, help='transpose resulting plate read data')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))
    args = parser.parse_args()
    
    # Read in the experiment data
    try:
        experiment = SoftMaxPro.Experiment(args.infile)
    except Execption as err:
        print("ERROR: Could not read input file:", err)
        sys.exit(1)
    else:
        plates_specified = False
        cuvettes_specified = False
        groups_specified = False

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

        target_cuvettes = []
        if args.cuvette and len(args.cuvette) > 0:
            cuvettes_specified = True

            for c in args.cuvette:
                try:
                    target_cuvettes.append(experiment.plates[c])
                except KeyError as err:
                    if c.lower() == "all":
                        target_cuvettes = list(experiment.cuvettes.values())
                    else:
                        print("Error: Cuvette '{0}' does not exist in experiment".format(c))
            target_cuvettes = set(target_cuvettes)
        else:
            target_cuvettes = list(experiment.cuvettes.values())

        target_groups = []
        if args.group and len(args.group) > 0:
            groups_specified = True
            for g in args.group:
                try:
                    target_groups.append(experiment.groups[g])
                except KeyError as err:
                    if g.lower() == "all":
                        target_groups = list(experiment.groups.values())
                    else:
                        print("Error: Group '{0}' does not exist in experiment".format(g))
            target_groups = set(target_groups)
        else:
            target_groups = list(experiment.groups.values())


        if not plates_specified and not cuvettes_specified and not groups_specified:
            num_notes = len(experiment.notes)
            print('{0}: {1}'.format('Number of Notes', num_notes))
            for n in experiment.notes:
                print('\t* {0}: {1}'.format('Note', n))

            num_plates = len(experiment.plates)
            print('{0}: {1}'.format('Number of Plates', num_plates))
            for k,v in experiment.plates.iteritems():
                print("\t* {0}".format(v))

            num_cuvettes = len(experiment.cuvettes)
            print('{0}: {1}'.format('Number of Cuvettes', num_cuvettes))
            for k,v in experiment.cuvettes.iteritems():
                print("\t* {0}".format(v))

            num_groups = len(experiment.groups)
            print('{0}: {1}'.format('Number of Groups', num_groups))
            for k,v in experiment.groups.iteritems():
                print("\t* {0}".format(v))

        else:
            if plates_specified:
                for p in target_plates:
                    if args.info:
                        p.print_information()
                    else:
                        for r in p.reads:
                            print("{0}".format(p))
                            if args.pretty:
                                r.pretty()
                            else:
                                r.csv(fp=args.outfile, transpose=args.transpose)
                            print("")
                        print("")

            if cuvettes_specified:
                for c in target_cuvettes:
                    # TODO: implement
                    if args.info:
                        pass
                    else:
                        pass

            if groups_specified:
                for g in target_groups:
                    if args.info:
                        g.print_information()
                    else:
                        g.print_data()

if __name__ == "__main__":
    main()
