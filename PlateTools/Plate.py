# -*- coding: utf-8 -*-
"""
Interface for Plates

A Plate represents a collection of Reads (see Read.py) and some additional
information.  This generic interface is intended to be used by classes that
provide more details and capabilities specific to the types of plates used or
other experimental methodologies.

"""

__author__ = "Brian Connelly <bdc@msu.edu>"
__credits__ = "Brian Connelly"

from math import sqrt
import re
import string

from PlateTools.Read import *

class Plate(object):
    """ Information about a plate in an experiment

    Properties:

    num_wells
        The number of wells in the plate
    num_rows
        The number of rows in the plate
    num_columns
        The number of columns in the plate
    info
        Dictionary of information related to the plate
    reads
        An array of Read objects describing reads of this plate

    """

    def __init__(self, num_wells, num_rows=None, num_columns=None):
        """ Initialize a Plate Object

        Arguments:

        num_wells
            The number of wells in the plate
        num_rows
            The number of rows in the plate (if not provided, this will be
            calculated based on the number of wells and assuming a 2:3 matrix)
        num_columns
            The number of columns in the plate (if not provided, this will be
            calculated based on the number of wells and assuming a 2:3 matrix)

        """
        self.num_wells = num_wells

        if num_rows:
            self.num_rows = num_rows
        else:
            self.num_rows = int(sqrt(self.num_wells / 6) * 2)

        if num_columns:
            self.num_columns = num_columns
        else:
            self.num_columns = int(sqrt(self.num_wells / 6) * 3)

        if self.num_rows * self.num_columns != self.num_wells:
            print("ERROR: Invalid number of wells, rows, or columns")

        self.info = {}
        self.reads = []

    def add_read(self, read):
        """ Add a new Read to this Plate's data set """
        self.reads.append(read)

    def num_reads(self):
        """ Return the number of reads for this plate """
        return len(self.reads)

    def get_well_index(self, wells):
        """ For each of the given wells, return its location as a (row, column)
        tuple

        Parameters:

        wells
            A list of wells, specified as e.g. 'A4', 'G10'

        """

        results = []
        for w in wells:
            m = re.match('^([A-ZA-z])([0-9]{1,2})$', w)
            if m:
                row = int(string.uppercase.find(m.group(1).upper()))
                col = int(m.group(2))-1

                if row >= self.num_rows or col >= self.num_columns:
                    # TODO: throw exception
                    print("ERROR: Invalid well '{w}'.  Skipping.".format(w=w))
                else:
                    results.append((row,col))

            else:
                print("ERROR: Invalid well '{w}'.  Skipping.".format(w=w))

        return results

    def get_well_names(self, well_indices):
        """ For each of the given well indices, return it's well name (e.g. 'A4')

        Parameters:

        well_indices:
            A list of well indices, specified as integers

        """

        results = []
        for wi in well_indices:
            well = "{r}{c}".format(r=chr(ord('A') + int(wi/self.num_columns)),
                                   c=wi % self.num_columns)
            results.append(well)

        return results
