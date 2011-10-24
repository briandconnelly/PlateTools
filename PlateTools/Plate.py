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
