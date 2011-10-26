# -*- coding: utf-8 -*-
"""
Interface for Plate Reads

A Read is a collection of data associated with one read of a Plate (see
Plate.py).  Data consists of an array of values associated with each well.
More specific subclasses of Read may include more information (e.g. time or
temperature) specific to a particular machine or experiment.

"""

__author__ = "Brian Connelly <bdc@msu.edu>"
__credits__ = "Brian Connelly"

import csv
import numpy
import sys

from PlateTools.Plate import *

class Read(object):
    """ Information from a plate read

    Properties:

    plate
        A reference to the Plate object being read
    data
        An array of data corresponding to the read values for each well
    info
        A dictionary of additional information.

    """

    def __init__(self, plate):
        """ Initialize a Read object

        Parameters:

        plate
            A reference to the plate read in this read

        """
        self.plate = plate
        self.data = None
        self.info = {}

    def well_values(self, coords):
        """ Retrieve the values at the given wells

        Parameters:

        coords
            A list of (row,column) tuples

        """

        results = []
        for c in coords:
            results.append(self.data[c[0]][c[1]])

        return results

    def csv(self, fp=sys.stdout, delimiter=',', transpose=False):
        """ Print well data in CSV format
        
        Parameters:
        
        fp:
            A file descriptor to print to
        delimiter
            The delimiter character for the outputted CSV (default: ',')
        transpose
            Whether or not to first transpose the data (default: False)

        """
        writer = csv.writer(fp, delimiter=delimiter)
        
        if transpose:
            data = numpy.transpose(self.data)
        else:
            data = self.data

        for row in data:
            writer.writerow(row)

    def pretty(self, fp=sys.stdout, transpose=False):
        """ Print a nice representation of the read

        Parameters:

        fp:
            A file descriptor to print to
        transpose
            Whether or not to first transpose the data (default: False)

        """

        if transpose:
            data = numpy.transpose(self.data)
        else:
            data = self.data
        # TODO: implement
