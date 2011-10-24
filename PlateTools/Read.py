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
        self.data_np = numpy.empty((self.plate.num_rows, self.plate.num_columns))
        print self.data_np
        self.data = []
        self.info = {}

    def csv(fp=sys.stdout):
        """ Print well data in CSV format """
        pass

    def pretty_print(fp=sys.stdout):
        """ Print a nice representation of the read """
        pass
