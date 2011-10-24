#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
from distutils.core import setup

try:
    import numpy
except ImportError as err:
    print("Error: PlateTools requires numpy.  Please install it first.")
    sys.exit(1)

try:
    import scipy
except ImportError as err:
    print("Error: PlateTools requires scipy.  Please install it first.")
    sys.exit(1)

try:
    import matplotlib
except ImportError as err:
    print("Error: PlateTools requires matplotlib.  Please install it first.")
    sys.exit(1)

import PlateTools as P

if os.path.exists('MANIFEST'): os.remove('MANIFEST')

if sys.argv[-1] == 'setup.py':
    print("To install, run 'python setup.py install'")
    print()

if sys.version_info[:2] < (2, 7):
    print("PlateTools requires Python version 2.7 or later (%d.%d detected)." %
    sys.version_info[:2])
    sys.exit(-1)

if __name__ == "__main__":

    with open('README.rst') as file:
        ldesc = file.read()

    setup(
        name = "PlateTools",
        version = P.__version__,
        packages = ['PlateTools','PlateTools.formats'],
        scripts = ['scripts/read_experiment.py'],
        license = P.__license__,
        author = "Brian Connelly",
        author_email = "bdc@msu.edu",
        maintainer = "Brian Connelly",
        maintainer_email = "bdc@msu.edu",
        url = "https://github.com/briandconnelly/PlateTools",
        download_url = P.__download_url__,
        keywords = ["microtitre", "microtiter", "microplate", "parser"],
        classifiers = [
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: Apache Software License",
            "Topic :: Scientific/Engineering",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Topic :: Scientific/Engineering :: Medical Science Apps."
        ],
        description = "Tools for reading, manipulating, and analyzing microplate data",
        long_description = ldesc
)

