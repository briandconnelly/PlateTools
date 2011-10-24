# -*- coding: utf-8 -*-
VERSION = (0, 0, 1)
__version__ = ".".join(map(str, VERSION[0:3])) + "".join(VERSION[3:])
__license__ = "Apache Version 2"
__download_url__ = "https://github.com/downloads/briandconnelly/PlateTools/PlateTools-%s.tar.gz" % (__version__)

import PlateTools.formats
