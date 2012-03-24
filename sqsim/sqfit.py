#! /usr/bin/env python
# -*- python -*-

import sys, astimage, bgfit
from bgfit import fitknown

stpath = sys.argv[1]
image = astimage.open (sys.argv[2], 'r')

if len (sys.argv) == 4:
    residpath = sys.argv[3]
else:
    residpath = None

fitknown (image, stpath, residpath, 1, 0, sys.stdout.write,
          fixpos=False, fixshape=False)
