#! /usr/bin/env python
# -*- python -*-

import sys, os.path, omega as om, numpy as np, srctable
import omega.pango

omega.pango.setFont (size=8)

## quickutil: arraygrower
#- snippet: arraygrower.py (2012 Feb 27)
#- SHA1: 8ae43ac24e7ea0fb6ee2cc1047cab1588433a7ec
class ArrayGrower (object):
    __slots__ = 'dtype ncols chunkSize _nextIdx _arr'.split ()

    def __init__ (self, ncols, dtype=None, chunkSize=128):
        if dtype is None:
            import numpy as np
            dtype = np.float

        self.dtype = dtype
        self.ncols = ncols
        self.chunkSize = chunkSize
        self.clear ()


    def clear (self):
        self._nextIdx = 0
        self._arr = None
        return self


    def __len__ (self):
        return self._nextIdx


    def addLine (self, line):
        from numpy import asarray, ndarray

        line = asarray (line, dtype=self.dtype)
        if line.size != self.ncols:
            raise ValueError ('line is wrong size')

        if self._arr is None:
            self._arr = ndarray ((self.chunkSize, self.ncols),
                                 dtype=self.dtype)
        elif self._arr.shape[0] <= self._nextIdx:
            self._arr.resize ((self._arr.shape[0] + self.chunkSize,
                               self.ncols))

        self._arr[self._nextIdx] = line
        self._nextIdx += 1
        return self


    def add (self, *args):
        self.addLine (args)
        return self


    def finish (self):
        if self._arr is None:
            from numpy import ndarray
            ret = ndarray ((0, self.ncols), dtype=self.dtype)
        else:
            self._arr.resize ((self._nextIdx, self.ncols))
            ret = self._arr

        self.clear ()
        return ret
## end

ag = ArrayGrower (3)
data = []
worstscore = 0
worstinfo = None

for d in sys.argv[1:]:
    for p in os.listdir (d):
        if not p.endswith ('.fk'):
            continue

        for s in srctable.readStream (os.path.join (d, p))[2]:
            ag.add (s.totflux_initial, s.totflux, s.totflux_uc)

            sc = np.abs ((s.totflux - s.totflux_initial) / s.totflux_uc)
            if sc > worstscore:
                worstscore = sc
                worstinfo = d, p, s.ra, s.dec

    true, fitted, uc = ag.finish ().T
    w = np.where (fitted > 0) ; y = np.log10 (fitted[w] / true[w])
    #y = fitted / true - 1
    #y = (fitted - true) / uc
    from scipy.stats import scoreatpercentile as sap
    print 'bounds:', d, sap (y, 2), sap (y, 98)
    print '25/75/iqr:', d, sap (y, 25), sap (y, 75), sap (y, 75) - sap (y, 25)
    print 'median:', d, np.median (y)
    data.append ((d, y))

import astutil
print 'worst:', worstinfo[0], worstinfo[1], worstscore, astutil.fmtradec (worstinfo[2],
                                                                          worstinfo[3])

tagmap = {'flux-ne': 'No effects', 'flux-sq': 'ATA-scale squint'}

p = om.RectPlot ()
for tag, y in data:
    usetag = tagmap.get (tag, tag)
    values, edges = np.histogram (y, bins=21, #range=(-0.0925, 0.0125),
                                  normed=True)
    edges = edges[:-1] # convert to old-style
    fp = om.rect.ContinuousSteppedPainter (keyText=usetag)
    fp.setFloats (edges, values)
    p.add (fp)
    p.rebound (False, True)

p.magicAxisPainters ('ltbr')
p.tpainter.paintLabels = p.rpainter.paintLabels = False
p.setLabels ('r = log<sub>10</sub>(f<sub>fit</sub>/f<sub>true</sub>)',
             'dP/dr')
p.show ()
#p.save ('sqfluxsim.pdf', style=om.styles.ColorOnWhiteVector (),
#        margins=(4,4,4,4), dims=(420,380))
