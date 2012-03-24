#! /usr/bin/env python
# -*- python -*-

"""
Simulate squinty data.
"""

import sys, numpy as np, miriad, progress
from numpy import add, dot, exp, multiply, square
from astutil import *
from mirtask import util

"""
Random note: the purest data are with makesqmag = lambda: 0, pbfreqdep
= False. In that case, the quality-limiting effect in the *UV* domain
is the 16-bit default MIRIAD storage format. In terms of imaging, it
is pixelization, if pxalignimg is not None. Of course fitting sources
and subtracting them in the UV plane should work equally well
regardless of the value of pxalignimg ...
"""

npsrc = 20
ngsrc = 0
makesqmag = lambda: np.random.triangular (50, 100, 500) * A2R
pbfreqdep = True
randseed = None
pxalignimg = None # '../sm-ne-casa.restored'

nant = 43
recmod = 20000
makesrcrad = lambda: np.random.triangular (0, 4000, 4000) * A2R
makesrcflux = lambda: np.random.uniform (1e-3, 3) #could do s**-1.5 distrib
groupdist = 300 * A2R
tooclosedist = 75 * A2R

twopii = (0+2j) * np.pi
gfac = -0.25 * np.pi**2 / np.log (2)
randangle = lambda: np.random.uniform (0, 2 * np.pi)

def pbalphas (freqs):
    sigmas = 3.5 / freqs * F2S * D2R
    return -0.5 * sigmas**-2


def basicsource (ra0, dec0, alpha, pximg=None):
    s = Source ()
    s.totflux = makesrcflux ()
    s.dec, s.ra = sphofs (dec0, ra0, makesrcrad (), randangle ())

    if pximg is not None:
        # Adjust to land precisely on a pixel if pixalign enabled
        y, x = pximg.topixel ([s.dec, s.ra])
        y = int (round (y))
        x = int (round (x))
        s.dec, s.ra = pximg.toworld ([y, x])

    s.basedist = sphdist (s.dec, s.ra, dec0, ra0)
    s.fadedflux = s.totflux * np.exp (alpha * s.basedist**2)

    l = np.sin (s.ra - ra0) * np.cos (s.dec)
    m = (np.sin (s.dec) * np.cos (dec0) -
         np.cos (s.ra - ra0) * np.cos (s.dec) * np.sin (dec0))
    n = np.sqrt (1 - l**2 - m**2)
    n -= 1 # fringe rot in effect; makes the work below easier
    s.lmn = np.asarray ([l, m, n])

    s.fitgroup = None
    s.delete = False

    return s


class Source (object):
    pass


class SquintSimulator (object):
    banner = 'squint simulator'


    def genSources (self, ra0, dec0, freq0):
        self.sources = []

        a = pbalphas (freq0)

        if pxalignimg is None:
            pximg = None
        else:
            import astimage
            pximg = astimage.open (pxalignimg, 'r').simple ()

        for i in xrange (npsrc):
            s = basicsource (ra0, dec0, a, pximg)
            self.sources.append (s)
            s.major = s.minor = s.pa = None

        for i in xrange (ngsrc):
            s = basicsource (ra0, dec0, a, pximg)
            self.sources.append (s)
            s.major = 300 * A2R
            s.minor = 80 * A2R
            s.pa = orientcen (randangle ())

        # Assign sources to fitgroups so that bgfit can handle very nearby
        # sources. We get rid of *extremely* close sources, on the grounds that
        # in practice we'd never identify them as two separate sources. Otherwise
        # we make sure that any two sources nearer than *groupdist* are in the
        # same fitgroup. This could have various pathological failure modes but
        # will be fine most of the time (and the only penalty to overgrouping is
        # a larger instantaneous fitting problem, I think).

        groupnum = 0

        for i in xrange (len (self.sources)):
            s1 = self.sources[i]
            if s1.delete:
                continue

            for j in xrange (i + 1, len (self.sources)):
                s2 = self.sources[j]
                if s2.delete:
                    continue

                d = sphdist (s1.dec, s1.ra, s2.dec, s2.ra)

                if d > groupdist:
                    continue
                if d < tooclosedist:
                    s2.delete = True
                    continue

                if s1.fitgroup is not None and s2.fitgroup is not None:
                    # Need to merge two existing fitgroups
                    newgroup = s1.fitgroup
                    oldgroup = s2.fitgroup
                    for s3 in self.sources:
                        if s3.fitgroup == oldgroup:
                            s3.fitgroup = newgroup
                elif s1.fitgroup is not None:
                    s2.fitgroup = s1.fitgroup
                elif s2.fitgroup is not None:
                    # This can happen, if sources A B and C appear in
                    # that order in the list, and A is near enough to
                    # C, and B is near enough to C, but A is not near
                    # enough to B.
                    s1.fitgroup = s2.fitgroup
                else:
                    s1.fitgroup = str (groupnum)
                    groupnum += 1
                    s2.fitgroup = s1.fitgroup

        self.sources = [s for s in self.sources if not s.delete]

        # Now we can output.

        import srctable, flatdb
        cols = srctable.columns ('ra dec totflux major minor pa')
        cols.append (flatdb.Column (name='fitgroup', kind=flatdb.K_STR, width=10))
        cols.append (srctable.floatcol ('fadedflux', 12, '%.5f'))
        flatdb.writeStreamedTable (sys.stdout.write, [], cols, self.sources)


    def genSquints (self):
        # Say that X's are all on-center, Y's are all offset
        self.squints = squints = {}

        for antnum in xrange (1, nant):
            squints[util.antpol2ap (antnum, util.FPOL_X)] = (0, 0)
            ymag = makesqmag ()
            ypa = randangle ()
            squints[util.antpol2ap (antnum, util.FPOL_Y)] = (ymag, ypa)

        for s in self.sources:
            s.hpbsqdists = {}


    def rewrite (self, ncorr, invis, outvis, **uvdargs):
        psources = [s for s in self.sources if s.major is None]
        gsources = [s for s in self.sources if s.major is not None]

        outhnd = outvis.open ('c')
        outhnd.setPreambleType ('uvw', 'time', 'baseline')

        first = True
        curhnd = None
        npol = prevnpol = 0
        neednpol = True
        npolvaried = False
        nproc = 0
        nrec = 0
        prevtime = None

        for inhnd, preamble, data, flags in invis.readLowlevel ('3', False, **uvdargs):
            if first:
                inhnd.copyItem (outhnd, 'history')
                outhnd.openHistory ()
                outhnd.writeHistory (self.banner)
                first = False

            if inhnd is not curhnd:
                inhnd.initVarsAsInput ('channel')
                outhnd.initVarsAsOutput (inhnd, 'channel')
                #outhnd.setCorrelationType ('r')
                curhnd = inhnd
                freqsupdated = True # assuming fixed freqs per dataset

            if freqsupdated:
                freqs = inhnd.getSkyFrequencies ()
                alphas = pbalphas (freqs)
                dwork = np.empty (freqs.shape, dtype=np.double)
                dwork2 = np.empty (freqs.shape, dtype=np.double)
                cwork = np.empty (freqs.shape, dtype=np.complex)
                freqsupdated = False

            if prevtime is None or preamble[3] != prevtime:
                # Need to determine new squint distances as PB rotates on the sky.
                lst = inhnd.getVarDouble ('lst')
                lat = inhnd.getVarDouble ('latitud')
                az, el = util.equToHorizon (inhnd.getVarDouble ('ra'),
                                            inhnd.getVarDouble ('dec'),
                                            lst, lat)

                for ap, (sqmag, sqpa) in self.squints.iteritems ():
                    oel, oaz = sphofs (el, az, sqmag, sqpa)
                    ora, odec = util.horizonToEqu (oaz, oel, lst, lat)

                    for s in self.sources:
                        # We take half of the squared distance because
                        # that's what's going to be relevant in the PB
                        # attenuation (half in exponent -> sqrt ->
                        # geometric mean)
                        s.hpbsqdists[ap] = 0.5 * sphdist (s.dec, s.ra, odec, ora)**2

                prevtime = preamble[3]

            if nrec % recmod == 0:
                print >>sys.stderr, '%5.1f%% (%d of %d)' % (100.* nproc / ncorr,
                                                            nproc, ncorr)

            if npol == 0:
                npol = inhnd.getNPol ()
                neednpol = True

            if neednpol:
                if npol != prevnpol:
                    outhnd.writeVarInt ('npol', npol)
                    npolvaried = npolvaried or (prevnpol != 0)
                    prevnpol = npol
                neednpol = False

            ap1, ap2 = util.mir2bp (inhnd, preamble)
            data.fill (0)
            flags.fill (1)

            for s in psources:
                # geometric mean of effective fluxes, attenuated by
                # (maybe) freq-dependent PB falloff.
                multiply (twopii * dot (s.lmn, preamble[:3]), freqs, cwork)
                if pbfreqdep:
                    multiply (s.hpbsqdists[ap1] + s.hpbsqdists[ap2], alphas, dwork)
                    add (dwork, cwork, cwork)
                else:
                    add ((s.hpbsqdists[ap1] + s.hpbsqdists[ap2]) * alphas[0], cwork, cwork)
                exp (cwork, cwork)
                multiply (s.totflux, cwork, cwork)
                data += cwork

            for s in gsources:
                # As above, with additional gaussian amplitude falloff term.
                cp = np.cos (s.pa)
                sp = np.sin (s.pa)
                g = gfac * ((s.minor * (preamble[0] * cp - preamble[1] * sp))**2 +
                            (s.major * (preamble[0] * sp + preamble[1] * cp))**2)
                square (freqs, dwork2)
                multiply (g, dwork2, dwork2)
                multiply (twopii * dot (s.lmn, preamble[:3]), freqs, cwork)
                if pbfreqdep:
                    multiply (s.hpbsqdists[ap1] + s.hpbsqdists[ap2], alphas, dwork)
                    add (dwork, cwork, cwork)
                else:
                    add ((s.hpbsqdists[ap1] + s.hpbsqdists[ap2]) * alphas[0], cwork, cwork)
                add (dwork2, cwork, cwork)
                exp (cwork, cwork)
                multiply (s.totflux, cwork, cwork)
                data += cwork

            inhnd.copyLineVars (outhnd)
            outhnd.writeVarInt ('pol', inhnd.getPol ())
            outhnd.write (preamble, data, flags)
            npol -= 1
            nrec += 1
            nproc += data.size

        if not npolvaried:
            outhnd.setScalarItem ('npol', np.int32, prevnpol)

        outhnd.closeHistory ()
        outhnd.close ()


if randseed is not None:
    np.random.seed (randseed)

vin = miriad.VisData (sys.argv[1])
vout = miriad.VisData (sys.argv[2])

tmphnd = vin.open ('rw')
# need process a record to get sky freqs
tmphnd.lowlevelRead (np.empty (5, dtype=np.float64),
                     np.empty (1024, dtype=np.complex64),
                     np.empty (1024, dtype=np.intc))
ra0 = tmphnd.getVarDouble ('ra') # very unsure about whether this or obsra is right
dec0 = tmphnd.getVarDouble ('dec') # need to be consistent with this and squint code
freq0 = tmphnd.getSkyFrequencies ()[0]
ncorr = tmphnd.getScalarItem ('ncorr')
tmphnd.close ()

ss = SquintSimulator ()
ss.genSources (ra0, dec0, freq0)
ss.genSquints ()
ss.rewrite (ncorr, vin, vout, select='-auto')

"""
Inefficient/readable point vis expression:
pb = alphas * (s.hpbsqdists[ap1] + s.hpbsqdists[ap2])
ph = twopii * dot (s.lmn, preamble[:3]) * freqs
data += s.totflux * np.exp (ph + pb)

Inefficient/readable gaussian vis expression:
pb = alphas * (s.hpbsqdists[ap1] + s.hpbsqdists[ap2])
ph = twopii * dot (s.lmn, preamble[:3]) * freqs
cp = np.cos (s.pa)
sp = np.sin (s.pa)
g = ((s.minor * (preamble[0] * cp - preamble[1] * sp) * freqs)**2 +
     (s.major * (preamble[0] * sp + preamble[1] * cp) * freqs)**2)
g *= gfac
data += s.totflux * np.exp (ph + pb + g)

Explicit expression for last direction cosine:
n = (np.sin (s.dec) * np.sin (dec0) +
     np.cos (s.ra - ra0) * np.cos (s.dec) * np.cos (dec0))
"""
