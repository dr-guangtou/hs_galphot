#!/usr/bin/env python
# encoding: utf-8
"""
Collection of useful tools .

 * Cosmology related procedures
 * Galactic extinction related
 * Angle normalization
 * Coordinate related shortcuts

"""

from __future__ import (absolute_import, division,
                        unicode_literals)
import os
import sys
import numbers
import collections

# Numpy & Scipy
import numpy as np
import scipy.ndimage
import scipy.interpolate
try:
    from scipy.stats import scoreatpercentile
except:
    scoreatpercentile = False
# Astropy
import astropy.units as u
from astropy.utils.misc import isiterable
from astropy.coordinates import SkyCoord


SIGMA1 = 0.3173
SIGMA2 = 0.0455
SIGMA3 = 0.0027


"""
Angle related functions:

Strongly based on: https://github.com/phn/angles/blob/master/angles.py
by Prasanth Nair

"""


def rad2deg(rad):
    """Convert radians into degrees."""
    return (rad * 180.0 / np.pi)


def deg2rad(deg):
    """Convert degrees into radians."""
    return (deg * np.pi / 180.0)


def hr2deg(deg):
    """Convert degrees into hours."""
    return (deg * (24.0 / 360.0))


def deg2hr(hr):
    """Convert hours into degrees."""
    return (hr * 15.0)


def normAngle(num, lower=0, upper=360, b=False):
    """
    Normalize number to range [lower, upper) or [lower, upper].

    Parameters
    ----------
    num : float
        The number to be normalized.
    lower : int
        Lower limit of range. Default is 0.
    upper : int
        Upper limit of range. Default is 360.
    b : bool
        Type of normalization. Default is False. See notes.
    Returns
    -------
    n : float
        A number in the range [lower, upper) or [lower, upper].
    Raises
    ------
    ValueError
      If lower >= upper.
    Notes
    -----
    If the keyword `b == False`, then the normalization is done in the
    following way. Consider the numbers to be arranged in a circle,
    with the lower and upper ends sitting on top of each other. Moving
    past one limit, takes the number into the beginning of the other
    end. For example, if range is [0 - 360), then 361 becomes 1 and 360
    becomes 0. Negative numbers move from higher to lower numbers. So,
    -1 normalized to [0 - 360) becomes 359.
    If the keyword `b == True`, then the given number is considered to
    "bounce" between the two limits. So, -91 normalized to [-90, 90],
    becomes -89, instead of 89. In this case the range is [lower,
    upper]. This code is based on the function `fmt_delta` of `TPM`.
    Range must be symmetric about 0 or lower == 0.
    Examples
    --------
    >>> normalize(-270,-180,180)
    90.0
    >>> import math
    >>> math.degrees(normalize(-2*math.pi,-math.pi,math.pi))
    0.0
    >>> normalize(-180, -180, 180)
    -180.0
    >>> normalize(180, -180, 180)
    -180.0
    >>> normalize(180, -180, 180, b=True)
    180.0
    >>> normalize(181,-180,180)
    -179.0
    >>> normalize(181, -180, 180, b=True)
    179.0
    >>> normalize(-180,0,360)
    180.0
    >>> normalize(36,0,24)
    12.0
    >>> normalize(368.5,-180,180)
    8.5
    >>> normalize(-100, -90, 90)
    80.0
    >>> normalize(-100, -90, 90, b=True)
    -80.0
    >>> normalize(100, -90, 90, b=True)
    80.0
    >>> normalize(181, -90, 90, b=True)
    -1.0
    >>> normalize(270, -90, 90, b=True)
    -90.0
    >>> normalize(271, -90, 90, b=True)
    -89.0
    """
    from math import floor, ceil
    # abs(num + upper) and abs(num - lower) are needed, instead of
    # abs(num), since the lower and upper limits need not be 0. We need
    # to add half size of the range, so that the final result is lower +
    # <value> or upper - <value>, respectively.
    res = num
    if not b:
        if lower >= upper:
            raise ValueError("Invalid lower and upper limits: (%s, %s)" %
                             (lower, upper))

        res = num
        if num > upper or num == lower:
            num = lower + abs(num + upper) % (abs(lower) + abs(upper))
        if num < lower or num == upper:
            num = upper - abs(num - lower) % (abs(lower) + abs(upper))

        res = lower if num == upper else num
    else:
        total_length = abs(lower) + abs(upper)
        if num < -total_length:
            num += ceil(num / (-2 * total_length)) * 2 * total_length
        if num > total_length:
            num -= floor(num / (2 * total_length)) * 2 * total_length
        if num > upper:
            num = total_length - num
        if num < lower:
            num = -total_length - num

        res = num

    res *= 1.0  # Make all numbers float, to be consistent

    return res


"""
Coordinate related shortcuts

    * Convert from (ra, dec) to (l, b)
    * Conversion between ICRS and FK5
"""


def radec2lb(ra, dec, radian=False, FK5=False):
    """
    Convert (ra, dec) into Galactic coordinate (l, b).

    Parameters
    ----------
    ra : float or list or array
        RA Coordinates in degree
    dec : float or list or array
        DEC Coordinates in degree

    Returns
    -------
    l : float or list or array
    b : float or list or array
    """
    """ See if the input is number or array"""
    if not (isiterable(ra) or isiterable(dec)):
        returnScalar = True
        if not FK5:
            raDec = [SkyCoord(ra, dec, frame='icrs', unit='deg')]
        else:
            raDec = [SkyCoord(ra, dec, frame='fk5', unit='deg')]
    else:
        returnScalar = False
        if not FK5:
            raDec = [SkyCoord(ra, dec, frame='icrs', unit='deg')
                     for rrr, ddd in zip(ra, dec)]
        else:
            raDec = [SkyCoord(ra, dec, frame='fk5', unit='deg')
                     for rrr, ddd in zip(ra, dec)]

    """ Convert to galactic coordinates
        Currently, coordinates do not support arrays; have to loop.
    """
    l = np.empty(len(raDec), dtype=np.float)
    b = np.empty(len(raDec), dtype=np.float)

    for ii, cc in enumerate(raDec):
        gg = cc.galactic
        # Hack to support both astropy v0.2.4 and v0.3.dev
        # TODO: remove this hack once v0.3 is out (and array-ify this
        # whole thing)
        if radian:
            l[ii] = gg.l.radian
            b[ii] = gg.b.radian
        else:
            l[ii] = gg.l.degree
            b[ii] = gg.b.degree

    if returnScalar:
        return l[0], b[0]
    else:
        return l, b


def icrs2fk5(ra, dec, radian=False):
    """Convert coordinates from ICRS to FK5 frame."""
    if not radian:
        raDec = SkyCoord(ra, dec, frame='icrs', unit='deg')
    else:
        raDec = SkyCoord(ra, dec, frame='icrs', unit='radian')

    raDecFK5 = raDec.transform_to('fk5')

    if not radian:
        return raDecFK5.ra.degree, raDecFK5.dec.degree
    else:
        return raDecFK5.ra.radian, raDecFK5.dec.radian


def fk52icrs(ra, dec, radian=False):
    """Convert coordinates from FK5 to ICRS frame."""
    if not radian:
        raDec = SkyCoord(ra, dec, frame='fk5', unit='deg')
    else:
        raDec = SkyCoord(ra, dec, frame='fk5', unit='radian')

    raDecFK5 = raDec.transform_to('icrs')

    if not radian:
        return raDecFK5.ra.degree, raDecFK5.dec.degree
    else:
        return raDecFK5.ra.radian, raDecFK5.dec.radian

"""
Image Visualization Related

    * zScale of image
"""


def zscale(img, contrast=0.25, samples=500):
    """
    Image scaling function.

    form http://hsca.ipmu.jp/hscsphinx/scripts/psfMosaic.html
    """
    ravel = img.ravel()
    ravel = ravel[np.isfinite(ravel)]

    if len(ravel) > samples:
        imsort = np.sort(np.random.choice(ravel, size=samples))
    else:
        imsort = np.sort(ravel)

    n = len(imsort)
    idx = np.arange(n)

    med = imsort[int(n / 2.0)]
    w = 0.25
    i_lo, i_hi = int((0.5-w)*n), int((0.5+w)*n)
    # BUG: Sometimes the polyfit could fail
    try:
        p = np.polyfit(idx[i_lo:i_hi], imsort[i_lo:i_hi], 1)
        slope, intercept = p
    except Exception:
        slope = 1.0

    z1 = med - (slope/contrast)*(n/2-n*w)
    z2 = med + (slope/contrast)*(n/2-n*w)

    return z1, z2


"""
Image Manipulation

    * Image resampling, like the rebin function in IDL
    * Image rotation
    * Image shift (sub-pixel accurarcy)
"""


def congrid(a, newdims, method='linear', centre=False, minusone=False):
    """
    From: http://wiki.scipy.org/Cookbook/Rebinning.

    Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).

    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    """
    if a.dtype not in [np.float64, np.float32]:
        a = np.cast[float](a)

    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array(a.shape)
    ndims = len(a.shape)
    if len(newdims) != ndims:
        print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
        return None
    newdims = np.asarray(newdims, dtype=float)
    dimlist = []

    if method == 'neighbour':
        for i in range(ndims):
            base = np.indices(newdims)[i]
            dimlist.append((old[i] - m1) / (newdims[i] - m1) * (
                           base + ofs) - ofs)
        cd = np.array(dimlist).round().astype(int)
        newa = a[list(cd)]
        return newa

    elif method in ['nearest', 'linear']:
        # calculate new dims
        for i in range(ndims):
            base = np.arange(newdims[i])
            dimlist.append((old[i] - m1) / (newdims[i] - m1) * (
                           base + ofs) - ofs)
        # specify old dims
        olddims = [np.arange(i, dtype=np.float) for i in list(a.shape)]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d(olddims[-1], a, kind=method)
        newa = mint(dimlist[-1])

        trorder = [ndims - 1] + range(ndims - 1)
        for i in range(ndims - 2, -1, -1):
            newa = newa.transpose(trorder)

            mint = scipy.interpolate.interp1d(olddims[i], newa, kind=method)
            newa = mint(dimlist[i])

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose(trorder)

        return newa
    elif method in ['spline']:
        # oslices = [slice(0, j) for j in old]
        # oldcoords = np.ogrid[oslices]
        nslices = [slice(0, j) for j in list(newdims)]
        newcoords = np.mgrid[nslices]

        newcoords_dims = range(np.rank(newcoords))
        # make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs
        deltas = (np.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas
        newcoords_tr -= ofs
        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print "Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported."
        return None


"""
File Manipulation

    * Save numpy array to cPickle file format
    * Save numpy array to hickle/HDF5 format
    * Save numpy array to csv file format
"""


def saveToPickle(array, name):
    """Save a numpy array to a cPickle/Pickle format binary file."""
    try:
        import cPickle as pickle
    except:
        import pickle

    output = open(name, 'w')
    pickle.dump(array, output, protocol=2)
    output.close()


def saveToHickle(array, name):
    """Save a numpy array to a hickle/HDF5 format binary file."""
    try:
        import hickle
    except:
        raise Exception("### The Hickle package is required!")

    output = open(name, 'w')
    hickle.dump(array, output, protocol=2)
    output.close()


def saveToCSV(array, name):
    """
    Save a numpy array to a CSV file.

    Use the dtype.name as column name if possible
    """
    output = open(name, 'w')
    colNames = array.dtype.names
    output.write("#" + ', '.join(colNames) + '\n')
    for item in array:
        line = ''
        for i in range(0, len(colNames)-1):
            col = colNames[i]
            line += str(item[col]) + ' , '
        line += str(item[colNames[-1]]) + '\n'
        output.write(line)
    output.close()


def parseRegEllipse(regName):
    """
    Parse a DS9 .reg files.

    convert the Ellipse or Circle regions
    into arrays of parameters for ellipse:
    x, y, a, b, theta
    """
    if os.path.isfile(regName):
        raise Exception("### Can not find the .reg file!")
    # Parse the .reg file into lines
    lines = [line.strip() for line in open(regName, 'r')]
    # Coordinate type of this .reg file: e.g. 'image'
    coordType = lines[2].strip()
    # Parse each region
    regs = [reg.split(" ") for reg in lines[3:]]

    xc = []
    yc = []
    ra = []
    rb = []
    theta = []

    for reg in regs:
        if reg[0].strip() == 'ellipse' and len(reg) is 6:
            xc.append(float(reg[1]))
            yc.append(float(reg[2]))
            ra.append(float(reg[3]))
            rb.append(float(reg[4]))
            theta.append(float(reg[5]) * np.pi / 180.0)
        elif reg[0].strip() == 'circle' and len(reg) is 4:
            xc.append(float(reg[1]))
            yc.append(float(reg[2]))
            ra.append(float(reg[3]))
            rb.append(float(reg[3]))
            theta.append(0.0)

    xc = np.array(xc, dtype=np.float32)
    yc = np.array(yc, dtype=np.float32)
    ra = np.array(ra, dtype=np.float32)
    rb = np.array(rb, dtype=np.float32)
    theta = np.array(theta, dtype=np.float32)

    return xc, yc, ra, rb, theta, coordType


"""
Cosmology Related

    * Get luminosity distance at redshift=z
    * Get angular diameter distance at redshift = z
    * Get pixel scale at redshift = z
    * Get distance module at redshift = z
    * Get comoving volume at redshift = z
    * Get differential volume at redsfhit = z
    * Get the age of the Universe at redshift = z
    * Get the look-back time at redshift = z

"""


def cosmoDL(redshift, WMAP9=True, H0=69.3, Om0=0.287,
            Planck15=True, kpc=False):
    """
    Get the Luminosity Distance at redshift=z.

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    dl = cosmo.luminosity_distance(redshift)

    if not kpc:
        return dl.value
    else:
        return dl.to(u.kpc).value


def cosmoDA(redshift, WMAP9=True, H0=69.3, Om0=0.287,
            Planck15=True, kpc=False):
    """
    Get the Angular Diameter Distance at redshift=z.

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    da = cosmo.angular_diameter_distance(redshift)

    if not kpc:
        return da.value
    else:
        return da.to(u.kpc).value


def cosmoScale(redshift, WMAP9=True, H0=69.3, Om0=0.287,
               Planck15=True):
    """
    Get the Angular Scale (kpc/") at redshift=z.

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    scale = cosmo.kpc_proper_per_arcmin(redshift).to(u.kpc / u.arcsec)

    return scale.value


def cosmoDistMod(redshift, WMAP9=True, H0=69.3, Om0=0.287,
                 Planck15=False):
    """
    Get the Distance Module at redshift=z.

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    dm = cosmo.distmod(redshift)

    return dm.value


def cosmoComVol(redshift, WMAP9=True, H0=69.3, Om0=0.287,
                Planck15=False, Gpc=False):
    """
    Get the Comoving Volume at redshift=z.

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    v = cosmo.comoving_volume(redshift)

    if not Gpc:
        return v.value
    else:
        return v.to(u.Gpc).value


def cosmodVol(redshift, WMAP9=True, H0=69.3, Om0=0.287,
              Planck15=False):
    """
    Get the Differential Comoving Volume at redshift=z.

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    dv = cosmo.differential_comoving_volume(redshift)

    return dv.value


def cosmoAge(redshift, WMAP9=True, H0=69.3, Om0=0.287,
             Planck15=False, Myr=False):
    """
    Get the Age of the Universe at redshift=z.

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    age = cosmo.age(redshift)

    if not Myr:
        return age.value
    else:
        return age.to(u.Myr).value


def cosmoLookBack(redshift, WMAP9=True, H0=69.3, Om0=0.287,
                  Planck15=False, Myr=False):
    """
    Get the Look-back Time at redshift=z.

    This is simply a wrapper of astropy.cosmology
    The input redsfhit can be an array
    """
    if WMAP9:
        from astropy.cosmology import WMAP9 as cosmo
    elif Planck15:
        from astropy.cosmology import Planck15 as cosmo
    else:
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)

    lbt = cosmo.lookback_time(redshift)

    if not Myr:
        return lbt.value
    else:
        return lbt.to(u.Myr).value


"""
Galactic Extinction Related

"""


def getExtinction(ra, dec, a_lambda=None):
    """
    Estimate the Galactic extinction for HSC filters.

    Parameters:
        ra, dec : The input coordinates can be arrays
    """
    # First try mwdust from Jo Bovy
    try:
        import mwdust
        sfd = mwdust.SFD(sf10=True)

        from astropy.coordinates import SkyCoord
        coords = SkyCoord(ra, dec, frame='icrs', unit='deg')
        galactic = coords.galactic
        l, b = galactic.l, galactic.b
        ebv = sfd(l, b, 0)
    except ImportError:
        try:
            # Then try sncosmo
            from sncosmo import SFD98Map
            dustDir = os.environ.get('DUST_DIR')
            if (not os.path.isfile(os.path.join(dustDir,
                                   'SFD_dust_4096_ngp.fits'))) or (
                not os.path.isfile(os.path.join(dustDir,
                                   'SFD_dust_4096_sgp.fits'))):
                print('# DUST_DIR : %s' % dustDir)
                raise Exception("# Can not find the SFD dust map!")
            else:
                sfd = SFD98Map(dustDir)
                ebv = sfd.get_ebv((ra, dec))
        except ImportError:
            raise Exception("# Both mwdust and sncosmo are not available")
    if a_lambda is not None:
        return (ebv * a_lambda)
    else:
        return ebv


"""
Geometry Related
"""


def ellipDist(x, y, x0, y0, pa=0.0, q=0.9):
    """Distance to center in elliptical coordinate."""
    theta = (pa * np.pi / 180.0)

    distA = ((x - x0) * np.cos(theta) + (y - y0) * np.sin(theta)) ** 2.0
    distB = (((y - y0) * np.cos(theta) - (x - x0) * np.sin(theta)) / q) ** 2.0

    return np.sqrt(distA + distB)


"""
Weighted mean and median

Based on https://github.com/tinybike/weightedstats
"""


def weighted_mean(data, weights=None):
    """Calculate the weighted mean of a list."""
    if weights is None:
        return np.mean(data)
    total_weight = float(sum(weights))
    weights = [weight / total_weight for weight in weights]
    w_mean = 0
    for i, weight in enumerate(weights):
        w_mean += weight * data[i]
    return w_mean


def numpy_weighted_mean(data, weights=None):
    """Calculate the weighted mean of an array/list using numpy."""
    weights = np.array(weights).flatten() / float(sum(weights))
    return np.dot(np.array(data), weights)


def weighted_median(data, weights=None):
    """Calculate the weighted median of a list."""
    if weights is None:
        return np.median(data)
    midpoint = 0.5 * sum(weights)
    if any([j > midpoint for j in weights]):
        return data[weights.index(max(weights))]
    if any([j > 0 for j in weights]):
        sorted_data, sorted_weights = zip(*sorted(zip(data, weights)))
        cumulative_weight = 0
        below_midpoint_index = 0
        while cumulative_weight <= midpoint:
            below_midpoint_index += 1
            cumulative_weight += sorted_weights[below_midpoint_index-1]
        cumulative_weight -= sorted_weights[below_midpoint_index-1]
        if cumulative_weight - midpoint < sys.float_info.epsilon:
            bounds = sorted_data[below_midpoint_index-2:below_midpoint_index]
            return sum(bounds) / float(len(bounds))
        return sorted_data[below_midpoint_index-1]


def numpy_weighted_median(data, weights=None):
    """Calculate the weighted median of an array/list using numpy."""
    if weights is None:
        return np.median(np.array(data).flatten())
    data, weights = np.array(data).flatten(), np.array(weights).flatten()
    if any(weights > 0):
        sorted_data, sorted_weights = map(np.array,
                                          zip(*sorted(zip(data, weights))))
        midpoint = 0.5 * sum(sorted_weights)
        if any(weights > midpoint):
            return (data[weights == np.max(weights)])[0]
        cumulative_weight = np.cumsum(sorted_weights)
        below_midpoint_index = np.where(cumulative_weight <= midpoint)[0][-1]
        if (cumulative_weight[below_midpoint_index] -
                midpoint) < sys.float_info.epsilon:
            return np.mean(sorted_data[below_midpoint_index:
                                       below_midpoint_index+2])
        return sorted_data[below_midpoint_index+1]


"""
PolyNomial Fitting
"""


def polyFit(x, y, order=4):
    """Fit polynomial."""
    if len(x) != len(y):
        raise Exception("### X and Y should have the same size")
    coefficients = np.polyfit(x, y, order)
    polynomial = np.poly1d(coefficients)
    fit = polynomial(x)

    return fit


"""
Random color map from Photoutils
"""


def random_cmap(ncolors=256, background_color='black', random_state=None):
    """
    Generate a matplotlib colormap consisting of random (muted) colors.

    A random colormap is very useful for plotting segmentation images.
    Parameters
    ----------
    ncolors : int, optional
        The number of colors in the colormap.  For use with segmentation
        images, ``ncolors`` should be set to the number of labels.  The
        default is 256.
    background_color : str, optional
        The name of the background (first) color in the colormap.  Valid
        colors names are defined by ``matplotlib.colors.cnames``.  The
        default is ``'black'``.
    random_state : int or `~numpy.random.RandomState`, optional
        The pseudo-random number generator state used for random
        sampling.  Separate function calls with the same
        ``random_state`` will generate the same colormap.
    Returns
    -------
    cmap : `matplotlib.colors.Colormap`
        The matplotlib colormap with random colors.
    """
    from matplotlib import colors

    prng = check_random_state(random_state)
    h = prng.uniform(low=0.0, high=1.0, size=ncolors)
    s = prng.uniform(low=0.2, high=0.7, size=ncolors)
    v = prng.uniform(low=0.5, high=1.0, size=ncolors)
    hsv = np.transpose(np.array([h, s, v]))
    rgb = colors.hsv_to_rgb(hsv)

    if background_color is not None:
        if background_color not in colors.cnames:
            raise ValueError('"{0}" is not a valid background color '
                             'name'.format(background_color))
        rgb[0] = colors.hex2color(colors.cnames[background_color])

    return colors.ListedColormap(rgb)


def check_random_state(seed):
    """
    Turn seed into a `numpy.random.RandomState` instance.

    Parameters
    ----------
    seed : `None`, int, or `numpy.random.RandomState`
        If ``seed`` is `None`, return the `~numpy.random.RandomState`
        singleton used by ``numpy.random``.  If ``seed`` is an `int`,
        return a new `~numpy.random.RandomState` instance seeded with
        ``seed``.  If ``seed`` is already a `~numpy.random.RandomState`,
        return it.  Otherwise raise ``ValueError``.

    Returns
    -------
    random_state : `numpy.random.RandomState`
        RandomState object.

    Notes
    -----
    This routine is from scikit-learn.  See
    http://scikit-learn.org/stable/developers/utilities.html#validation-tools.
    """
    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (numbers.Integral, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    raise ValueError('%r cannot be used to seed a numpy.random.RandomState'
                     ' instance' % seed)


"""
Boostrap Resampling Method for Confidence Interval.

See: https://github.com/ptweir/pyresampling/blob/master/bootstrap.py
And: http://peterthomasweir.blogspot.jp/2014/03/
     statistics-based-on-resampling-in.html
"""


def confidence_interval_1d(A, alpha=SIGMA1, metric=np.mean,
                           numResamples=10000, interpolate=True):
    """Calculate confidence interval along one dimensional array."""
    if not isinstance(alpha, collections.Iterable):
        alpha = np.array([alpha])

    N = len(A)
    resampleInds = np.random.randint(0, N, (numResamples, N))
    metricOfResampled = metric(A[resampleInds], axis=-1)
    confidenceInterval = np.zeros(2*len(alpha), dtype='float')

    if interpolate:
        for thisAlphaInd, thisAlpha in enumerate(alpha):
            percenPos = (thisAlpha * 100 / 2.0)
            samplePos = scoreatpercentile(metricOfResampled, percenPos)
            confidenceInterval[2*thisAlphaInd] = samplePos
            percenNeg = (100 - thisAlpha * 100 / 2.0)
            sampleNeg = scoreatpercentile(metricOfResampled, percenNeg)
            confidenceInterval[2*thisAlphaInd+1] = sampleNeg
    else:
        sortedMetricOfResampled = np.sort(metricOfResampled)
        for thisAlphaInd, thisAlpha in enumerate(alpha):
            percenPos = int(round(thisAlpha*numResamples / 2.0))
            samplePos = sortedMetricOfResampled[percenPos]
            confidenceInterval[2*thisAlphaInd] = samplePos
            percenNeg = int(round(numResamples -
                                  (thisAlpha * numResamples / 2.0)))
            sampleNeg = sortedMetricOfResampled[percenNeg]
            confidenceInterval[2*thisAlphaInd+1] = sampleNeg
    return confidenceInterval


def ma_confidence_interval_1d(A, alpha=.05, metric=np.mean,
                              numResamples=1000, interpolate=True):
    """
    Confidence interval for 1-D array.

    Parameters
    """
    A = np.ma.masked_invalid(A, copy=True)
    A = A.compressed()
    confidenceInterval = confidence_interval_1d(A, alpha, metric,
                                                numResamples, interpolate)
    return confidenceInterval


def confidence_interval(A, axis=None, alpha=.05, metric=np.mean,
                        numResamples=1000, interpolate=True):
    """
    Bootstrap confidence interval.

    Return the bootstrap confidence interval of an array or along an axis
    ignoring NaNs and masked elements.

    Parameters
    ----------
    A : array_like
        Array containing numbers whose confidence interval is desired.
    axis : int, optional
        Axis along which the confidence interval is computed.
        The default is to compute the confidence interval of the flattened
        array.
    alpha: float or array, optional
        confidence level of confidence interval. 100.0*(1-alpha) percent
        confidence interval will be returned.
        If length-n array, n confidence intervals will be computed
        The default is .05
    metric : numpy function, optional
        metric to calculate confidence interval for.
        The default is numpy.mean
    numResamples : int, optional
        number of bootstrap samples. The default is 10000.
    interpolate: bool, optional
        uses scipy.stats.scoreatpercentile to interpolate between
        bootstrap samples if alpha*numResamples/2.0 is not integer.
        The default is True

    Returns
    -------
    confidenceInterval : ndarray
    An array with the same shape as `A`, with the specified axis replaced by
    one twice the length of the alpha

    If `A` is a 0-d array, or if axis is None, a length-2 ndarray is returned.
    """
    if interpolate is True and scoreatpercentile is False:
        print("need scipy to interpolate between values")
        interpolate = False
    A = A.copy()
    if axis is None:
        A = A.ravel()
        outA = ma_confidence_interval_1d(A, alpha, metric, numResamples,
                                         interpolate)
    else:
        outA = np.apply_along_axis(ma_confidence_interval_1d, axis, A, alpha,
                                   metric, numResamples, interpolate)

    return outA


def songPlotSetup(ax, border=4.5,
                  xlabel=30, ylabel=30,
                  majorTickL=12, minorTickL=8,
                  majorTickW=4.5, minorTickW=4.0):
    """Setup the format of the figure."""
    # Axes setup
    #  Minor Ticks on
    ax.minorticks_on()
    #  Axes Thickness
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(border)

    #  Tick Label Size
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(xlabel)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(ylabel)

    #  Tick Length and Width
    ax.tick_params('both', length=majorTickL, width=majorTickW,
                   which='major')
    ax.tick_params('both', length=minorTickL, width=minorTickW,
                   which='minor')

    return ax


"""
* HSC Catalog related functions.
"""


def removeIsNullCol(cat, output=None, catHdu=1,
                    string='isnull'):
    """Remove the xxx_isnull columns from the catalog."""
    from astropy.table import Table
    if not os.path.isfile(cat):
        raise Exception("Can not find catalog: %s" % cat)

    if output is None:
        output = cat.replace('.fits', '_clean.fits')

    data = Table.read(cat, format='fits')
    print "Reading the data"
    colnames = data.colnames

    colRemove = [col for col in colnames if string in col]

    data.remove_columns(colRemove)

    data.write(output, format='fits', overwrite=True)

    print "Saving new data"
