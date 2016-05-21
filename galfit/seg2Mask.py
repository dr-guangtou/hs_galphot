#!/usr/bin/env python
# encoding: utf-8
"""Convert Segmentation to Mask Image."""

import os
import argparse
import collections

import numpy as np
from astropy.io import fits

# Scipy
import scipy.ndimage as ndimage


def run(segFile, sigma=6.0, mskThr=0.01, objs=None, removeCen=True):
    """Convert segmentation map to mask image."""
    segFile = args.segFile
    mskFile = segFile.replace('.fits', '_msk.fits')
    if not os.path.isfile(segFile):
        raise Exception("## Can not find the segmentation image : %s" %
                        segFile)

    """Load in the segmentation image"""
    segImg = fits.open(segFile)[0].data
    xSize, ySize = segImg.shape

    if removeCen:
        """Find out the value of the central pixel"""
        cenVal = segImg[int(xSize/2), int(ySize/2)]
        print "# Segment %d is the central object" % cenVal
        """Clear the central object"""
        segImg[segImg == cenVal] = 0

    if (objs is not None) and isinstance(objs, collections.Iterable):
        for obj in objs:
            try:
                segImg[segImg == int(obj)] = 0
            except Exception:
                print "# Can not find object: %d" % obj
                continue

    """Making a mask array"""
    segImg[segImg > 0] = 1
    segImg = segImg.astype(int)

    # Convolve the mask image with a gaussian kernel
    mskConv = ndimage.gaussian_filter(segImg * 1.0, sigma=sigma, order=0)
    mskBool = mskConv > mskThr
    mskInt = mskBool.astype(np.int16)

    """Save the array to fits file"""
    hdu = fits.PrimaryHDU(mskInt)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(mskFile, clobber=True)

    return segImg


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("segFile", help="Name of the segmentation file")
    parser.add_argument("-s", "--sigma", dest="sigma",
                        type=float, default=2.0,
                        help="Sigma of the Gaussian kernel for convolution")
    parser.add_argument("-t", "--threshold", dest="mskThr",
                        type=float, default=0.02,
                        help="Lower value cut")

    args = parser.parse_args()

    run(args.segFile, sigma=args.sigma, mskThr=args.mskThr)
