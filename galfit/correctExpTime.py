#!/usr/bin/env python
# encoding: utf-8
"""Convert Segmentation to Mask Image."""

import os
import argparse

from astropy.io import fits


def run(imgFile, expTime, hdu=0):
    """Find and read in the image."""
    if os.path.isfile(imgFile):
        imgData = fits.open(imgFile)[hdu].data
        imgHeader = fits.open(imgFile)[hdu].header
    else:
        raise Exception("# Can not find input image")

    """ Correct for the exposure time """
    imgData *= expTime
    imgHeader['EXPTIME'] = expTime

    """ Save a new file """
    newFits = imgFile.replace('.fits', '_cor.fits')
    hdu = fits.PrimaryHDU(imgData)
    hdu.header = imgHeader
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(newFits, clobber=True)

    return imgData


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("imgFile", help="Name of the image file")
    parser.add_argument("expTime", type=float, help="Exposure time")

    args = parser.parse_args()

    run(args.imgFile, args.expTime)
