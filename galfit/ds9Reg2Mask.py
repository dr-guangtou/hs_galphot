#!/usr/bin/env python
# encoding: utf-8
"""Convert DS9 regions into mask image."""

from __future__ import division

import os
import copy
import numpy as np
import argparse

from astropy.io import fits
import cubehelix
# For high-contrast image
cmap1 = cubehelix.cmap(start=0.5, rot=-0.8, gamma=1.0,
                       minSat=1.2, maxSat=1.2,
                       minLight=0.0, maxLight=1.0)
cmap1.set_bad('k', 1.)

# import coaddCutoutPrepare as ccp
# Matplotlib related
import matplotlib as mpl
mpl.use('Agg')
mpl.rcParams['figure.figsize'] = 12, 10
mpl.rcParams['xtick.major.size'] = 8.0
mpl.rcParams['xtick.major.width'] = 1.5
mpl.rcParams['xtick.minor.size'] = 4.0
mpl.rcParams['xtick.minor.width'] = 1.5
mpl.rcParams['ytick.major.size'] = 8.0
mpl.rcParams['ytick.major.width'] = 1.5
mpl.rcParams['ytick.minor.size'] = 4.0
mpl.rcParams['ytick.minor.width'] = 1.5
mpl.rc('axes', linewidth=2)
import matplotlib.pyplot as plt
plt.ioff()


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

    med = imsort[int(n/2)]
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


def showSEPImage(image, contrast=0.2, size=10, cmap=cmap1,
                 title='Image', pngName='sep.png', titleInside=True,
                 ellList1=None, ellList2=None, ellList3=None,
                 ellColor1='b', ellColor2='r', ellColor3='g',
                 ell1=None, ell2=None, ell3=None, ellColor4='k',
                 ax=None, mask=None, mskAlpha=0.4):
    """
    Visualization of the results.

    Parameters:
    """
    fig = plt.figure(figsize=(size, size))
    fig.subplots_adjust(hspace=0.0, wspace=0.0,
                        bottom=0.08, left=0.08,
                        top=0.92, right=0.98)
    ax = fig.add_axes([0.000, 0.002, 0.996, 0.996])
    fontsize = 16
    ax.minorticks_on()

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)

    ax.set_title(title, fontsize=28, fontweight='bold', color='r')
    if not titleInside:
        ax.title.set_position((0.5, 1.01))
    else:
        ax.title.set_position((0.5, 0.90))

    imcopy = copy.deepcopy(image)
    imin, imax = zscale(imcopy, contrast=contrast, samples=500)

    ax.imshow(np.arcsinh(imcopy), interpolation="none",
              vmin=imin, vmax=imax, cmap=cmap, origin='lower')

    if mask is not None:
        # imcopy[mask > 0] = np.nan
        ax.imshow(mask, interpolation="none", vmin=0, vmax=1, origin='lower',
                  alpha=mskAlpha, cmap='gray_r')

    if ellList1 is not None:
        for e in ellList1:
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_alpha(0.8)
            e.set_edgecolor(ellColor1)
            e.set_facecolor('none')
            e.set_linewidth(1.5)

    if ellList2 is not None:
        for e in ellList2:
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_alpha(0.8)
            e.set_edgecolor(ellColor2)
            e.set_facecolor('none')
            e.set_linewidth(1.5)

    if ellList3 is not None:
        for e in ellList3:
            ax.add_artist(e)
            e.set_clip_box(ax.bbox)
            e.set_alpha(0.8)
            e.set_edgecolor(ellColor3)
            e.set_facecolor('none')
            e.set_linewidth(1.5)

    if ell1 is not None:
        ax.add_artist(ell1)
        ell1.set_clip_box(ax.bbox)
        ell1.set_alpha(0.8)
        ell1.set_edgecolor('r')
        ell1.set_facecolor('none')
        ell1.set_linewidth(2.0)
        ell1.set_linestyle('dashed')

    if ell2 is not None:
        ax.add_artist(ell2)
        ell2.set_clip_box(ax.bbox)
        ell2.set_alpha(0.8)
        ell2.set_edgecolor(ellColor4)
        ell2.set_facecolor('none')
        ell2.set_linewidth(2.5)
        ell2.set_linestyle('dashed')

    if ell3 is not None:
        ax.add_artist(ell3)
        ell3.set_clip_box(ax.bbox)
        ell3.set_alpha(0.8)
        ell3.set_edgecolor(ellColor4)
        ell3.set_facecolor('none')
        ell3.set_linewidth(2.5)
        ell3.set_linestyle('dashed')

    fig.savefig(pngName)
    plt.close(fig)


def saveFits(img, fitsName, head=None, clobber=True):
    """
    Save an image to FITS file.

    Parameters:
    """
    imgHdu = fits.PrimaryHDU(img)
    if head is not None:
        imgHdu.header = head
    imgHdu.writeto(fitsName, clobber=clobber)


def reg2Mask(imgFile, regFile, mskFile=None, hdu=0, show=False,
             save=True, imgHead=None, reverse=False):
    """
    Mask out the regions in a DS9 region file.

    Parameters:
    """
    try:
        import pyregion
    except Exception:
        raise Exception("### Please have pyregion installed first")

    if imgHead is None:
        if not os.path.isfile(imgFile):
            raise Exception("### Can not find the Image: %s" % imgFile)
        else:
            img = fits.open(imgFile)[0].data
            head = fits.open(imgFile)[0].header
    else:
        img = imgFile
        head = imgHead

    if not os.path.isfile(regFile):
        raise Exception("### Can not find the Region file: %s" % regFile)
    else:
        reg = pyregion.open(regFile).as_imagecoord(head)

    imgX, imgY = img.shape
    regMask = reg.get_mask(shape=(imgX, imgY))
    if not reverse:
        intMask = regMask.astype(int)
    else:
        intMask = np.invert(regMask).astype(int)

    if save:
        if mskFile is None:
            if not reverse:
                mskFile = regFile.replace('.reg', '_msk.fits')
            else:
                mskFile = regFile.replace('.reg', '_invmsk.fits')
        saveFits(intMask, mskFile, head=head, clobber=True)

    if show:
        pngName = mskFile.replace('.fits', '.png')
        showSEPImage(img, pngName=pngName, mask=intMask, cmap=cmap1,
                     title=mskFile.replace('.fits', ''))

    return intMask


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("imgFile", help="Name of the image")
    parser.add_argument("regFile", help="Name of the DS9 region")
    parser.add_argument("-m", "--mskFile", dest='mskFile', default=None,
                        help="Name of the DS9 region")
    parser.add_argument("-s", "--show", dest='show', action="store_true",
                        default=False,
                        help="Whether to show the mask in PNG figure")
    parser.add_argument("--hdu", dest='hdu', default=0,
                        help="The HDU to be used", type=int)
    parser.add_argument("-r", "--reverse", dest='reverse',
                        action="store_true", default=False,
                        help="Reverse the mask")
    args = parser.parse_args()

    reg2Mask(args.imgFile, args.regFile, mskFile=args.mskFile, hdu=args.hdu,
             show=args.show, reverse=args.reverse)
