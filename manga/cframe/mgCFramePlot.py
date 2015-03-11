#!/usr/bin/env python
# Filename : mgCFramePlot.py

import numpy
import argparse
from mgCFrameRead import CFrame as cf
from matplotlib import pyplot as plt

def HighlightOutput(width=40, symbol='#'):
    """
    Print something out in a highlight mode
    """
    bar = symbol * width
    print bar
    return None

def getFrameStr(name):
    if name[-1] is 'z':
        frameStr = name.replace(".fits.gz", "")
    elif name[-1] is 's':
        frameStr = name.replace(".fits", "")
    else:
        raise Exception("Wrong file format for mgCFrame !!")
    return frameStr

def plotSize19Avg1(wave, fluxarr, slitmap, name):
    slitindex = numpy.argwhere(slitmap['IFUSIZE'] == 19).reshape(2,21)
    fluxUse = fluxarr[slitindex[0]]
    fluxAvg = numpy.mean(fluxUse, axis=0)

    use = numpy.argwhere(numpy.logical_or(wave <= 5565, wave >=5590))
    minFlux = min(fluxAvg[use])
    maxFlux = max(fluxAvg[use])

    #------------------------------------------------------------
    # Plot the spectrum
    #------------------------------------------------------------
    plt.clf()
    ax1 = plt.axes()
    ax1.plot(wave, fluxAvg, '-', color='gray', label='error')
    ax1.set_xlabel(r'$\lambda (\AA)$')
    ax1.set_ylabel('Flux')

    ax1.set_xlim(min(wave), max(wave))
    ax1.set_ylim(minFlux, maxFlux)

    plt.savefig(name + "_ifusize19a_avg.png")

    return None

def plotSize19Avg2(wave, fluxarr, slitmap, name):
    slitindex = numpy.argwhere(slitmap['IFUSIZE'] == 19).reshape(2,21)
    fluxUse = fluxarr[slitindex[1]]
    fluxAvg = numpy.mean(fluxUse, axis=0)

    use = numpy.argwhere(numpy.logical_or(wave <= 5565, wave >=5590))
    minFlux = min(fluxAvg[use])
    maxFlux = max(fluxAvg[use])

    #------------------------------------------------------------
    # Plot the spectrum
    #------------------------------------------------------------
    plt.clf()
    ax2 = plt.axes()
    ax2.plot(wave, fluxAvg, '-', color='gray', label='error')
    ax2.set_xlabel(r'$\lambda (\AA)$')
    ax2.set_ylabel('Flux')

    ax2.set_xlim(min(wave), max(wave))
    ax2.set_ylim(minFlux, maxFlux)

    plt.savefig(name + "_ifusize19b_avg.png")

    return None

def plotSize37Avg(wave, fluxarr, slitmap, name):
    slitindex = numpy.argwhere(slitmap['IFUSIZE'] == 37).reshape(4,39)
    fluxUse = fluxarr[slitindex[3]]
    fluxAvg = numpy.mean(fluxUse, axis=0)

    use = numpy.argwhere(numpy.logical_or(wave <= 5565, wave >=5590))
    minFlux = min(fluxAvg[use])
    maxFlux = max(fluxAvg[use])

    #------------------------------------------------------------
    # Plot the spectrum
    #------------------------------------------------------------
    plt.clf()
    ax3 = plt.axes()
    ax3.plot(wave, fluxAvg, '-', color='gray', label='error')
    ax3.set_xlabel(r'$\lambda (\AA)$')
    ax3.set_ylabel('Flux')

    ax3.set_xlim(min(wave), max(wave))
    ax3.set_ylim(minFlux, maxFlux)

    plt.savefig(name + "_ifusize37_avg.png")

    return None


def main(name, fiber):

    #------------------------------------------------------------
    # Load in the mgCFrame file;
    # cframe is a HDUlist object from astropy.io.fits
    cframe = cf(name)
    wave    = cframe.getWave()
    fluxarr = cframe.getFlux()
    ferrarr = cframe.getError()
    fskyarr = cframe.getSky()
    maskarr = cframe.getMask()
    header  = cframe.getHeader()
    slitmap = cframe.getSlitMap()

    frameStr = getFrameStr(name)

    #plateRa, plateDec = cframe.getCenPos()
    plateType  = header["PLATETYP"]
    surveyMode = header["SRVYMODE"]
    exptime    = header["EXPTIME"]
    slitfile   = header["SLITFILE"]

    #------------------------------------------------------------
    # Check the index for fiber
    if (fiber < 1 or fiber >= cframe.nfiber):
        raise Exception("FIBER should be between 1 and %d"%(cframe.nfiber))
    flux = fluxarr[fiber-1]
    ferr = ferrarr[fiber-1]
    fsky = fskyarr[fiber-1]
    finfo = slitmap[fiber-1]
    #mask = maskarr[fiber-1]

    holetype   = finfo["HOLETYPE"]
    fiberid    = finfo["FIBERID"]
    fibertype  = finfo["FIBERTYPE"]
    plugstatus = finfo["PLUGSTATUS"]
    targettype = finfo["TARGETTYPE"]
    ifuinblock = finfo["IFUINBLOCK"]
    ifusize    = finfo["IFUSIZE"]
    ifudesign  = finfo["IFUDESIGN"]
    mangaid    = finfo["MANGAID"]

    use = numpy.argwhere(numpy.logical_or(wave <= 5565, wave >=5590))
    minFlux = min(min(flux[use]), min(fsky[use]))
    maxFlux = max(max(flux[use]), max(fsky[use]))

    #------------------------------------------------------------
    # Meng's plot
    #------------------------------------------------------------
    plotSize19Avg1(wave, fluxarr, slitmap, frameStr)
    plotSize19Avg2(wave, fluxarr, slitmap, frameStr)
    # TODO : Check the reshape part
    plotSize37Avg(wave, fluxarr, slitmap, frameStr)

    #------------------------------------------------------------
    # Display some useful information about thefile
    HighlightOutput(width=60)
    print " Read in mgCFrame file: " + name
    print "   This PLATE has %d fibers on it"%(cframe.nfiber)
    print "   Each spectrum has %d data points in it"%(cframe.npixel)
    print "   The wavelength range is %7.1d to %7.1d Angstrom"%(cframe.minwave,
                                                          cframe.maxwave)
    HighlightOutput(width=60)
    #print "   The central (RA,DEC) of the plate is (%d,%d)"%(plateRa,plateDec)
    print "   The type of this plate is : " + plateType
    print "   The survey mode for this plate is : " + surveyMode
    print "   The exposure time for this plate is : %6.1d seconds"%(exptime)
    print "   The associated slitmap name is :" + slitfile
    HighlightOutput(width=60)
    HighlightOutput(width=60)
    print " Will plot the spectra for fiber : " + str(fiber)
    HighlightOutput(width=60)

    #------------------------------------------------------------
    # Plot the spectrum
    #------------------------------------------------------------
    plt.clf()
    ax = plt.axes()
    ax.plot(wave, ferr, '-', color='gray', label='error')
    ax.plot(wave, fsky, '-', color='cyan', label='sky')
    ax.plot(wave, flux, '-k', label='spectrum')
    ax.set_title('File = ' + name + ', Fiber = %(fiber)i' % locals())

    ax.text(0.03, 0.95, "Holetype: " + holetype, size=12,
        ha='left', va='top', transform=ax.transAxes)
    ax.text(0.03, 0.90, "TargetType: " + targettype, size=12,
        ha='left', va='top', transform=ax.transAxes)
    ax.text(0.03, 0.85, "IFUinBlock: %d"%ifuinblock, size=12,
        ha='left', va='top', transform=ax.transAxes)
    ax.text(0.03, 0.80, "IFUSize: %d"%ifusize, size=12,
        ha='left', va='top', transform=ax.transAxes)
    ax.text(0.03, 0.75, "IFUDesign: %d"%ifudesign, size=12,
        ha='left', va='top', transform=ax.transAxes)
    ax.text(0.03, 0.70, "MaNGAID: " + mangaid, size=12,
        ha='left', va='top', transform=ax.transAxes)

    ax.set_xlabel(r'$\lambda (\AA)$')
    ax.set_ylabel('Flux')

    ax.set_xlim(cframe.minwave, cframe.maxwave)
    ax.set_ylim(minFlux, maxFlux)

    plt.savefig(frameStr + "-%d.png"%(fiber))

    return None

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("name", help="Name or location of the mgCFrame file")
    parser.add_argument("fiber", type=int, help="Number of the fiber")
    args = parser.parse_args()

    main(args.name, args.fiber)
