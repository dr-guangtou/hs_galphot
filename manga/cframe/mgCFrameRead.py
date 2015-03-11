#!/usr/bin/env python
# Filename : mgCFrameRead.py

import numpy
import os
from astropy.io import fits
from matplotlib import pyplot as plt

"""
Data model for mgCFrame file:

HDU[0] = Empty, only used to store header
HDU[1] = Flux [NPixels,NFiber], in Unit of 10^(-17) erg/s/cm^2/Ang/Fiber
HDU[2] = InverseVariance of the flux; [NPixels,NFiber][=1/sigma^2]
HDU[3] = Pixel Mask [NPixels,NFiber]
HDU[4] = Wavelength [NPixels] [Vacuum]
HDU[5] = WavelengthDispersion [NPixels,NFiber]
HDU[6] = SlitMap
HDU[7] = SkyFlux: [NPixels,Nfiber]

Filename: mgCFrame-00177379-LIN.fits.gz
No.    Name         Type      Cards   Dimensions   Format
0    PRIMARY     PrimaryHDU     173   ()
1    FLUX        ImageHDU        12   (6732, 1423)   float32
2    IVAR        ImageHDU        12   (6732, 1423)   float32
3    MASK        ImageHDU        12   (6732, 1423)   int32
4    WAVE        ImageHDU        11   (6732,)      float64
5    DISP        ImageHDU        12   (6732, 1423)   float64
6    SLITMAP     BinTableHDU    175   1423R x 33C   [12A, J, J, 8A, J, J, 8A, J, J, J, E, J, K, J, 5A, J, 3A, J, J, J, D, D, D, D, D, D, E, E, E, J, J, 5E, 9A]
7    SKY         ImageHDU        12   (6732, 1423)   float32
"""

class CFrame:

    def __init__(self, name):
        if not os.path.exists(name):
            raise Exception("Can not find the mgCFrame file: " + name)
        self.name = name.strip()
        self.hdulist = fits.open(self.name)
        self.npixel = -1
        self.nfiber = -1
        self.minwave = -1
        self.maxwave = -1

    def listHDU(self):
        self.hdulist.info()

    def listHeader(self):
        header = self.hdulist[0].header
        print header
        return header

    def getHeader(self):
        return self.hdulist[0].header

    def getWave(self):
        wave = self.hdulist[4].data
        self.minwave = wave.min()
        self.maxwave = wave.max()
        return wave

    def getFlux(self):
        flux = self.hdulist[1].data
        self.npixel = flux.shape[1]
        self.nfiber = flux.shape[0]
        return flux

    def getIvar(self):
        return self.hdulist[2].data

    def getError(self):
        return numpy.sqrt(1.0 / self.hdulist[2].data)

    def getMask(self):
        return self.hdulist[3].data

    def getWdisp(self):
        return self.hdulist[5].data

    def getSlitMap(self):
        return self.hdulist[6].data

    def getSky(self):
        return self.hdulist[7].data

    def getFluxOri(self):
        return (self.hdulist[1].data + self.hdulist[7].data)

