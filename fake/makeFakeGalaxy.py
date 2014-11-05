# File: makeFakeGalaxy.py

import numpy as np
import galsim
import pyfits as fits

def arrayToGSObj(imgArr, scale=1.0, norm=False):
    if norm:
        gsObj = galsim.InterpolatedImage(galsim.image.Image(imgArr), scale=scale,
                                         normalization="flux")
    else:
        gsObj = galsim.InterpolatedImage(galsim.image.Image(imgArr), scale=scale)
    return gsObj


def galSimDrawImage(galObj, size=0, method="auto", addPoisson=False):

    # Generate an "Image" object for the model
    # TODO: The default setting is "auto"; However, when there is trunctation,
    # it will use the "real_space" method, which is very time consuming; Using
    # "fft" method will be less accurate, but much faster
    # TODO: Also, simply apply truncated profile will not reduce the size of the
    # image, we can manually truncate the image by draw it into a smaller image
    if size > 0:
        imgTemp = galsim.image.Image(size, size)
        galImg = galObj.drawImage(imgTemp, method=method)
    else:
        galImg = galObj.drawImage(method=method)

    # Just an option for test
    if addPoisson:
        galImg.addNoise(galsim.PoissonNoise())

    # Return the Numpy array version of the image
    return galImg.array


def galSimConvolve(galObj, psfObj, size=0, method="auto", returnObj=False):
    """
    Just do convolution using GalSim
    Make sure the inputs are both GalSim GSObj
    The galaxy model should be the first one, and the PSF object is the second
    one; Returns a imgArr or GSObj
    """
    outObj = galsim.Convolve([galObj, psfObj])

    if returnObj:
        return outObj
    else:
        outArr = galSimDrawImage(galObj, size=size, method=method)
        return outArr


def galSimAdd(galObjList, size=0, method="auto", returnArr=False):
    """
    Just add a list of GSObjs together using GalSim
    Make sure all elements in the input list are GSObjs
    """
    if len(galObjList) < 2:
        raise Exception("Should be more than one GSObjs to add !")

    outObj = galsim.Add(galObjList)

    if returnArr:
        outArr = galSimDrawImage(outObj, size=size, method=method)
        return outArr
    else:
        return outObj


def plotFakeGalaxy(galObj, galID=None, size=0, addPoisson=False):

    import matplotlib.pyplot as plt

    if galID is None:
        outPNG = 'fake_galaxy.png'
    else:
        outPNG = 'fake_galaxy_%i.png' % galID

    plt.figure(1, figsize=(8,8))
    # Use "fft" just to be fast
    plt.imshow(np.arcsinh(galSimDrawImage(galObj, size=size, method="fft",
                                         addPoisson=addPoisson)))
    plt.savefig(outPNG)


def galSimFakeSersic(flux, gal, psfImage=None, scaleRad=False, returnObj=True,
                     expAll=False, devAll=False, plotFake=False, trunc=0,
                     psfNorm=False, drawMethod="auto", addPoisson=False):

    # TODO: The real input parser should be here
    # TODO: So...numpy.float32 is different from just float data type
    galID     = int(gal["ID"])
    nSersic   = float(gal["sersic_n"])
    reffPix   = float(gal["reff_pix"])
    axisRatio = float(gal["b_a"])
    posAng    = float(gal["theta"])

    # TODO: Decide how to do truncation or simply cut the image
    # Right now, we truncate at trunc * reffPix
    if trunc > 0:
        trunc = int(trunc * reffPix)

    # Make sure Sersic index is not too large
    if nSersic > 6.0:
        raise ValueError("Sersic index is too large! Should be <= 6.0")
    # Check the axisRatio value
    if axisRatio <= 0.15:
        raise ValueError("Axis Ratio is too small! Should be >= 0.15")

    # Make the Sersic model based on flux, re, and Sersic index
    if nSersic == 1.0 or expAll:
        if scaleRad:
            serObj = galsim.Exponential(scale_radius=reffPix)
        else:
            serObj = galsim.Exponential(half_light_radius=reffPix)
    elif nSersic == 4.0 or devAll:
        serObj = galsim.DeVaucouleurs(half_light_radius=reffPix, trunc=trunc)
    else:
        serObj = galsim.Sersic(nSersic, half_light_radius=reffPix, trunc=trunc)

    # If necessary, apply the Axis Ratio (q=b/a) using the Shear method
    if axisRatio < 1.0:
        serObj = serObj.shear(q=axisRatio, beta=0.0*galsim.degrees)

    # If necessary, apply the Position Angle (theta) using the Rotate method
    if posAng != 0.0 or posAng != 180.0:
        serObj = serObj.rotate(posAng*galsim.degrees)

    # Convolve the Sersic model using the provided PSF image
    if psfImage is not None:
        # Convert the PSF Image Array into a GalSim Object
        if psfNorm:
            psfObj = arrayToGSObj(psfImage, norm=True)
        else:
            psfObj = arrayToGSObj(psfImage)
        serFinal = galsim.Convolve([serObj, psfObj])
    else:
        serFinal = serObj

    # Pass the flux to the object
    # TODO: Pass the flux here or before convolution
    serFinal = serFinal.withFlux(float(flux))
    #serObj = serObj.withFlux(float(flux))
    print " With Flux : ", serFinal.getFlux()

    # Make a PNG figure of the fake galaxy to check if everything is Ok
    # TODO: For test, should be removed later
    if plotFake:
        plotFakeGalaxy(serFinal, galID=galID, size=trunc)

    # Now, by default, the function will just return the GSObj
    if returnObj:
        return serFinal
    else:
        if trunc > 0:
            galArray = galSimDrawImage(serFinal, size=trunc, method=drawMethod,
                                       addPoisson=addPoisson)
        else:
            galArray = galSimDrawImage(serFinal, method=drawMethod,
                                       addPoisson=addPoisson)
        return galArray



def galSimFakeDoubleSersic(flux1, reffPix1, nSersic1, axisRatio1, posAng1,
                           flux2, reffPix2, nSersic2, axisRatio2, posAng2,
                           psfImage, noConvole=False, addPoisson=False,
                           plotFake=False, galID=None):

    serModel1 = galSimFakeSersic(flux1, reffPix1, nSersic1, axisRatio1, posAng1,
                                 psfImage, noConvole=True, returnObj=True)
    serModel2 = galSimFakeSersic(flux2, reffPix2, nSersic2, axisRatio2, posAng2,
                                 psfImage, noConvole=True, returnObj=True)
    serDouble = serModel1 + serModel2

    # Convert the PSF Image Array into a GalSim Object
    psfObj = arrayToGSObj(psfImage)

    # Convolve the Sersic model using the provided PSF image
    # TODO: Make sure we understand the normalization
    if not noConvole:
        serFinal = galsim.Convolve([serDouble, psfObj])
    else:
        serFinal = serDouble

    # Generate an "Image" object for the convolved model
    serImg = serFinal.drawImage()

    # Return the Numpy array version of the image
    galArray = serImg.array

    # Make a PNG figure of the fake galaxy to check if everything is Ok
    # TODO: Should be removed later
    if plotFake:
        plotFakeGalaxy(galArray, galID=galID)

    return galArray


def testMakeFake(galList, asciiTab=False):

    # Make a fake Gaussian PSF
    psfGaussian = galsim.Gaussian(fwhm=2.0)
    psfImage    = psfGaussian.drawImage().array

    if asciiTab:
        galData = np.loadtxt(galList, dtype=[('ID','int'),
                                             ('mag','float'),
                                             ('sersic_n','float'),
                                             ('reff_pix','float'),
                                             ('b_a','float'),
                                             ('theta','float')])
    else:
        galData = fits.open(galList)[1].data

    # Test SingleSersic
    for igal, gal in enumerate(galData):

        flux = 10.0 ** ((27.0 - gal['mag']) / 2.5)

        print '---------------------------------'
        print " Input Flux : ", flux
        print " Input Parameters : ", gal["sersic_n"], gal["reff_pix"], gal["b_a"], gal["theta"]

        galArray = galSimFakeSersic(flux, gal, psfImage=psfImage, plotFake=True,
                                    returnObj=False, psfNorm=True, trunc=10.0,
                                    drawMethod="auto")

        print " Output Flux : ", np.sum(galArray)
        print " Shape of the Output Array : ", galArray.shape
        print '---------------------------------'

    # Test DoubleSersic
    flux1 = 10.0 ** ((27.0 - galData['mag'][0]) / 2.5)
    flux2 = 10.0 ** ((27.0 - galData['mag'][1]) / 2.5)
    re1, re2 = galData['reff_pix'][0], galData['reff_pix'][1]
    ns1, ns2 = galData['sersic_n'][0], galData['sersic_n'][1]
    ba1, ba2 = galData['b_a'][0], galData['b_a'][1]
    pa1, pa2 = galData['theta'][0], galData['theta'][1]

    #doubleArray = galSimFakeDoubleSersic(flux1, re1, ns1, ba1, pa1,
    #                                     flux2, re2, ns2, ba2, pa2,
    #                                     psfImage, plotFake=True, galID=4)

    #print (flux1 + flux2), np.sum(doubleArray)
    #print doubleArray.shape

