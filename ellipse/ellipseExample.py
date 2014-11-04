#!/usr/bin/env python
import os
import os.path

from pyraf import iraf

"""
To get familiar with Ellipse:

1. Check the help file for Ellipse, controlpar, samplepar, magpar, geompar
> ecl
> stsdas.analysis.isophote
> help ellipse
> help controlpar

2. See the examples on this page:
http://www.ast.uct.ac.za/~sarblyth/TullyFisher/ellipseEg/EllipseEg.html

3. Read the relevant section of Li, Ho, et al. 2011(CGS-II), and try to run
   ellipse in interactive mode on any data
"""

# Define the name of the input and output file
inputImg = "/home/song/work/ellipse/NGC1600_r.fit"
outBin = inputImg.replace(".fit", "_ellipse_1.bin")
outTab = inputImg.replace(".fit", "_ellipse_1.tab")
outCdf = inputImg.replace(".fit", "_ellipse_1.cdf")
# TODO: Check the .pl mask file, which should be
# inputMsk = inputImg.replace(".fit", ".pl")

# Call the STSDAS.ANALYSIS.ISOPHOTE package
iraf.stsdas()
iraf.analysis()
iraf.isophote()

# Define parameters for the ellipse run
# 1. Initial guess of the central X, Y (need to be as accurate as possible)
iraf.ellipse.geompar.x0 = 460.526
iraf.ellipse.geompar.y0 = 464.399
# 2. Initial guess of the ellipticity and PA of the first ISOPHOTE
#    Do not need to be very accurate, unless you want to fix them for all
#    isophotes and only derive surface brightness
iraf.ellipse.geompar.ellip0 = 0.6003035
iraf.ellipse.geompar.pa0 = -12.10127
# 3. Initial radius for ellipse fitting (The major axis length of the first
#    elliptical isophote); Can not be too small, and can not be too large
iraf.ellipse.geompar.sma0 = 40.48682917785644
# 4. The minimum and maximum radius for the ellipse fitting
iraf.ellipse.geompar.minsma = 0.5571857376098632
iraf.ellipse.geompar.maxsma = 94.98832999420166
# 5. Parameters about the stepsize during the fitting.
#    Unless you know what you what, normally should use log-stepsize instead of
#    linear one; and step=0.05 will generate more isophotes than step=0.1, but
#    may not help if you want a robust surface brightness profile.
iraf.ellipse.geompar.linear = "no"
iraf.ellipse.geompar.step = 0.1
# 6. Do you want to allow the ellipse to decide the galaxy center during the
#    fitting.  In general, it's a good idea to turn this on.  If the center you
#    provide is accurate enough, ELlipse results will not deviate from it.
iraf.ellipse.geompar.recenter = "yes"
# 7. The next three parameters control the behavior of the fit
#    hcenter = yes/no : Do all the isophotes have the same central X, Y?
#    hellip  = yes/no : Do all the isophotes have the same ellipticity?
#    hpa     = yes/no : Do all the isophotes have the same position angle?
# Based on our experience, the formal Ellipse fitting should be done in three
# separate runs
#    1) hcenter=no, hellip=no, hpa=no : Give Ellipse the total freedom to fit
#       the isophotes; And take the median/mean central X,Y from inner N
#       isophotes, then use these X,Y as the center of the galaxy
#    2) hcenter=yes, hellip=no, hpa=yes : Hold the central X, Y to the
#       previously determined ones; Let the ellipticity and position angle to be
#       free, then extract an appropriate average ellipticity and PA from this
#       run
#    3) hcenter=yes, hellip=yes, hpa=yes : Hold the center, and hold the
#       ellipticity and PA to the average values decided from the previous run.
#       Just extracted an robust surface brightness profile using the average
#       geometry
iraf.ellipse.controlpar.hcenter = "no"
iraf.ellipse.controlpar.hellip = "no"
iraf.ellipse.controlpar.hpa = "no"
# 8. Parameters about the iterations
#    minit/maxit: minimun and maximum number of the iterations
iraf.ellipse.controlpar.minit = 10
iraf.ellipse.controlpar.maxit = 100
# 9. Threshold for the object locator algorithm
#    By lowering this value, the locator become less strict.
iraf.ellipse.controlpar.olthresh = 1.00000
# 10. Make sure the Interactive Mode is turned off
iraf.ellipse.controlpar.interactive = "no"

# Check and remove outputs from the previous Ellipse run, or Ellipse will report
# error (Quite stupid!)
if os.path.exists(outBin):
	os.remove(outBin)
if os.path.exists(outTab):
	os.remove(outTab)
if os.path.exists(outCdf):
	os.remove(outCdf)

# Start the fitting
iraf.ellipse(input=inputImg, output=outTab)
# TODO: Demonstrate the direct photometry mode using input catalog
# inBin = input_bin_file
## The inBin is a Binary result from previous Ellipse run, and the isophote
## stored in it will overwrite all the above settings. Ellipse will simply
## extract surface brightness profile using these isophote instead of doing any
## fitting
# iraf.ellipse(input=inputImg, output=outBin, inellip=inBin)

# The Ellipse output is a binary table file, which is very hard to deal with
#  "Dump" it into a nice ASCII table
iraf.tdump(table=outBin, datafile=outTab, cdfile=outCdf)
os.remove(outCdf)
