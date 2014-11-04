#!/usr/bin/env python
""" @ hs_run_ellipse.py
Various utilities for running IRAF.ELLIPSE to extract 1-D photometric
information.
"""
import os
import sys
import numpy
from pyraf import iraf
from iraf import tables, ttools

#
def fits_to_pl( in_image, out_image=None ):
    """
    Convert a fits image to .pl file for IRAF
    Useful for feeding a object mask image to ELLIPSE
    """
    try:
        if out_image is None:
            input_ext = os.path.splitext( in_image )[1]
            out_image = in_image.replace( input_ext, '.pl' )
        """
        Delete the old file, if existed
        """
        if os.path.exists( out_image ):
            os.remove( out_image )
        """ imcopy """
        iraf.imcopy( in_image, out_image )
        """ Return the path of the output """
        return out_image
    except:
        print "<Error>: ", sys.exc_info()[1]
        raise

#
def bin_to_asc( in_table, out_table=None, cdfile=None, pfile=None ):
    """
    Convert the IRAF binary table to ASCII table
    Using IRAF.TABLES.TTOOLS.TDUMP
    """
    try:
        if out_table is None:
            input_ext = os.path.splitext( in_table )[1]
            out_table = in_table.replace( input_ext, '.asc' )
        """
        Delete the old table, if existed
        """
        if os.path.exists( out_table ):
            os.remove( out_table )
        if cdfile is None:
            cdfile = 'asc_column'
        if pfile  is None:
            pfile  = 'asc_header'
        """ tdump """
        tables.ttools.tdump( table=in_table, datafile=out_table,
                            cdfile=cdfile, pfile=pfile )
        return out_table
    except:
        print "<Error>: ", sys.exc_info()[1]
        raise

#
def parse_ellipse_out( bin_table, NaN='NaN', scale=None, magzpt=None ):
    """
    Parse the output binary table from IRAF.ELLIPSE

    @param NaN      The string to replace INDEF in the output table
    @param scale    Pixel scale in unit of arcsec/pixel to convert the value
                    of radius from pixel to arcse
    @param magzpt   Magnitude zeropoint to convert intensity into magnitude
    """
    if os.path.exists( bin_table ):
        asc_tab = bin_to_asc( bin_table)
        """ Replace the INDEF with NaN """
        NaN = str( NaN )
        os.system( 'sed -i "s/INDEF/' + NaN + '/g" ' + asc_tab )
        """ Read in the table using Astropy.IO.ASCII """
        from astropy.io import ascii
        data = ascii.read( asc_tab )

    else:
        raise NameError( '<ERROR>: Can not find input table!!' )


"""
Sanity checks and Test of the functions
"""

""" Test fits_to_pl """
#im_fits = '/home/hs/astro1/hsc_test/example/cosmos_example/temp/test1a.fits'
#out_pl  = '/home/hs/astro1/hsc_test/example/cosmos_example/temp/test1a_2.pl'
#im_1 = fits_to_pl( im_fits, out_image=out_pl )
#im_2 = fits_to_pl( im_fits )
#im_fits = '/home/hs/astro1/hsc_test/example/cosmos_example/temp/test1b.fits'
#im_3 = fits_to_pl( im_fits )

""" Test bin_to_asc """
#tab_bin = '/home/hs/astro1/hsc_test/example/cosmos_example/temp/test1.bin'
#tab_out = '/home/hs/astro1/hsc_test/example/cosmos_example/temp/test1a.asc'
#tab_1 = bin_to_asc( tab_bin )
#tab_2 = bin_to_asc( tab_bin, out_table=tab_out )

