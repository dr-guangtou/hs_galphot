import sqlcl
import os
import numpy as np
from astropy.coordinates import FK5, ICRS
from astropy import units as u
from astropy.io import fits

"""
Search for SDSS images around certain area, collect information about these
images, and download them

By Song Huang; 2014-10-28
"""

def download_cutout(ra, dec, size, grid=False, label=False, specobj=False,
                   invert=False, name=None, suffix=None, pix=0.396,
                    dirout=None):
    """
    Download a JPEG color cutout image
    """
    return None

def download_frame(run, rerun, camcol, field, band='r', dirout=None, silent=True):
    """
    Download corrected frame image from SDSS archive
    """
    sdss_sas = "http://data.sdss3.org/sas/dr10/boss/photoObj/frames/"

    fits_file = "frame-{}-{}-{}-{}.fits.bz2".format(band, str(run).zfill(6),
                                              str(camcol), str(field).zfill(4))
    link = sdss_sas + str(rerun) + "/" + str(run) + "/" + str(camcol) +"/" + \
            fits_file

    if dirout is not None:
        dirout = dirout.strip()
        if dirout[-1] is not '/':
            dirout += '/'
        if not os.path.exists(dirout):
            os.system("mkdir -p " + dirout)
        os.system("wget -P " + dirout + " " + link)
    else:
        dirout = ""
        os.system("wget " + link)
    if os.path.exists(dirout + fits_file):
        os.system("bunzip2 " + dirout + fits_file)
    else:
        raise Exception("Can not download file : " + link + " !!")

    if not silent:
        highlight_output("Download image: " + fits_file + " !")

    return fits_file

def ra_dec_distance(ra1,dec1,ra2,dec2, use_icrs=False, arcsec=False):
    # conversion between degree and radians
    pi  = 3.1415926536
    d2r = pi/180.00

    if use_icrs:
        coord1 = ICRS(ra=ra1, dec=dec1, unit=(u.degree,u.degree))
        coord2 = ICRS(ra=ra2, dec=dec2, unit=(u.degree,u.degree))
    else:
        coord1 = FK5(ra=ra1, dec=dec1, unit=(u.degree,u.degree))
        coord2 = FK5(ra=ra2, dec=dec2, unit=(u.degree,u.degree))
    # The separation is in unit of arcsecond or degree
    if arcsec:
        sep = coord1.separation(coord2).arcsecond
    else:
        sep = (coord1.separation(coord2).radian / d2r )

    return sep


def highlight_output(string, width=30, symbol='#'):
    """
    Print something out in a highlight mode
    """
    bar = symbol * width
    print bar
    print string
    print bar

    return None


def define_query(ra,dec,radius,use_all=False,extra=None):
    """
    Define the SQL search string based on input, can be used by sqlcl

    !!! NOTE !!!
    This is still a very rough search,
    """
    radius += 0.12
    if extra:
        radius += extra
    # convert into arcmin
    radius *= 60.0

    if (use_all):
        query = "SELECT fieldID, run, rerun, camcol, field, quality, score, " + \
                "ra, dec, primaryArea from Field " + \
                "WHERE dbo.fDistanceArcMinEq({0},{1}".format(str(ra),str(dec)) + \
                ",ra,dec) <= {0}".format(str(radius))
    else:
        query = "SELECT fieldID, run, rerun, camcol, field, quality, score, " + \
                "ra, dec, primaryArea from Field " + \
                "WHERE dbo.fDistanceArcMinEq({0},{1}".format(str(ra),str(dec)) + \
                ",ra,dec) <= {0}".format(str(radius)) + \
                " AND primaryArea > 0.0"

    return query


def remote_search_field(ra,dec,radius,use_all=False,extra=None):
    """
    Given the Ra, DEC for the center of the field, the radius of the search area
    in unit of degree, returns the information for the FIELDS that cover this
    search area.

    The remote search use SQLCL to access the SDSS on-line database;
    The search is not very accurate

    @use_all : If True, all fields, including the ones with primaryArea==0, will
    be returned.
    """
    if (ra<0.0) or (ra>360.0):
        raise Exception("RA should be between 0 and 360 degree!")
    if (dec<-90.0) or (dec>90.0):
        raise Exception("Dec should be between -90 and 90 degree!")
    if radius > 1.0:
        warning = "Radius is too large!! Be careful!! (Radius<1.0deg)"
        highlight_output(warning)

    query = define_query(ra,dec,radius,use_all=use_all,extra=extra)
    result = sqlcl.query(query).readlines()

    n_field = (len(result)-2)
    if n_field == 0:
        raise Exception("No useful field is returned!! Check!!")

    data = []
    for ii in result[2:]:
        line = ii.replace("\n"," ")
        data.append(line.split(','))
    result = np.recarray(data, dtype=[('fieldid', int), ('run', int),
                                      ('rerun', int), ('camcol', int),
                                      ('field', int), ('quality', str),
                                      ('score', float), ('ra', float),
                                      ('dec', float), ('primaryarea', float)])
    print result.shape

    return result

def local_search_field(ra,dec,radius,use_all=False,dircat=None,extra=None):
    """
    Given the Ra, DEC for the center of the field, the radius of the search area
    in unit of degree, returns the information for the FIELDS that cover this
    search area.

    @use_all : If True, all fields, including the ones with primaryArea==0, will
    be returned.
    """
    field_cat = 'sdss_field_sum_coord.fit'
    if dircat:
        field_cat += dircat
    if os.path.exists(field_cat):
        catalog = fits.open(field_cat)
        data = catalog[1].data
        if use_all is False:
            data = data[data["primaryArea"] > 0.0]
        highlight_output("Load in :" + field_cat)
    else:
        raise Exception("Can not find FIELD_SUM catalog: " + field_cat + '!!')

    radius += 0.12
    if extra:
        radius += extra
    # Convert into arcsec
    radius *= 3600.0

    cen0 = FK5(ra=data["ra"], dec=data["dec"], unit=(u.degree, u.degree))
    cen1 = FK5(ra=ra, dec=dec, unit=(u.degree, u.degree))
    overlap = data[cen0.separation(cen1).arcsec <= radius]

    return overlap

if __name__ == '__main__':

    #data = remote_search_field(143.0,15.0,0.3)
    #print data.ra, data.dec

    data = local_search_field(143.0,15.0,0.3)
    print data.shape

    run    = data["run"][0]
    rerun  = data["rerun"][0]
    camcol = data["camcol"][0]
    field  = data["field"][0]
    fits = download_frame(run,rerun,camcol,field,silent=False,dirout="temp")
