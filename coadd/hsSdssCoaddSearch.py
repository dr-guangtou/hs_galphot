#!/usr/bin/env python

import sqlcl
import os
import numpy as np
from StringIO import StringIO
from astropy.coordinates import FK5, ICRS
from astropy import units as u
from astropy.io import fits

"""
Search for SDSS images around certain area, collect information about these
images, and download them

By Song Huang; 2014-10-28
"""

def download_cutout(ra, dec, size, grid=True, label=True, specobj=True,
                   invert=False, name=None, suffix=None, pix=0.396,
                   dirout=None):
    """
    Download a JPEG color cutout image
    """
    # Image option
    option = ""
    if grid:
        option += "G"
    if label:
        option += "L"
    if specobj:
        option += "S"
    if invert:
        option += "I"

    # Parameters for the image
    ra_str   = str(ra).strip()
    dec_str  = str(dec).strip()
    size_str = str(size).strip()
    pix_str  = str(pix).strip()

    # Name of the output image
    if name is not None:
        output = name.strip()
    else:
        output = ra_str + '_' + dec_str
    if suffix is not None:
        output = output + "_" + suffix.strip()
    output += ".jpeg"

    # Link to the image on SDSS server
    link = "'http://skyservice.pha.jhu.edu/DR10/ImgCutout/getjpeg.aspx?" + \
            "ra=" + ra_str + "&dec=" + dec_str + "&scale=" + \
            pix_str + "&width=" + size_str + "&height=" + size_str + \
            "&opt=" + option + "&query='"

    if dirout is not None:
        dirout = dirout.strip()
        if dirout[-1] is not '/':
            dirout += '/'
        if not os.path.exists(dirout):
            os.system("mkdir -p " + dirout)
    else:
        dirout = ""
    output = dirout + output
    os.system("wget " + link + " -O " + output)

    if not os.path.exists(output):
        raise Exception("Can not find dowloaded image :" + output + " !!")

    return output


def download_frame(run, rerun, camcol, field, band='r', dirout=None,
                   silent=True):
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


def get_frame_name(run, rerun, camcol, field, band='r'):
    """
    Just returns the frame file name
    """

    return "frame-{}-{}-{}-{}.fits".format(band, str(run).zfill(6), str(camcol),
                                           str(field).zfill(4))


def get_frame_link(run, rerun, camcol, field, band='r'):
    """
    Just returns the frame file link
    """

    sdss_sas = "http://data.sdss3.org/sas/dr10/boss/photoObj/frames/"
    fits_file = "frame-{}-{}-{}-{}.fits.bz2".format(band, str(run).zfill(6),
                                              str(camcol), str(field).zfill(4))
    link = sdss_sas + str(rerun) + "/" + str(run) + "/" + str(camcol) +"/" + \
            fits_file

    return link


def ra_dec_distance(ra1, dec1, ra2, dec2, use_icrs=False, arcsec=False):
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


def define_query(ra, dec, radius, use_all=False, extra=None):
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
                "ra, dec, ramin, ramax, decmin, decmax, ngalaxy, nstars, " + \
                "primaryArea from Field " + \
                "WHERE dbo.fDistanceArcMinEq({0},{1}".format(str(ra),str(dec)) + \
                ",ra,dec) <= {0}".format(str(radius))
    else:
        query = "SELECT fieldID, run, rerun, camcol, field, quality, score, " + \
                "ra, dec, ramin, ramax, decmin, decmax, ngalaxy, nstars, " + \
                "primaryArea from Field " + \
                "WHERE dbo.fDistanceArcMinEq({0},{1}".format(str(ra),str(dec)) + \
                ",ra,dec) <= {0}".format(str(radius)) + \
                " AND primaryArea > 0.0"

    return query

def photo_query(ra, dec, radius, band='r', use_all=True, extra=None):
    """
    Define the SQL search strings for photometry parameters of objects within
    certain distance with input RA, DEC

    !!! NOTE !!!
    This is still a very rough search
    TODO: Deal with Time Out Error
    """
    radius += 0.12
    if extra:
        radius += extra
    ra1, ra2 = (ra - radius), (ra + radius)
    dec1, dec2 = (dec - radius), (dec + radius)

    if (use_all):
        query = "SELECT objId, ra, dec, type, clean, nChild, petroR90_" + \
                band + ", psfMag_" + band + ", cModelMag_" + band + \
                ", expAB_" + band + ", expPhi_" + band + \
                ", devAB_" + band + ", devPhi_" + band + \
                " FROM PhotoObj WHERE ra between " + \
                "%10.6f and %10.6f AND dec between %10.6f and %10.6f " % (ra1, ra2, dec1, dec2)
    else:
        query = "SELECT objId, ra, dec, type, clean, nChild, petroR90_" + \
                band + ", psfMag_" + band + ", cModelMag_" + band + \
                ", expAB_" + band + ", expPhi_" + band + \
                ", devAB_" + band + ", devPhi_" + band + \
                " FROM PhotoObj WHERE ra between " + \
                "%10.6f and %10.6f AND dec between %10.6f and %10.6f AND mode = 1" % (ra1, ra2, dec1, dec2)

    return query


def photoobj_search(ra, dec, radius, band='r', use_all=True, extra=0.05):
    """
    Searh for photometric information of all objects within certain distance
    between the central RA, DEC

    """

    if (ra < 0.0) or (ra > 360.0):
        raise Exception("RA should be between 0 and 360 degree!")
    if (dec < -90.0) or (dec > 90.0):
        raise Exception("Dec should be between -90 and 90 degree!")
    if radius > 1.0:
        warning = "Radius is too large!! Be careful!! (Radius < 1.0deg)"
        highlight_output(warning)

    query = photo_query(ra, dec, radius, band=band, use_all = use_all,
                        extra=extra)
    result = sqlcl.query(query).readlines()

    n_field = (len(result) - 2)
    if n_field <= 0:
        raise Exception("No useful field is returned!! Check!!")
    else:
        result = result[2:]

    data = []
    for ii in result:
        line = ii.replace("\n", "")
        temp = np.genfromtxt(StringIO(line), delimiter=",", dtype=None)
        data.append(temp)

    dtype = [('objID', int), ('ra', float), ('dec', float), ('type', int),
             ('clean', int), ('nChild', int), ('petroR90', float),
             ('psfMag', float), ('cModelMag', float), ('expAB', float),
             ('expPhi', float), ('devAB', float), ('devPhi', float)]
    table = np.array(data, dtype)

    return table


def remote_search_field(ra, dec, radius, use_all=False, extra=0.05):
    """
    Given the Ra, DEC for the center of the field, the radius of the search area
    in unit of degree, returns the information for the FIELDS that cover this
    search area.

    The remote search use SQLCL to access the SDSS on-line database;
    The search is not very accurate

    @use_all : If True, all fields, including the ones with primaryArea==0, will
    be returned.
    """
    if (ra < 0.0) or (ra > 360.0):
        raise Exception("RA should be between 0 and 360 degree!")
    if (dec < -90.0) or (dec > 90.0):
        raise Exception("Dec should be between -90 and 90 degree!")
    if radius > 1.0:
        warning = "Radius is too large!! Be careful!! (Radius < 1.0deg)"
        highlight_output(warning)

    query = define_query(ra, dec, radius, use_all = use_all, extra=extra)
    result = sqlcl.query(query).readlines()

    n_field = (len(result) - 2)
    if n_field <= 0:
        raise Exception("No useful field is returned!! Check!!")
    else:
        result = result[2:]

    data = []
    for ii in result:
        line = ii.replace("\n"," ")
        temp = np.genfromtxt(StringIO(line), delimiter=",", dtype=None)
        data.append(temp)

    dtype=[('fieldID', int), ('run', int), ('rerun', int), ('camcol', int),
           ('field', int), ('quality', int), ('score', float), ('ra', float),
           ('dec', float), ('raMin', float), ('raMax', float),
           ('decMin', float), ('decMax', float), ('nGalaxy', int),
           ('nStars', float), ('primaryArea', float)]

    result = np.array(data, dtype=dtype)

    return result


def local_search_field(ra, dec, radius, use_all=False, dircat=None, extra=None):
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
    result = data[cen0.separation(cen1).arcsec <= radius]

    return result


def dirty_coadd(imgList, outImg, ra, dec, radius, name=None):
    """
    Using SWarp to generate a quick-look quality coadd image
    """
    # TODO

    return None


def result_to_csv(results, listname, dirout=None ):
    """
    Output the search results to a ASCII file
    """
    #  TODO : Very very clumsy method, should be improved
    comma = " , "
    lineb = " \n"
    n_field = results.shape[0]

    if dirout:
        listname = dirout.strip() + listname.strip()
    else:
        listname = listname.strip()
    csv = open(listname, "w")

    # Header line
    csv.write("#RUN , RERUN , CAMCOL , FIELD , RA , DEC , RA_MIN , RA_MAX ," + \
             " DEC_MIN , DEC_MAX , SCORE , PRIMARYAREA , NGALAXY , NSTARS ," + \
             " FIELDID \n")

    for ii in range(n_field):
        str_1 = "%6i" % (results["run"][ii])
        str_2 = "%6i" % (results["rerun"][ii])
        str_3 = "%6i" % (results["camcol"][ii])
        str_4 = "%6i" % (results["field"][ii])
        str_5 = "%10.6d" % (results["ra"][ii])
        str_6 = "%10.6d" % (results["dec"][ii])
        str_7 = "%10.6d" % (results["raMin"][ii])
        str_8 = "%10.6d" % (results["raMax"][ii])
        str_9 = "%10.6d" % (results["decMin"][ii])
        str_10 = "%10.6d" % (results["decMax"][ii])
        str_11 = "%8.2d" % (results["score"][ii])
        str_12 = "%8.2d" % (results["primaryArea"][ii])
        str_13 = "%5i" % (results["nGalaxy"][ii])
        str_14 = "%5i" % (results["nStars"][ii])
        str_15 = "%13i" % (results["fieldID"][ii])
        csv.write(str_1 + comma + str_2 + comma + str_3 + comma + str_4 + comma +
                 str_5 + comma + str_6 + comma + str_7 + comma + str_8 + comma +
                 str_9 + comma + str_10 + comma + str_11 + comma + str_12 +
                 comma + str_13 + comma + str_14 + comma + str_15 + lineb)
    csv.close()

    return None


def search_and_down(ra, dec, radius, band='r', use_all=False, extra=None,
                    name=None, jpeg=False, remote=False):
    """
    Shortcut: Search and Download
    """

    # Ra and Dec of the center of the field
    ra_str  = str(ra).strip()
    dec_str = str(dec).strip()
    # Radius for the search, in degree
    rad_str = "%3.1d" % radius
    rad_str = rad_str.strip()

    # Check if band belongs to [u,g,r,i,z]
    sdss_filters = ['u', 'g', 'r', 'i', 'z']
    band = band.strip()
    if band not in sdss_filters:
        raise Exception("The band should belong to u, g, r, i, z")

    # The name of the main object
    if name is not None:
        dirout = name
    else:
        dirout = ra_str + "_" + dec_str

    # The output list of overlapped fields
    listname = dirout + "_" + rad_str
    if use_all:
        listname += "_all.csv"
    else:
        listname += "_pri.csv"

    # Search for the overlapped field
    if remote is True:
        results = remote_search_field(ra, dec, radius, use_all=use_all,
                                      extra=extra)
    else:
        results = local_search_field(ra, dec, radius, use_all=use_all,
                                     extra=extra)
    result_to_csv(results, listname, dirout=dirout)

    # Number of the returned fields
    n_field = results.shape[0]

    frame_list = []
    for ii in range(n_field):
        frame = download_frame(results["run"][ii], results["rerun"][ii],
                               results["camcol"][ii], results["field"][ii],
                               dirout=dirout, silent=False, band=band)
        frame_list.append(frame)

    if jpeg:
        jpeg_1 = download_cutout(ra, dec, 4096, name=name, suffix="large",
                                    dirout=dirout, grid=True, label=True,
                                    specobj=True)
        jpeg_2 = download_cutout(ra, dec, 1024, name=name, suffix="small",
                                    dirout=dirout, grid=True, label=True,
                                    specobj=True)

    return frame_list


if __name__ == '__main__':

    ra  = 226.62209
    dec = 55.76272
    rad = 0.20
    name = "NGC5866"
    band = "r"

    #data = remote_search_field(143.0,15.0,0.3)
    #print data.ra, data.dec

    #data = local_search_field(143.0,15.0,0.05)
    #print data.shape[0]
    #local_result2csv(data, "test.csv", dirout=None )

    #frame_list = local_search_and_down(ra, dec, rad, band=band, use_all=True,
    #                               extra=None, name=name, jpeg=True)
    #print len(frame_list)

