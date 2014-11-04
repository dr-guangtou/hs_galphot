## SDSS Image Coadd 

* Eventually we want to download all the SDSS single frame images around certain
   object, and generate the coadd image by ourselves.  We should: 
   - Take care of the PSF matching 
   - Be careful about the background matching 
   - Mask out problematic or suspicious regions: e.g. saturated stars 

### The First Step: Search and Download Images 

* We want to download all the images that covers a certain area around the 
    object we are interested in: 
    - We need to the RA, DEC of the galaxy (center) 
    - The size of the coadd image we want to build: a RADIUS in degree 
    - The filter we want to build coadd image for 

#### Useful information about SDSS images: 

* For more detailed information, please read: 
    - [SDSS-III Terminology](https://www.sdss3.org/dr10/help/glossary.php)
    - [Data Model for SDSS Corrected Frame File](http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html) 

* SDSS is a "drift-scan" survey, which generate 6 long fits image in a single
    exposure (called "run"); Later, the fits image from each of the camera
    column (called "camcol") is separated into a series of smaller images
    (called "field"; each has a dimension of 2048x1489 pixels); The same data
    have been reduced by different version of pipeline, each time they are run
    through a different pipeline, a unique ID is provided (called "RERUN") 
* So, a specific combination of RUN-RERUN-CAMCOL-FIELD corresponds to a specific
    region of sky 
    - For the DR10 archive, all the images have the same RERUN: 301 
    - For each field, all the detected objects have been given a unique index,
       "OBJ"; Hence, RUN-RERUN-CAMCOL-FIELD-OBJ can be used to find all the
       detected object on SDSS images 
    - With a "FILTER" identifier (from "u", "g", "r", "i", and "z"), any SDSS
       image file can be selected (as file, it is usually called a "FRAME")  
    - For now (and probably forever...), SDSS has 938046 fields 

* Once you have the RUN-RERUN-CAMCOL-FIELD-FILTER combination you want, you can
    simply download the FITS file for the corrected, calibrated images from: 
    - http://data.sdss3.org/sas/dr10/boss/photoObj/frames/RERUN/RUN/CAMCOL/frame-FILTER-RUN6-CAMCOL-FIELD.fits.bz2
    - Where RUN6 is the zero-padded version of RUN: e.g. RUN=1234, RUN6=001234

* For each of the field, the SDSS pipeline measures a lot of useful information,
    and keep them in a catalog called "Field" in the on-line database. 
    - The description of this catalog can be found [here](http://skyserver.sdss3.org/public/en/help/browser/browser.aspx#&&history=description+Field+U)
    - It can be accessed through on-line SQL searched from Casjob
    - Another easy way to do a quick search is to use the "sqlcl.py" Python
        script, which you can find
        [here](http://skyserver.sdss3.org/dr10/en/help/download/sqlcl/sqlcl.aspx)
    - However, to use that, certain amount of knowledge about SDSS database, and
        the SQL language is still required

* The information we are interested in this catalog includes: 
    1. FIELDID, RUN, RERUN, CAMCOL, FIELD 
    2. nTotal, nObjects, nChild, nGalaxy, nStars 
    3. nBrightObj_[FILTER] 
    4. Quality: 1:Bad; 2:Acceptable; 3:Good; 4:Missing; 5:Hole
    5. mjd_[FILTER] 
    6. airmass_[FILTER] 
    7. ra, dec: RA, DEC of the center 
    8. raMin, raMax, decMin, decMax 
    9. sky_[FILTER]: in unit of nmgy/arcsec^2
    10. psfNStar_[FILTER] 
    11. psfWidth_[FILTER]: in unit of arcsec
    12. gain_[FILTER]: in unit of electrons/DN 
    13. darkVariance_[FILTER] 
    14. score: between [0-1], also quality of the field 
    15. aterm_[FILTER]; kterm_[FILTER]; kdot_[FILTER]; Parameter for photometry
          calibration 
    16. calibStatus_[FILTER]: "Photometric", "Unphot_overlap", 
        "Unphot_Extrap_Clear", "Unphot_Extrap_Cloudy", "Unphot_Disjoint"
        - Use "fCalibStatusN()" function to get the strings 
    17. imageStatus_[FILTER]: "Bad_Astrom", "Bad_Focus", "Bad_Rotator", "Clear" 
        "Cloudy", "Dead_CCD", "Noisy_CCD", "FF_Petals", "Shutters" 
        - Use "fImageStatusN()" function to get the strings
    18. nMgyPerCount_[FILTER]: in unit of nmgy/count
    19. primaryArea: in unit of deg^2
        - if primaryArea == 0, then this field is not primary anywhere
    20. saturationLevel_[FILTER]: in counts

* Some basic information for each field has been stored in the
    "sdss_field_sum.fit" file; it a simple FITS ASCII catalog




