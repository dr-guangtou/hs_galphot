;+
; NAME:
;   ml_mgframeread
;
; PURPOSE:
;   Read data from MaNGA data files; either mgFrame, mgSFrame, mgFFrame, or
;   mgCFrame files.  This code started life as spframe_read.pro and
;   has been adapted to the MaNGA data format, in addition to some
;   other tweaks.
;
; CALLING SEQUENCE:
;   ml_mgframeread, filename, [ indx, objflux=, objivar=, mask=, $
;    wset=, loglam=, dispset=, dispimg=, ximg=, $
;    slitmap=, slithdr=, skyflux=, superflat=, hdr=, adderr= ]
;
; INPUTS:
;   filename   - Input file name
;
; OPTIONAL INPUTS:
;   indx       - Optional 0-indexed row numbers; default to all
;   adderr     - Additional error to add to the formal errors, as a
;                fraction of the flux; default to none
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   objflux    - Object flux
;   objivar    - Object inverse variance (units of 1/OBJFLUX^2)
;   mask       - Pixel bit mask
;   wset       - Trace-set for wavelength solution
;   loglam     - Wavelength image (vacuum log-10 Ang)
;   dispset    - Trace-set for dispersion solution
;   dispimg    - Dispersion image (per native pixel)
;   slitmap    - Slitmap for the exposure
;   slithdr    - Slitmap header
;   ximg       - X position on CCD image
;   skyflux    - Sky flux (same units as OBJFLUX)
;   superflat  - Superflat vector from quartz lamps
;   hdr        - FITS header for HDU#0
;
; COMMENTS:
;   The mgFrame, mgSFrame, mgFFrame, and mgCFrame files contain similar HDUs,
;   except for the wavelength and dispersion extensions which are
;   trace-sets in mgFrame, mgSFrame, and mgFFrame and more easily-interpreted
;   2-dimensional images in mgCFrame.  All HDUs are called by their EXTNAME
;   rather than their HDU number.
; 
;   Note that the mgFrame, mgSFrame, and mgFFrame files are per-ccd, so they have
;   about 700 rows and native wavelength columns.  In contrast,
;   mgCFrame combines cameras and channels (b1/b2,r1/r2) into a single
;   file, so it has 1423 rows and common wavelength columns.
;
;   The slitmap FITS extension header contains the yanny-format slitmap
;   header within it.  We parse out the extra FITS header stuff and
;   return just the original yanny-format FITS header.  This relies on the
;   convention that slitmap header cards start with lower case letters.
;
;   mgFrame: Extracted science frame
;   HDU #0 []:          Empty except for global header
;   HDU #1 [FLUX]:      Flux in flatfielded e- [FLOAT, CCDROW x NFIBER]
;   HDU #2 [IVAR]:      Inverse variance [FLOAT, CCDROW x NFIBER]
;   HDU #3 [MASK]:      Pixel mask [LONG, CCDROW x NFIBER]
;   HDU #4 [WSET]:      Wavelength solution traceset in log10 Angstroms (vacuum) [BINARY FITS TABLE]
;   HDU #5 [DISPSET]:   Wavelength dispersion in units of (10^-4 log10 wavelength) [BINARY FITS TABLE]
;   HDU #6 [SLITMAP]:   Slitmap [BINARY FITS TABLE]
;   HDU #7 [XPOS]:      Xpositions of traces on CCD [FLOAT, CCDROW x NFIBER]
;   HDU #8 [SUPERFLAT]: Superflat vector from quartz lamps [FLOAT, CCDROW x NFIBER]
;
;   mgSFrame: Sky-subtracted science frame
;   HDU #0 []:          Empty except for global header
;   HDU #1 [FLUX]:      Sky-subtracted flux in flatfielded e- [FLOAT, CCDROW x NFIBER]
;   HDU #2 [IVAR]:      Inverse variance [FLOAT, CCDROW x NFIBER]
;   HDU #3 [MASK]:      Pixel mask [LONG, CCDROW x NFIBER]
;   HDU #4 [WSET]:      Wavelength solution traceset in log10 Angstroms (vacuum) [BINARY FITS TABLE]
;   HDU #5 [DISPSET]:   Wavelength dispersion in units of (10^-4 log10 wavelength) [BINARY FITS TABLE]
;   HDU #6 [SLITMAP]:   Slitmap [BINARY FITS TABLE]
;   HDU #7 [XPOS]:      Xpositions of traces on CCD [FLOAT, CCDROW x NFIBER]
;   HDU #8 [SUPERFLAT]: Superflat vector from quartz lamps [FLOAT, CCDROW x NFIBER]
;   HDU #9 [SKY]:       Sky flux in flatfielded e- [FLOAT, CCDROW x NFIBER]
;
;   mgFFrame: Flux-calibrated science frame
;   HDU #0 []:          Empty except for global header
;   HDU #1 [FLUX]:      Sky-subtracted flux in 1e-17 erg/s/cm2/Ang [FLOAT, CCDROW x NFIBER]
;   HDU #2 [IVAR]:      Inverse variance [FLOAT, CCDROW x NFIBER]
;   HDU #3 [MASK]:      Pixel mask [LONG, CCDROW x NFIBER]
;   HDU #4 [WSET]:      Wavelength solution traceset in log10 Angstroms (vacuum) [BINARY FITS TABLE]
;   HDU #5 [DISPSET]:   Wavelength dispersion in units of (10^-4 log10 wavelength) [BINARY FITS TABLE]
;   HDU #6 [SLITMAP]:   Slitmap [BINARY FITS TABLE]
;   HDU #7 [XPOS]:      Xpositions of traces on CCD [FLOAT, CCDROW x NFIBER]
;   HDU #8 [SUPERFLAT]: Superflat vector from quartz lamps [FLOAT, CCDROW x NFIBER]
;   HDU #9 [SKY]:       Sky flux in 1e-17 erg/s/cm2/Ang [FLOAT, CCDROW x NFIBER]
;
;   mgCFrame: Flux Calibrated and Camera Combined frame on a Common wave grid
;   HDU #0 []:          Empty except for global header
;   HDU #1 [FLUX]:      Sky-subtracted flux in 1e-17 erg/s/cm2/Ang [FLOAT, CCDROW x NFIBER]
;   HDU #2 [IVAR]:      Inverse variance [FLOAT, CCDROW x NFIBER]
;   HDU #3 [MASK]:      Pixel mask [LONG, CCDROW x NFIBER]
;   HDU #4 [WAVE]:      Wavelength solution in log10(Ang) [FLOAT, NWAVE]
;   HDU #5 [DISP]:      Wavelength dispersion in log10(Ang)[FLOAT, NWAVE x NFIBER?]
;   HDU #6 [SLITMAP]:   Slitmap [BINARY FITS TABLE]
;   HDU #7 [SKY]:       Sky flux in 1e-17 erg/s/cm2/Ang [FLOAT, NWAVE x NFIBER]
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   fileandpath()
;   headfits()
;   lookforgzip()
;   mrdfits()
;   traceset2xy
;   traceset_trim()
;
; REVISION HISTORY:
;   05-Feb-2004  Written by D. Schlegel, Princeton
;   04-Sep-2013  Adapted for MaNGA by D. Law, Dunlap Institute
;   07-May-2014  Added mgFFrame support (D. Law)
;-
;------------------------------------------------------------------------------
pro ml_mgframeread, filename, indx, objflux=objflux, objivar=objivar, $
 mask=mask, wset=wset, loglam=loglam, dispset=dispset, dispimg=dispimg, $
 ximg=ximg, slitmap=slitmap, slithdr=slithdr, skyflux=skyflux, superflat=superflat, $
 hdr=hdr, adderr=adderr

  ; Was it requested to trim the input data to select rows?
  qtrim = n_elements(indx) GT 0

  ; If no filename specified, error.
  if (~keyword_set(filename)) then begin
    splog,'ERROR: Must specify FILENAME'
    return
  endif

  ; Check to see if the desired file is gzipped, and change name if so
  thisfile = lookforgzip(filename[0])

  ; If file not found, error
  if (~file_test(thisfile)) then begin
    splog,'ERROR: FILENAME not found'
    return
  endif

  ; Is this an extracted (mgFrame) file?
  qeframe = strmatch(fileandpath(filename),'mgFrame*')
  ; Is this a sky-subtracted (mgSFrame) file?
  qsframe = strmatch(fileandpath(filename),'mgSFrame*')
  ; Is this a flux-calibrated (mgFFrame) file?
  qfframe = strmatch(fileandpath(filename),'mgFFrame*')
  ; Is this a camera-combined (mgCFrame) file?
  qcframe = strmatch(fileandpath(filename),'mgCFrame*')

  ; If filename is neither of the above, error
  if ((~qeframe)and(~qsframe)and(~qfframe)and(~qcframe)) then begin
    splog, 'ERROR: Invalid FILENAME, must specify mgFrame, mgSFrame, mgFFrame, or mgCFrame.'
    return
  endif

  ; If a mgFrame file was chosen, check that requested extensions are valid
  if ((qeframe)and(arg_present(skyflux))) then begin
    splog,'ERROR: Cannot request SKY extension from an mgFrame file.'
  endif

  ; If a mgSFrame file was chosen, check that requested extensions are valid
  ; (All currently valid)

  ; If a mgFFrame file was chosen, check that requested extensions are valid
  ; (All currently valid)

  ; If a mgCFrame file was chosen, check that requested extensions are valid
  if ((qcframe)and((arg_present(ximg))or(arg_present(superflat)) $
   or(arg_present(wset))or(arg_present(dispset)))) then begin
     splog,'ERROR: Cannot request XIMG, SUPERFLAT, DISPSET, or WSET extensions from an mgCFrame file.'
  endif

  ; Get the image header (whatever is in extension 0)
  if (arg_present(hdr)) then hdr = headfits(thisfile[0])

  ; Get the flux extension.  Also read it into this function if both objivar
  ; and adderr keywords are set, as it is required for them.
  if (arg_present(objflux) $
   OR (arg_present(objivar) AND keyword_set(adderr))) then begin
     objflux = mrdfits(thisfile[0], 'FLUX', /silent)
     ; If desired, trim the rows returned
     if (qtrim) then objflux = objflux[*,indx]
  endif

  ; Get the inverse variance extension
  if (arg_present(objivar)) then begin
    objivar = mrdfits(thisfile[0], 'IVAR', /silent)
    ; If desired, trim the rows returned
    if (qtrim) then objivar = objivar[*,indx]
    ; If desired, add in error proportional to the flux
    if (keyword_set(adderr)) then begin
      gmask = objivar NE 0 ; =1 for good points
      objivar = 1.0 / ( 1.0/(objivar + (1-gmask)) $
       + (adderr * (objflux>0))^2 ) * gmask
    endif
  endif

  ; Get the pixel mask
  if (arg_present(mask)) then begin
    mask = mrdfits(thisfile[0], 'MASK', /silent)
    ; If desired, trim the rows returned
    if (qtrim) then mask = mask[*,indx]
  endif

  ; Get the wavelength solution if reading from mgFrame, mgSFrame, mgFFrame
  if (((qeframe)or(qsframe)or(qfframe))and((arg_present(wset))or(arg_present(loglam)))) then begin
    wset = mrdfits(thisfile[0], 'WSET', /silent)
    ; If desired, trim the rows returned
    if (qtrim) then wset = traceset_trim(wset, indx)
    ; Compute the wavelength image if desired
    if (arg_present(loglam)) then begin
      traceset2xy, wset, xtmp, loglam
      ; Reset the xtmp dummy variable
      xtmp=0
    endif
  endif

  ; Get the wavelength image if reading from mgCFrame
  if ((qcframe)and(arg_present(loglam))) then begin
    loglam = mrdfits(thisfile[0], 'WAVE', /silent)
    ; If desired, trim the rows returned
    if (qtrim) then loglam = loglam[*,indx]
  endif

  ; Get the dispersion solution if reading from mgFrame, mgSFrame, mgFFrame
  if (((qeframe)or(qsframe)or(qfframe))and((arg_present(dispset))or(arg_present(dispimg)))) then begin
    dispset = mrdfits(thisfile[0], 'DISPSET', /silent)
    ; If desired, trim the rows returned
    if (qtrim) then dispset = traceset_trim(dispset, indx)
    ; Compute the dispersion image if desired
    if (arg_present(dispimg)) then begin
      traceset2xy, dispset, xtmp, dispimg
      ; Reset the xtmp dummy variable
      xtmp=0
    endif
  endif

  ; Get the dispersion image if reading from mgCFrame
  if ((qcframe)and(arg_present(dispimg))) then begin
    dispimg = mrdfits(thisfile[0], 'DISP', /silent)
    ; If desired, trim the rows returned
    if (qtrim) then dispimg = dispimg[*,indx]
  endif

  ; Get the slitmap
  if ((arg_present(slitmap))or(arg_present(slithdr))) then begin
    slitmap = mrdfits(thisfile[0], 'SLITMAP', slithdr, /silent)
    ; If desired, trim the rows returned
    if (qtrim) then slitmap = slitmap[indx]
    ; Trim the slithdr to eliminate FITS addenda
    ; (lines that start with upper case, not pound sign)
    ml_cleanslithdr,slithdr
  endif

  ; Get the sky extension from mgSFrame, mgFFrame, or mgCFrame
  if (((qsframe)or(qcframe)or(qfframe))and(arg_present(skyflux))) then begin
    skyflux = mrdfits(thisfile[0], 'SKY', /silent)
    ; If desired, trim the rows returned
    if (qtrim) then skyflux = skyflux[*,indx]
  endif

  ; Get the x position image from mgFrame, mgSFrame, mgFFrame
  if (((qeframe)or(qsframe)or(qfframe))and(arg_present(ximg))) then begin
    ximg = mrdfits(thisfile[0], 'XPOS', /silent)
    ; If desired, trim the rows returned
    if (qtrim) then ximg = ximg[*,indx]
  endif

  ; Get the superflat from mgFrame, mgSFrame, mgFFrame
  if (((qeframe)or(qsframe)or(qfframe))and(arg_present(superflat))) then begin
    superflat = mrdfits(thisfile[0], 'SUPERFLAT', /silent)
    ; If desired, trim the rows returned
    if (qtrim) then superflat = superflat[*,indx]
  endif

return
end
;------------------------------------------------------------------------------
