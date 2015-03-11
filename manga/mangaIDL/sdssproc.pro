;+
; NAME:
;   sdssproc
;
; PURPOSE:
;   Read in raw SDSS files, and process with opConfig, opECalib, opBC par files.
;
; CALLING SEQUENCE:
;   sdssproc, infile, [image, invvar, indir=, $
;    outfile=, nsatrow=, fbadpix=, $
;    hdr=hdr, configfile=, ecalibfile=, bcfile=, $
;    /applybias, /applypixflat, /silent, /do_lock, minflat=, maxflat=, $
;    spectrographid=, color=, camname=, /applycrosstalk, ccdmask= ]
;
; INPUTS:
;   infile     - Raw SDSS file name
;
; OPTIONAL KEYWORDS:
;   indir      - Input directory for INFILE
;   outfile    - Calibrated 2D frame, HDU#0 is the image and #1 is invvar;
;                if set as /OUTFILE, then default name is constructed from
;                INFILE as sdProc-XX-XXXXXXXX.fits, excluding any path info
;   nsatrow    - Number of saturated rows, assuming that a row is saturated
;                if at least 20 of its pixels are above saturation level
;   fbadpix    - Fraction of bad pixels, not including bad columns
;   configfile - Default to "opConfig*par", selecting the file with the
;                appropriate MJD.
;   ecalibfile - Default to "opECalib*par", selecting the file with the
;                appropriate MJD.
;   bcfile     - Default to "opBC*par", selecting the file with the
;                appropriate MJD.
;   applybias  - Apply 2-D bias image.
;   applypixflat- Apply 2-D pixel-to-pixel flat (after subtracting bias).
;   silent     - If set, then don't output any text.
;   do_lock    - If set, then lock the "sdHdrFix-$MJD.par" file
;                using DJS_LOCKFILE().
;   minflat    - Minimum values allowed for pixflat; pixels with the
;                flat out of range are masked; default to 0.
;   maxflat    - Maximum values allowed for pixflat; pixels with the
;                flat out of range are masked; default to 1e10.
;   applycrosstalk - choose whether or not to apply crosstalk correction
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   image      - Processed 2d image
;   invvar     - Associated inverse variance
;   hdr        - Processed FITS header
;   spectrographid - Return spectrograph ID (1 or 2)
;   color      - Return spectrograph color ('red' or 'blue')
;   camname    - Return camera name: 'b1', 'r1', 'b2', or 'r2'
;   ccdmask    - Image with 'NODATA' bit set for parts of BOSS data
;
; COMMENTS:
;   Only the header is read from the image if IMAGE, INVVAR, OUTFILE and
;   VARFILE are all not set.
;
;   Required header keywords: EXPTIME.
;
;   The returned image is in electrons, not ADU.
;
;   The signal-to-noise is limited to never exceed 100, by adding 1.e-4
;   times the flux squared to the variance term.
;
;   Change the CAMERAS keyword to the camera as specified by the file name.
;
;   Rename 'target' to 'science', and 'calibration' to 'arc' in the
;   header keyword FLAVOR.
;
;   Determine the exposure number from the file name itself.
;
; BUGS:
;   The open-shutter correction SMEARIMG will include smeared data from
;   any cosmic rays, which is wrong.  At the minimum, I could interpolate
;   over A/D saturations (in ADMASK) before constructing SMEARIMG.
;
; PROCEDURES CALLED:
;   djs_filepath()
;   djs_iterstat
;   fileandpath()
;   findopfile()
;   fits_purge_nans
;   fits_wait()
;   headfits()
;   idlspec2d_version()
;   idlutils_version()
;   lookforgzip()
;   ml_mghdrfix
;   mwrfits
;   rdss_fits()
;   splog
;   sxaddpar
;   sxpar()
;   yanny_free
;   yanny_read
;
; INTERNAL SUPPORT ROUTINES:
;   make_badcolumn_mask()
;
; DATA FILES:
;   $MANGADRP_DIR/examples/opConfig*par
;   $MANGADRP_DIR/examples/opECalib*par
;   $MANGADRP_DIR/examples/opBC*par
;   $SPECFLAT_DIR/biases/pixbias*.fits*
;   $SPECFLAT_DIR/flats/pixflat*.fits*
;
; REVISION HISTORY:
;   13-May-1999  Written by Scott Burles & David Schlegel, Apache Point.
;   08-Sep-1999  Modified to read Yanny param files instead of FITS
;                versions of the same (DJS).
;   01-Dec-1999  Added version stamping (DJS).
;   07-Dec-1999  Mask neighbors of pixels that saturated the A/D converter.
;                Identify blead trails and mask from that row up (DJS).
;   10-Dec-1999  Test if the shutter was open during readout, and try
;                to correct the light for that (DJS).
;   04-Feb-2000  Declare that the shutter was open if it is a >640 sec
;                exposure taken before MJD=51570 (DJS).
;   26-Jul-2000  Added fix for "dropped pixel" problem for data on or after
;                MJD=51688 (23 May 2000).  Should disable this code for later
;                MJD's once this problem is fixed in the electronics.
;   26-Jul-2000  Added fix for more severe "shifted row" electronics problem
;                for data taken on MJD=51578 to 51580 (3 to 5 Feb 2000).
;   26-Aug-2000  Horrible kludge for this night MJD=51779 (22/23 Aug 2000).
;                Add a noise term of 100 ADU to the left amplifier of r2.
;                That amplifier had random bits being set incorrectly,
;                in particular the 32 bit and 256 bit.
;   04-Nov-2000  Measure the bias values using DJS_ITERSTAT instead of MEDIAN,
;                since the median is always only an integer value.
;   31-Jan-2001  Determine the exposure number from the file name itself,
;                since the counting got off by one on MJD=51882.
;   08-Aug-2011: Changing bias-subtraction recipe for BOSS (A. Bolton, Utah)
;      Jul-2013  Stripped Out SDSSI stuff for MaNGA (B. Cherinka, Toronto)
;-
;------------------------------------------------------------------------------
pro sdssproc_badformat, image, camname=camname, mjd=mjd

  on_error,0
  compile_opt idl2
  compile_opt idl2, hidden
  
   dims = size(image,/dimens)
   nx = dims[0]
   ny = dims[1]

   case camname of
   'b1': begin
       xs = [77,4274,77,4274]
       xw = [1,1,1,1]
       thresh = [15,15,0,0]
       end
   'b2': begin
       xs = [77,4274,77,4274]
       xw = [1,1,1,1]
       thresh = [5,10,10,0]
       end
   'r1': begin
       ; The transients on the right of r1 are really 2 columns wide
       xs = [111,4240,111,4240]
       xw = [1,2,1,2]
       thresh = [30,20,60,15]
       end
   'r2': begin
          ;mjd gt 55300
          xs = [112,4239,111,4240]
          xw = [1,1,1,1]
          thresh = [5,10,12,5]
       end
   endcase

   nsub = 7
   reg0 = image[xs[0]:xs[0]+nsub-1,0:ny/2-1]
   reg1 = rotate(image[xs[1]-nsub+1:xs[1],0:ny/2-1],5)
   reg2 = rotate(image[xs[2]:xs[2]+nsub-1,ny/2:ny-1],7)
   reg3 = rotate(image[xs[3]-nsub+1:xs[3],ny/2:ny-1],2)

   if (xw[0] GT 1) then reg0 = convol(reg0, [intarr(xw[0]-1),1,intarr(xw[0]-1)+1], /edge_wrap)
   if (xw[1] GT 1) then reg1 = convol(reg1, [intarr(xw[1]-1),1,intarr(xw[1]-1)+1], /edge_wrap)
   if (xw[2] GT 1) then reg2 = convol(reg2, [intarr(xw[2]-1),1,intarr(xw[2]-1)+1], /edge_wrap)
   if (xw[3] GT 1) then reg3 = convol(reg3, [intarr(xw[3]-1),1,intarr(xw[3]-1)+1], /edge_wrap)

   shiftvec = lonarr(ny/2) ; Pixel shift of each row (same in all quadrants)
   nshift = lonarr(nsub) ; Keep track of # of rows shifted each pixel amount

   peak = intarr(4)
   for iy=0L, ny/2-1L do begin
      ; Find the peak value relative to where it should be in this row
      max0 = max(reg0[0:nsub-1-xw[0],iy], peak0)
      max1 = max(reg1[0:nsub-1-xw[1],iy], peak1)
      max2 = max(reg2[0:nsub-1-xw[2],iy], peak2)
      max3 = max(reg3[0:nsub-1-xw[3],iy], peak3)
      ; Test that the peak is at least "thresh" larger than the next pix, where that threshold is different in each amplifier
      peak0 *= (reg0[peak0,iy] GE reg0[peak0+xw[0],iy] + thresh[0])
      peak1 *= (reg1[peak1,iy] GE reg1[peak1+xw[1],iy] + thresh[1])
      peak2 *= (reg2[peak2,iy] GE reg2[peak2+xw[2],iy] + thresh[2])
      peak3 *= (reg3[peak3,iy] GE reg3[peak3+xw[3],iy] + thresh[3])
      peak = [peak0,peak1,peak2,peak3]
      bpeak = peak[(sort(peak))[1]]
      ; If at least 3 quadrants appear shifted the same amount, then use that shift
      if (total(peak EQ bpeak) GE 3) then begin
         shiftvec[iy] = bpeak
         nshift[bpeak]++
      endif else begin
         ; If no shifts detected, then assume it is the same as the previous row
         if (iy GT 0) then shiftvec[iy] = shiftvec[iy-1]
      endelse
   endfor

   ; Do the actual image shifts
   for iy=0L, ny/2-1L do begin
      if (shiftvec[iy] GT 0) then begin
         image[0:nx/2-1,iy] = shift(image[0:nx/2-1,iy], -shiftvec[iy])
         image[0:nx/2-1,ny-1-iy] = shift(image[0:nx/2-1,ny-1-iy], -shiftvec[iy])
         image[nx/2:nx-1,iy] = shift(image[nx/2:nx-1,iy], shiftvec[iy])
         image[nx/2:nx-1,ny-1-iy] = shift(image[nx/2:nx-1,ny-1-iy], shiftvec[iy])
      endif
   endfor

   nunknown = ny - long(2*total(nshift))
   if (nunknown GT 0) then splog, 'Warning: Serial transient not detected in ', nunknown, ' rows'

   for i=1, nsub-1 do $
    if (nshift[i] GT 0) then splog, 'WARNING: Electronics shifted ', 2*nshift[i], ' rows by ', i, ' pix'

   return
end
;------------------------------------------------------------------------------
;  Create the bad column mask (1 for a masked pixel) with image size nc,nr
;  If the operation doesn't work just return 0 for no masked pixels

function make_badcolumn_mask, bcfile, camrow, camcol, nc=nc, nr=nr, silent=silent

  on_error,0
  compile_opt idl2
  compile_opt idl2, hidden
      
  if ~keyword_set(nc) then nc=2048L
  if ~keyword_set(nr) then nr=2048L
  
  yanny_read, bcfile, pdata
  if (size(pdata,/tname)) EQ 'INT' then begin
    if (~keyword_set(silent)) then splog, 'WARNING: Could not read BC file ' + fileandpath(bcfile)
    return, 0
  endif
  
  bc = *pdata[0]
  yanny_free, pdata
  
  ibc = where(bc.camrow EQ camrow AND bc.camcol EQ camcol, nbc)
  if (~keyword_set(nbc)) then begin
    if (~keyword_set(silent)) then splog,'WARNING: Could not find this camera info in BC file ' + fileandpath(bcfile)
    return, 0
  endif
  
  bc = bc[ibc]
  bcmask = bytarr(nc, nr)
  
  ; Mask out bad columns
  
  if (nbc GT 0) then begin
    bcsc = (bc.dfcol0 > 0) < (nc-1)
    bcec = (bc.dfcol0 + bc.dfncol - 1 < (nc-1)) > bcsc
    bcsr = (bc.dfrow0 > 0) < (nr-1)
    bcer = (bc.dfrow0 + bc.dfnrow - 1 < (nr-1)) > bcsr
    
    for i=0, nbc-1 do bcmask[bcsc[i]:bcec[i],bcsr[i]:bcer[i]] = 1
  endif
  
  return, bcmask
end

;------------------------------------------------------------------------------
pro sdssproc, infile1, image, invvar, indir=indir, $
 outfile=outfile1, nsatrow=nsatrow, fbadpix=fbadpix, $
 hdr=hdr, configfile=configfile, ecalibfile=ecalibfile, bcfile=bcfile, $
 applybias=applybias, applypixflat=applypixflat, silent=silent, $
 do_lock=do_lock, minflat=minflat, maxflat=maxflat, $
 spectrographid=spectrographid, color=color, camname=camname, $
 applycrosstalk=applycrosstalk, ccdmask=ccdmask, VISUAL=VISUAL, survey=survey

  on_error,0
  compile_opt idl2

  common com_sdssproc, vers2d, versutils, versflat, verslog
    
  if (N_params() LT 1) then begin
    doc_library, 'sdssproc'
    return
  endif

  infile = infile1[0]
  readimg = arg_present(image) OR keyword_set(outfile1)
  readivar = arg_present(invvar) OR keyword_set(outfile1) OR arg_present(nsatrow) OR arg_present(fbadpix)

  fullname = djs_filepath(infile, root_dir=indir)
  fullname = (lookforgzip(fullname, count=ct))[0]
  if (ct NE 1) then message, 'Cannot find image ' + infile
    
  if (keyword_set(outfile1)) then begin
     if (size(outfile1,/tname) EQ 'STRING') then begin
        outfile = outfile1
     endif else begin
        outfile = 'sdProc-' + strmid(fileandpath(infile),4,11)+'.fits'
     endelse
  endif

  if (readimg OR readivar) then rawdata = rdss_fits(fullname, hdr, /nofloat, silent=silent) $
    else hdr = headfits(fullname)

  sxdelpar, hdr, 'BZERO'
  sxdelpar, hdr, 'BSCALE'
  sxdelpar, hdr, 'CHECKSUM'
  sxdelpar, hdr, 'DATASUM'

  ;-----------
  ; Determine the exposure number from the file name itself.
  ; Very bad form, but this information is sometimes wrong in the header.
  ; In particular, the counting got off by 1 on MJD=51882.
  
  i = strpos(infile, '-', /reverse_search)
  expnum = long( strmid(infile, i+1, 8) )
  if (~keyword_set(expnum)) then message, 'Cannot determine exposure number from file name ' + infile

  ;-----------
  ; Fix the headers with any hand-edits that we have determined.

  ml_mghdrfix, infile, hdr, silent=silent
  
  ;-----------
  ; Replace exposure number with that found in the file name.
  hdrexp = sxpar(hdr, 'EXPOSURE')
  if (expnum NE hdrexp) then begin
    if (~keyword_set(silent)) then $
      splog, 'WARNING: Exposure number in header (' + strtrim(string(hdrexp),2) + ') disagrees w/filename (' + strtrim(string(expnum),2) + ') !!'
    sxaddpar, hdr, 'EXPOSURE', expnum
  endif
  
  ;-----------
  ; Determine which CCD from the file name itself, using either the numbering scheme (01,02,03,04) or naming scheme (b1,r2,b2,r1).
  ; Very bad form, but this information is not in the header since the CAMERAS keyword is sometimes wrong.
  i = strpos(infile, '-', /reverse_search)
  if (i[0] EQ -1 OR i-2 LT 0) then message, 'Cannot determine CCD number from file name ' + infile
    
  camnames = ['b1', 'r2', 'b2', 'r1']
  camnums = ['01', '02', '03', '04']
  
  ; First try to match a camera name (e.g., 'b1'), then try to match a camera number (e.g., '01').  If both fail, then abort.
  indx = where(strmid(infile, i-2, 2) EQ camnames, ct)
  if (ct NE 1) then indx = where(strmid(infile, i-2, 2) EQ camnums, ct)
  if (ct NE 1) then message, 'Cannot determine CCD number from file name ' + infile
    
  ; Do not read the camera from the CAMERAS keyword, since this was often wrong in the early days!
  camname = camnames[indx[0]]
  case camname of
    'b1': begin
      spectrographid = 1
      color = 'blue'
    end
    'r1': begin
      spectrographid = 1
      color = 'red'
    end
    'b2': begin
      spectrographid = 2
      color = 'blue'
    end
    'r2': begin
      spectrographid = 2
      color = 'red'
    end
  endcase
  camcol = indx[0] + 1
  camrow = 0

  ; Cache these versions in a common block, since the calls to spawn are slow
  if (~keyword_set(vers2d)) then begin
     vers2d = mangadrp_version(/simple)
     versutils = idlutils_version()
     spawn, 'specflat_version', versflat, err, /noshell
     if (~keyword_set(versflat)) then versflat = 'Unknown'
  endif
  
  sxaddpar, hdr, 'CAMROW', camrow
  sxaddpar, hdr, 'CAMCOL', camcol
  sxaddpar, hdr, 'TELESCOP', 'SDSS 2.5-M', ' Sloan Digital Sky Survey'
  hdr=mdrp_makeheader(head=hdr,/drp2d)
    
  ;-----------
  ; Rename 'target' -> 'science', and 'calibration' -> 'arc'
    
  mjd = sxpar(hdr, 'MJD')
  flavor = strtrim(sxpar(hdr, 'FLAVOR'),2)
  if (flavor EQ 'target') then flavor = 'science'
  if (flavor EQ 'calibration') then flavor = 'arc'
  if (sxpar(hdr, 'COLBIN') GT 1 OR sxpar(hdr, 'ROWBIN') GT 1) then flavor = 'unknown'
  sxaddpar, hdr, 'FLAVOR', flavor
  sxaddpar, hdr, 'CAMERAS', camname

  ;--------------
  ; Flag to trigger bolton bias subtraction for survey-quality BOSS data:
  bossgood = 1B 

  ; Mark when the BOSS red CCDs switched from 1-phase to 2-phase readout
  if (mjd GE 55415 AND strmatch(camname,'r*')) then sxaddpar, hdr, 'TWOPHASE', 'T' else sxaddpar, hdr, 'TWOPHASE', 'F'
  
  ;-----------
  ; Deal with the first light BOSS data on MJD 55052

if (mjd GE 55052) then begin
   ; Declare any BOSS exposures 'excellent' if that keyword missing, unless it was a Hartmann exposure
   junk = sxpar(hdr, 'QUALITY', count=ct)
   if (ct EQ 0) then begin
      hartmann = strtrim(sxpar(hdr, 'HARTMANN'),2)
      if (hartmann EQ 'Left' OR hartmann EQ 'Right') then sxaddpar, hdr, 'QUALITY', 'test' else sxaddpar, hdr, 'QUALITY', 'excellent'
   endif
   if (readimg OR readivar) then begin
      rawdata += 32768.
      sdssproc_badformat, rawdata, camname=camname, mjd=mjd

      case strmid(camname,0,1) of
      'b': begin
         case spectrographid of
            1: gain = [1.048, 1.048, 1.018, 1.006] ; b1 gain
            2: gain = [1.040, 0.994, 1.002, 1.010] ; b2 gain
         end
         ; Do bolton bias subtraction for survey-quality BOSS dates:
         ; (Note that these lines are identical between b and r cams.)
         if bossgood then begin
            if (keyword_set(silent) EQ 0) then splog, 'BOSS survey-quality MJD: applying pixbias whether you like it or not!'
            pp = filepath('', root_dir=getenv('SPECFLAT_DIR'), subdirectory='biases')
            pixbiasname = findopfile('boss_pixbias-*-'+camname+'.fits*', mjd, pp, silent=silent)
            image = bolton_biassub(rawdata, pp+pixbiasname, rnoise=rnoise, cam=camname, sigthresh=3.0)
            xwid = (size(image))[1]
            ywid = (size(image))[2]
            rdnoise = rnoise[*] * 1.015 * gain ; account for sigma-clipping
            if (keyword_set(silent) EQ 0) then splog, 'Read noise = ', rdnoise, ' electrons'
            image[0:xwid/2-1,0:ywid/2-1] *= gain[0]
            image[xwid/2:xwid-1,0:ywid/2-1] *= gain[1]
            image[0:xwid/2-1,ywid/2:ywid-1] *= gain[2]
            image[xwid/2:xwid-1,ywid/2:ywid-1] *= gain[3]
         endif 
         
         if (readivar) then begin
            invvar = 0.*image
            invvar[0:2047,0:2055] = 1. / (rdnoise[0]^2 + (image[0:2047,0:2055]>0))
            invvar[0:2047,2056:4111] = 1. / (rdnoise[1]^2 + (image[0:2047,2056:4111]>0))
            invvar[2048:4095,0:2055] = 1. / (rdnoise[2]^2 + (image[2048:4095,0:2055]>0))
            invvar[2048:4095,2056:4111] = 1. / (rdnoise[3]^2 + (image[2048:4095,2056:4111]>0))
            mask = 0.*image
            mask[*,696:3516] = 1
            invvar *= mask
         endif
         end
      'r': begin
         case spectrographid of
            1: gain = [1.9253, 1.5122, 1.4738, 1.5053]  ; R1 replaced summer 2011
            2: gain = [1.598, 1.656, 1.582, 1.594] ; R2 replaced April 2011
         end
          ; Do bolton bias subtraction for survey-quality BOSS dates:
         ; (Note that these lines are identical between b and r cams.)
         if bossgood then begin
            if (keyword_set(silent) EQ 0) then splog, 'BOSS survey-quality MJD: applying pixbias whether you like it or not!'
            pp = filepath('', root_dir=ml_getenv('SPECFLAT_DIR'), subdirectory='biases')
            pixbiasname = findopfile('boss_pixbias-*-'+camname+'.fits*', mjd, pp, silent=silent)
            image = bolton_biassub(rawdata, pp+pixbiasname, rnoise=rnoise, cam=camname, sigthresh=3.0)
            xwid = (size(image))[1]
            ywid = (size(image))[2]
            rdnoise = rnoise[*] * 1.015 * gain ; account for sigma-clipping
            if (keyword_set(silent) EQ 0) then splog, 'Read noise = ', rdnoise, ' electrons'
            image[0:xwid/2-1,0:ywid/2-1] *= gain[0]
            image[xwid/2:xwid-1,0:ywid/2-1] *= gain[1]
            image[0:xwid/2-1,ywid/2:ywid-1] *= gain[2]
            image[xwid/2:xwid-1,ywid/2:ywid-1] *= gain[3]
         endif 
         
         if (readivar) then begin
            invvar = 0.*image
            invvar[0:2056,0:2063] = 1. / (rdnoise[0]^2 + (image[0:2056,0:2063]>0))
            invvar[0:2056,2064:4127] = 1. / (rdnoise[1]^2 + (image[0:2056,2064:4127]>0))
            invvar[2057:4113,0:2063] = 1. / (rdnoise[2]^2 + (image[2057:4113,0:2063]>0))
            invvar[2057:4113,2064:4127] = 1. / (rdnoise[3]^2 + (image[2057:4113,2064:4127]>0))
            mask = 0.*image
            mask[*,28:3668] = 1
            invvar *= mask
         endif
         end
      endcase

      if (arg_present(ccdmask)) then begin
         ccdmask = lonarr(size(image,/dimens))
         ccdmask += (mask EQ 0) * sdss_flagval('SPPIXMASK','NODATA')
      endif

      ; Add read-noise to header
      sxaddpar, hdr, 'RDNOISE0', rdnoise[0], 'CCD read noise amp 0 [electrons]'
      sxaddpar, hdr, 'RDNOISE1', rdnoise[1], 'CCD read noise amp 1 [electrons]'
      sxaddpar, hdr, 'RDNOISE2', rdnoise[2], 'CCD read noise amp 2 [electrons]'
      sxaddpar, hdr, 'RDNOISE3', rdnoise[3], 'CCD read noise amp 3 [electrons]'

      ; Identify CRs, and grow by 1 pix
      if (keyword_set(invvar)) then begin
         psfvals = [0.496,0.246] ; for FWHM=2.0 pix
         reject_cr, image, invvar, psfvals, rejects, c2fudge=0.8, niter=6, nrejects=nrejects
         if (nrejects GT 0) then begin
            crmask = 0. * image
            crmask[rejects] = 1
            crmask = dilate(crmask, replicate(1,3,3))
            invvar *= (crmask EQ 0)
         endif
      endif
      if (arg_present(fbadpix)) then fbadpix = mean(invvar EQ 0 AND mask EQ 1)
   end

endif 
;
;---------------------------------------------------------------------------
; Correct image with pixel-to-pixel flat-field
;---------------------------------------------------------------------------

if (keyword_set(applypixflat) AND (readimg OR readivar)) then begin
  pp = filepath('', root_dir=getenv('SPECFLAT_DIR'), subdirectory='flats')
  ; First search for files "pixflatave-*.fits*", and if not found then look for "pixflat-*.fits*".
  pixflatname = findopfile('pixflatave-*-'+camname+'.fits*', mjd, pp, silent=silent)
  if (~keyword_set(pixflatname)) then pixflatname = findopfile('pixflat-*-'+camname+'.fits*', mjd, pp, silent=silent)
   
  if (~keyword_set(pixflatname)) then begin
    if (~keyword_set(silent)) then splog, 'WARNING: Pixel flat not found for camera ' + camname
  endif else begin
    if (~keyword_set(silent)) then splog, 'Correcting with pixel flat ' + pixflatname
      
    pixflatimg = mrdfits(djs_filepath(pixflatname, root_dir=pp), /fscale, silent=silent)
    if (total(size(pixflatimg,/dimens) NE size(image,/dimens)) GT 0) then message, 'Dimensions of image and pixel flat differ!'

    ; now get bad pixel mask
    badpixelname = findopfile('badpixels-*-'+camname+'.fits*', mjd, pp, silent=silent)
   
    if (~keyword_set(badpixelname)) then begin
      if (~keyword_set(silent)) then splog, 'WARNING: Badpixels not found for camera ' + camname
    endif else begin
      if (~keyword_set(silent)) then splog, 'Correcting with badpixels ' + badpixelname
  
      badpixelimg = mrdfits(djs_filepath(badpixelname, root_dir=pp), /fscale, silent=silent)
      if (total(size(badpixelimg,/dimens) NE size(image,/dimens)) GT 0) then message, 'Dimensions of image and badpixels differ!'
    
      ; include badpixels into pixflat
      badpixuse=where(badpixelimg ne 0,ct)
      if (ct ne 0) then pixflatimg[badpixuse]=0.0
    endelse

                                ; now check for saturated pixels/bad 
                                ; columns on individual exposure
                                ; compare to neighbor columns, add
                                ; columns with hotpixel trail to pixflatimg 
    
    if (camname eq 'b1' or camname eq 'b2') and flavor eq 'science' then begin
        nxsat=n_elements(image[*,0])
        nysat=n_elements(image[0,*])
        rowsat=20       ;  number of rows to average
        threshsat=10    ;  threshold for column comparison
        satrat=0.05     ;  ratio of column brightness for trail
        for i=2,nxsat-3 do begin
            hp=where(image[i,*] gt 60000 and (pixflatimg[i,*] ge 0.5 or pixflatimg[i,*] eq 0),ct) 
                                ;  column search for hotpixels gt 60000 not caused by pixflat
            if ct ne 0 then begin
                if hp[0] lt nysat/2. and hp[0] +10 +rowsat lt nysat/2.$
                then begin      ;bottom half compare to neighbor columns
                    colleft=mean(image[i-1,hp[0]+10:hp[0]+10+rowsat])
                    colcen=mean(image[i,hp[0]+10:hp[0]+10+rowsat])
                    colright=mean(image[i+1,hp[0]+10:hp[0]+10+rowsat])
                                ; if difference greater than
                                ; threshold, mask to center
                    if hp[0] lt 3 then hp[0]=3 ;fix edge problem
                    if colcen-colleft gt threshsat and colcen-colright gt $ 
                      threshsat and colleft gt 0 and colcen gt 0 and $
                      colright gt 0 and (colcen-colleft)/colcen gt satrat and $
                      (colcen-colright)/colcen gt satrat $ 
                      then pixflatimg[i,hp[0]-3:nysat/2.]=0.0
                endif  
                                ; if hotpixel near center, mask to center
                if hp[0] lt nysat/2. and hp[0] +10 +rowsat ge nysat/2. then pixflatimg[i,hp[0]-3:nysat/2.]=0.0
                
                if max(hp) gt nysat/2. and max(hp)-10-rowsat gt nysat/2. $
                  then begin    ;top half same as above
                    colleft=mean(image[i-1,max(hp)-10-rowsat:max(hp)-10])
                    colcen=mean(image[i,max(hp)-10-rowsat:max(hp)-10])
                    colright=mean(image[i+1,max(hp)-10-rowsat:max(hp)-10])
                    if max(hp) gt nysat-4 then hp=nysat-4 ;fix edge problem
                    if colcen-colleft gt threshsat and colcen-colright gt $ 
                      threshsat and colleft gt 0 and colcen gt 0 and $
                      colright gt 0 and (colcen-colleft)/colcen gt satrat and $
                      (colcen-colright)/colcen gt satrat $ 
                      then pixflatimg[i,nysat/2.:max(hp)+3]=0.0
                endif 
                if max(hp) gt nysat/2. and max(hp)-10 -rowsat le nysat/2. then pixflatimg[i,nysat/2.:max(hp)+3]=0.0
            endif
        endfor
    endif

    if (readimg) then image /= (pixflatimg + (pixflatimg LE 0))
    if (~keyword_set(minflat)) then minflat = 0.0
    if (~keyword_set(maxflat)) then maxflat = 1.0e10
    if (readivar) then invvar = invvar * pixflatimg^2 * (pixflatimg GT minflat) * (pixflatimg LT maxflat)
    pixflatimg = 0 ; clear memory
    badpixelimg = 0 ; clear memory

    ; Add pixflatname to header since it has just been applied
    sxaddpar, hdr, 'PIXFLAT', pixflatname
    if (keyword_set(badpixelname)) then sxaddpar, hdr, 'BADPIXEL', badpixelname

  endelse
endif

;---------------------------------------------------------------------------
; Check for NaN's
;---------------------------------------------------------------------------

; This should never happen, but just in case...

if (readimg OR readivar) then begin
  inan = where(finite(image) EQ 0, nnan)
  if (nnan GT 0) then begin
    if (~keyword_set(silent)) then splog, 'WARNING: Replacing ', nnan, ' NaN values'
    image[inan] = 0
    invvar[inan] = 0
  endif
endif

;---------------------------------------------------------------------------
; Write output files
;---------------------------------------------------------------------------

if keyword_set(bcfile) then sxaddpar, hdr, 'OPBC', bcfile
if keyword_set(configfile) then sxaddpar, hdr, 'OPCONFIG', configfile
if keyword_set(ecalibfile) then sxaddpar, hdr, 'OPECALIB', ecalibfile
sxdelpar, hdr, 'UNSIGNED'

if (keyword_set(outfile)) then begin
  mwrfits, image, outfile, hdr, /create
  mwrfits, invvar, outfile
endif

return
end
;------------------------------------------------------------------------------
