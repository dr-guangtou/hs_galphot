;+
;; NAME:
;   mdrp_combinecameras
;
; PURPOSE:
;   Taken in flux-calibrated individual mgFFrames and combines them
;   together into mgCFrame files.
;
;   Based on the manga prototype routine mlfluxcal.pro
;   for which the relevant part was based on BOSS code
;   spcoadd_v5.
;
; CALLING SEQUENCE:
;   mdrp_combinecameras,filenames,[/linwave]
;
; INPUTS:
;   filenames- 4 element array of filenames pointing to the mgFFrame files
;              for b1/r1/b2/r2 cameras (MUST be in this order)
;
; OPTIONAL INPUTS:
;   /linwave- Set this flag to use MaNGA linear wavelength solution
;             (default is log wavelength steps)
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   ml_mgframeread
;   yanny_par()
;   ml_strreplace()
;   ml_readslitmap()
;   traceset2xy
;   correct_dlam
;   mdrp_divideflat
;   ml_makelabel()
;   combine1fiber
;   fileandpath()
;   ml_mwrfits
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   25-Mar-2014  First written (D. Law)
;   22-Apr-2014  Modified output keywords, handle log or linear solns (D. Law)
;   23-Apr-2014  Inserted explicit wavelength call for ivar apodization
;                in the dichroic region (D. Law)
;   02-May-2014  Output SN2 in each camera to FITS header (D. Law)
;   07-May-2014  Input now has flux calibration applied already (D. Law)
;   15-May-2014  Added quality control bitmask logic (D. Law)
;
;-
pro mdrp_combinecameras,filenames,linwave=linwave

on_error,0
compile_opt idl2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Quality control checks.  If fails return from function
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; mgsframenames must be a 4-element array
nfiles=n_elements(filenames)
if (n_elements(filenames) ne 4) then begin
  splog,strcompress('Not passed 4 frames, skipping!')
  return
endif

; Assumes that 0,1,2,3 is in order b1,r1,b2,r2
camnames=['b1','r1','b2','r2']
for i=0,nfiles-1 do begin
  ml_mgframeread,filenames[i],hdr=hdr

  ; Check that input files read properly
  if (~keyword_set(hdr)) then begin
    splog,strcompress('File '+filenames[i]+' not found, skipping!')
    return
  endif

  ; Check that input file is for the right camera
  if (strtrim(sxpar(hdr,'CAMERAS'),2) ne camnames[i]) then begin
    splog,strcompress('File '+filenames[i]+' is wrong camera, skipping!')
    return
  endif

  ; Check that no DRP status flags set that would prohibit processing 
  failout=0
  drpqual=long(fxpar(hdr,'DRPQUAL'))
  ; If the VALIDFILE bit wasn't set, skip frame
  if ((drpqual and sdss_flagval('MANGA_DRPQUALFLAG','VALIDFILE')) eq 0) then failout=1
  ; If the SKYSUBBAD bit was set, skip frame
  if ((drpqual and sdss_flagval('MANGA_DRPQUALFLAG','SKYSUBBAD')) ne 0) then failout=1

  if (failout) then begin
    splog,strcompress('File '+infile+' has quality problems ('+sdss_flagname('MANGA_DRPQUALFLAG',drpqual,/concat)+'), skipping!')
    return
  endif
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End quality control checks
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


ncam=n_elements(camnames)
nexpvec=lonarr(ncam)
exptimevec=fltarr(ncam)

splog,'Reading science + wavelength data from files:'
splog,'b1 File: '+filenames[0]
splog,'r1 File: '+filenames[1]
splog,'b2 File: '+filenames[2]
splog,'r2 File: '+filenames[3]

zeropoint = 3.5D

; Read the slitmap from the first file
ml_mgframeread,filenames[0],slitmap=slitmap,slithdr=slithdr

; Reading slithdr from the mgFrame is buggy and doesn't work so well
; when trying to write back out again.  Read it from mangacore
plugname=yanny_par(slithdr,'plugmap')
slitname=ml_strreplace(plugname,'plPlugMapM','slitmap')
temp=ml_readslitmap(slitname,hdr=slithdr)

; Figure out number of spectra in each spectrograph
; (nobj1, nobj2)
temp=where(slitmap.spectrographid eq 1,nobj1)
partialslit1=slitmap[temp]
temp=where(slitmap.spectrographid eq 2,nobj2)
partialslit2=slitmap[temp]
; Assign this to a vector for the input frames
nobj=lonarr(ncam)
nobj[0]=nobj1
nobj[1]=nobj1
nobj[2]=nobj2
nobj[3]=nobj2

; Make a big vector of fiberid corresponding to the big arrays
; that will be defined later.  This has length 2*(nobj1+nobj2)
; as is ordered b1,r1,b2,r2
bigfiberid=lonarr(2*(nobj1+nobj2))
for i=0,nobj1-1 do bigfiberid[i]=partialslit1[i].fiberid
for i=0,nobj1-1 do bigfiberid[i+nobj1]=partialslit1[i].fiberid
for i=0,nobj2-1 do bigfiberid[i+nobj1+nobj1]=partialslit2[i].fiberid
for i=0,nobj2-1 do bigfiberid[i+nobj1+nobj1+nobj2]=partialslit2[i].fiberid

; Start by determining sizes of all files
npixarr = lonarr(nfiles)
for ifile=0, nfiles-1 do begin
  ml_mgframeread, filenames[ifile], hdr=objhdr
  npixarr[ifile] = sxpar(objhdr,'NAXIS1')
endfor
npixmax = max(npixarr)

; Loop over the files
for ifile=0, nfiles-1 do begin
  ; Read in all data from this input file.
  ; Note that these are already flux calibrated and in units of
  ; flux/Angstrom
  splog, 'Reading file #', ifile, ': ', filenames[ifile]
  ml_mgframeread, filenames[ifile], objflux=tempflux, objivar=tempivar, $
    mask=temppixmask, wset=tempwset, dispset=tempdispset, $
    skyflux=tempsky, loglam=tempwave, dispimg=tempdispersion, $
    hdr=hdr, adderr=adderr

  ; hdrarr is an array of the headers
  if (ifile EQ 0) then $
    hdrarr = ptr_new(hdr) $
  else $
    hdrarr = [hdrarr, ptr_new(hdr)]

  ; Read header info
  thismjd = sxpar(hdr, 'MJD')
  if (NOT keyword_set(mjdlist)) then mjdlist = thismjd $
  else mjdlist = [mjdlist, thismjd]
  cameras = strtrim(sxpar(hdr, 'CAMERAS'),2)
  expstr = string(sxpar(hdr, 'EXPOSURE'), format='(i8.8)')

  ; Here is the correct conversion from pixels to log-lambda dispersion.
  ; We are converting from the dispersion in units of mgFrame pixel sizes
  ; to the dispersion in units of the new rebinned pixel size, which is
  ; 1e-4 in log-lambda units.
  ; DRL- note that if /linear was set then we'll do the conversion to linear
  ; units at the very end when everything on a common wave grid.
      
  ; this probably should be fixed elsewhere but limit where the fit range
  tempxmax=tempwset.xmax
  tempwset.xmax=(size(tempdispersion,/dimens))[0]-1
  correct_dlam, tempdispersion, 0, tempwset, /inverse
  tempwset.xmax=tempxmax

  dims = size(tempflux, /dimens)
  npix = dims[0]
  nfib = dims[1]

  ; Determine if this is a blue or red spectrum
  icam = (where(cameras EQ camnames))[0]
  if (icam EQ -1) then $
    message, 'Unknown camera ' + cameras
  ;nexpvec[icam] = nexpvec[icam] + 1
  ;exptimevec[icam] = exptimevec[icam] + sxpar(hdr, 'EXPTIME')

  expnum = sxpar(hdr, 'EXPOSURE')

  ; Apodize the errors for the dichroic overlap region.
  ; DRL- BOSS previously did this for a window of the first 100 pixels,
  ; which was pretty useless for us (and not the dichroic region for b1/b2 anyway)
  ; Rewriting to do this as a function of wavelength explicitly.
  if ((cameras eq 'b1')or(cameras eq 'b2')) then begin
    lam_max=alog10(6300.)
    index=where(tempwave gt lam_max,nindex)
    if (nindex gt 0) then tempivar[index]=0.
  endif
  if ((cameras eq 'r1')or(cameras eq 'r2')) then begin
    lam_min=alog10(5900.)
    index=where(tempwave lt lam_min,nindex)
    if (nindex gt 0) then tempivar[index]=0.
  endif

  ; Concatenate data from all images
  if (ifile EQ 0) then begin
    ; Construct the image arrays
    nbig=2*(nobj1+nobj2)
    flux = make_array(npixmax,nbig,type=size(tempflux,/type))
    fluxivar = make_array(npixmax,nbig,type=size(tempivar,/type))
    wave = make_array(npixmax,nbig,type=size(tempwave,/type))
    dispersion = make_array(npixmax,nbig,type=size(tempdispersion,/type))
    pixelmask = make_array(npixmax,nbig,type=size(temppixmask,/type))
    skyflux = make_array(npixmax,nbig,type=size(tempsky,/type))

    ; Append as vectors...
    camerasvec = cameras
    label = ml_makelabel(hdr)
    filenum = lonarr(nfib) + ifile
  endif else begin
    ; Append as vectors...
    camerasvec = [camerasvec, cameras]
    label = [label, ml_makelabel(hdr)]
    filenum = [filenum, lonarr(nfib) + ifile]
  endelse

  ; Figure out where we should start in y space for sticking
  ; info into the big arrays
  if (ifile eq 0) then ystart=0
  if (ifile eq 1) then ystart=nobj[0]
  if (ifile eq 2) then ystart=nobj[0]+nobj[1]
  if (ifile eq 3) then ystart=nobj[0]+nobj[1]+nobj[2]

  flux[0:npixarr[ifile]-1,ystart:ystart+nobj[ifile]-1] = tempflux
  fluxivar[0:npixarr[ifile]-1,ystart:ystart+nobj[ifile]-1] = tempivar
  wave[0:npixarr[ifile]-1,ystart:ystart+nobj[ifile]-1] = tempwave
  ; Pad the wavelengths with reasonable values
  if (npixarr[ifile] LT npixmax) then begin
    dwave = tempwave[npixarr[ifile]-1,*] - tempwave[npixarr[ifile]-2,*]
    addwave = tempwave[npixarr[ifile]-1,*] $
      ## (1+lonarr(npixmax-npixarr[ifile])) $
      + dwave ## (1+lindgen(npixmax-npixarr[ifile]))
    wave[npixarr[ifile]:npixmax-1,ystart:ystart+nobj[ifile]-1] = addwave
  endif
  dispersion[0:npixarr[ifile]-1,ystart:ystart+nobj[ifile]-1] = tempdispersion
  pixelmask[0:npixarr[ifile]-1,ystart:ystart+nobj[ifile]-1] = temppixmask
  skyflux[0:npixarr[ifile]-1,ystart:ystart+nobj[ifile]-1] = tempsky
endfor

tempflux = 0
tempivar = 0
tempwave = 0
tempdispersion = 0
temppixmask = 0
tempsky = 0

; Construct output data structures, including the wavelength scale
totalpix = (size(flux, /dimens))[0]
nonzero = where(fluxivar GT 0.0)
minfullwave = min(wave[nonzero])
maxfullwave = max(wave[nonzero])

; Set output wavelengths
; LOG spacing unless /linwave keyword set
if (~keyword_set(linwave)) then finalwave=alog10(ml_setwcalib()) $
else finalwave=alog10(ml_setwcalib(/linear))

nfinalpix=n_elements(finalwave)
nfiber=nobj1+nobj2

finalflux = fltarr(nfinalpix, nfiber)
finalivar = fltarr(nfinalpix, nfiber)
finalandmask = lonarr(nfinalpix, nfiber)
finalormask = lonarr(nfinalpix, nfiber)
finaldispersion = fltarr(nfinalpix, nfiber)
finalsky = fltarr(nfinalpix, nfiber)

; Combine each fiber, one at a time along the slit
for ifiber=0, nfiber-1 do begin
  ; Find where fiberid=ifiber+1 is in the big structure
  indx=where(bigfiberid eq ifiber+1)

  if (indx[0] NE -1) then begin
    temppixmask = pixelmask[*,indx]
    ; DRL- call this with binsz=1d-4
    ; It isn't really true, but it's effectively what BOSS did,
    ; and while it gives a few more artifacts on the B-spline most
    ; of these are covered by the mask.  NOT doing this sometimes
    ; gives bad behaviour on emission line features.
    combine1fiber, wave[*,indx], flux[*,indx], fluxivar[*,indx], $
      finalmask=temppixmask, indisp=dispersion[*,indx], $
      skyflux=skyflux[*,indx], $
      newloglam=finalwave, newflux=bestflux, newivar=bestivar, $
      andmask=bestandmask, ormask=bestormask, newdisp=bestdispersion, $
      newsky=bestsky, $
      nord=nord, binsz=1.0d-4, bkptbin=bkptbin, maxsep=maxsep, $
      maxiter=50, upper=3.0, lower=3.0, maxrej=1

    finalflux[*,ifiber] = bestflux
    finalivar[*,ifiber] = bestivar
    finalandmask[*,ifiber] = bestandmask
    finalormask[*,ifiber] = bestormask
    finaldispersion[*,ifiber] = bestdispersion
    finalsky[*,ifiber] = bestsky

    ; The following adds the COMBINEREJ bit to the input pixel masks
    pixelmask[*,indx] = temppixmask
  endif else begin
    splog, 'Fiber', ifiber+1, ' NO DATA'
    finalandmask[*,ifiber] = pixelmask_bits('NODATA')
    finalormask[*,ifiber] = pixelmask_bits('NODATA')
  endelse
endfor

; Modify the 1st file's header to use for the combined plate header.
bighdr = *hdrarr[0]

; Create a combined mask array out of the AND and OR mask
; arrays, in some way physically meaningful for MaNGA.
; In the dichroic region some important flags are lost by
; the AND mask, but the OR mask flags way too much.
; Added by DRL on April 2 2014
mangamask=finalandmask
; Define the FULLREJECT and COMBINEREJ masks
; Initialize them to zero
fullreject=mangamask
fullreject[*]=0.
combinerej=fullreject
; Identify where fullreject but NOT lowflat are set in the OR mask
; (this is because we don't want to fullreject things in the dichroic region
; that don't actually contribute to a given pixel b/c the flat is so low)
index=where( ((finalormask and sdss_flagval('SPPIXMASK','FULLREJECT')) ne 0) and ((finalormask and sdss_flagval('SPPIXMASK','LOWFLAT')) eq 0),nindex)
if (nindex ne 0) then fullreject[index]=pixelmask_bits('FULLREJECT')
; Identify where combinerej is set in the OR mask
index=where((finalormask and sdss_flagval('SPPIXMASK','COMBINEREJ')) ne 0,nindex)
if (nindex ne 0) then combinerej[index]=pixelmask_bits('COMBINEREJ')
; Add these bits to the final mangamask
mangamask = mangamask OR fullreject
mangamask = mangamask OR combinerej

; Clear memory
wave = 0
flux = 0
fluxivar = 0
temppixmask = 0
dispersion = 0
skyflux = 0
fullreject = 0
combinerej = 0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Write the corrected spCFrame files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Concatenate the SN2, DRPQUAL,BADIFU, PSTAT information
drpqual=0
badifu=''
for ifile=0,nfiles-1 do begin
  cameras = strtrim(sxpar(*hdrarr[ifile], 'CAMERAS'),2)
  thisdrpqual=long(sxpar(*hdrarr[ifile], 'DRPQUAL'))

  case cameras of
    'b1': begin
            sxaddpar,bighdr,'B1SN2',sxpar(*hdrarr[ifile],'FRAMESN2'),' SN2 in B1 CAMERA'
            sxaddpar,bighdr,'B1PSTAT',sxpar(*hdrarr[ifile],'PSTAT'),' PSTAT in B1 CAMERA'
          end
    'r1': begin
            sxaddpar,bighdr,'R1SN2',sxpar(*hdrarr[ifile],'FRAMESN2'),' SN2 in R1 CAMERA',after='B1SN2'
            sxaddpar,bighdr,'R1PSTAT',sxpar(*hdrarr[ifile],'PSTAT'),' PSTAT in R1 CAMERA',after='B1PSTAT'
          end
    'b2': begin
            sxaddpar,bighdr,'B2SN2',sxpar(*hdrarr[ifile],'FRAMESN2'),' SN2 in B2 CAMERA',after='B1SN2'
            sxaddpar,bighdr,'B2PSTAT',sxpar(*hdrarr[ifile],'PSTAT'),' PSTAT in B2 CAMERA',after='B1PSTAT'
          end
    'r2': begin
            sxaddpar,bighdr,'R2SN2',sxpar(*hdrarr[ifile],'FRAMESN2'),' SN2 in R2 CAMERA',after='R1SN2'
            sxaddpar,bighdr,'R2PSTAT',sxpar(*hdrarr[ifile],'PSTAT'),' PSTAT in R2 CAMERA',after='R1PSTAT'
          end
  endcase

  drpqual=drpqual or thisdrpqual

  if (cameras eq 'b1') then badifu+= sxpar(*hdrarr[ifile],'BADIFU')
  if (cameras eq 'b2') then badifu+= sxpar(*hdrarr[ifile],'BADIFU')
endfor
sxaddpar, bighdr, 'DRPQUAL', drpqual, 'DRP quality bitmask'
sxaddpar, bighdr, 'BADIFU', badifu, 'Harness names of missing/bad ifus'

; Remove header cards that were specific to this first exposure
; (where we got the header) or otherwise irrelevant.

ncoeff = sxpar(bighdr, 'NWORDER')
for i=2, ncoeff-1 do sxdelpar, bighdr, 'COEFF'+strtrim(string(i),2)

sxdelpar, bighdr, ['SPA', 'IPA', 'IPARATE']
sxdelpar, bighdr, 'REQTIME'
sxdelpar, bighdr, 'QUALITY'
sxdelpar, bighdr, 'FILENAME'
sxdelpar, bighdr, 'SEQID'
sxdelpar, bighdr, 'POINTING'
sxdelpar, bighdr, 'DARKTIME'
sxdelpar, bighdr, 'CAMERAS'
sxdelpar, bighdr, 'PLUGMAPO'
for i=1, 4 do sxdelpar, bighdr, 'GAIN'+strtrim(string(i),2)
for i=1, 4 do sxdelpar, bighdr, 'RDNOISE'+strtrim(string(i),2)
sxdelpar, bighdr, ['CAMCOL', 'CAMROW']
sxdelpar, bighdr, ['AMPLL', 'AMPLR', 'AMPUL', 'AMPUR']
sxdelpar, bighdr, ['FFS', 'FF', 'NE', 'HGCD']
sxdelpar, bighdr, ['SPEC1', 'SPEC2']
sxdelpar, bighdr, 'NBLEAD'
sxdelpar, bighdr, 'PIXFLAT'
sxdelpar, bighdr, 'PIXBIAS'
sxdelpar, bighdr, 'FLATFILE'
sxdelpar, bighdr, 'ARCFILE'
sxdelpar, bighdr, 'OBJFILE'
sxdelpar, bighdr, 'FRAMESN2'
sxdelpar, bighdr, 'PSTAT'
sxdelpar, bighdr, 'DEREDSN2'

; Average together some of the fields from the individual headers.

cardname = [ 'AZ', 'ALT', 'TAI', 'WTIME', 'AIRTEMP', 'DEWPOINT', $
    'DEWDEP', 'DUSTA', 'DUSTB', 'DUSTC', 'DUSTD', 'GUSTS', 'HUMIDITY', $
    'HUMIDOUT', 'PRESSURE', 'WINDD', 'WINDS', 'TEMP01', 'TEMP02', $
    'TEMP03', 'TEMP04', 'HELIO_RV', 'SEEING20', 'SEEING50', 'SEEING80', $
    'RMSOFF20', 'RMSOFF50', 'RMSOFF80', 'XCHI2', 'SKYCHI2', $
    'WSIGMA', 'XSIGMA' ]
sxcombinepar, hdrarr, cardname, bighdr, func='average'
sxcombinepar, hdrarr, 'TAI-BEG', bighdr, func='min'
sxcombinepar, hdrarr, 'TAI-END', bighdr, func='max'
sxcombinepar, hdrarr, 'XCHI2', bighdr, func='max', outcard='XCHI2MAX', $
    after='XCHI2'
sxcombinepar, hdrarr, 'XCHI2', bighdr, func='min', outcard='XCHI2MIN', $
    after='XCHI2'
sxcombinepar, hdrarr, 'SKYCHI2', bighdr, func='max', outcard='SCHI2MAX', $
    after='SKYCHI2'
sxcombinepar, hdrarr, 'SKYCHI2', bighdr, func='min', outcard='SCHI2MIN', $
    after='SKYCHI2'
sxcombinepar, hdrarr, 'WSIGMA', bighdr, func='max', outcard='WSIGMAX', $
    after='WSIGMA'
sxcombinepar, hdrarr, 'WSIGMA', bighdr, func='min', outcard='WSIGMIN', $
    after='WSIGMA'
sxcombinepar, hdrarr, 'XSIGMA', bighdr, func='max', outcard='XSIGMAX', $
    after='XSIGMA'
sxcombinepar, hdrarr, 'XSIGMA', bighdr, func='min', outcard='XSIGMIN', $
    after='XSIGMA'

; Add the NGUIDE keywords for all headers of one flavor of CAMERAS
; (e.g., for all the 'b1' exposures if the first frame is 'b1'.)
cardname = 'NGUIDE'
sxcombinepar, hdrarr[0], cardname, bighdr, func='total'
cameras0 = sxpar(*(hdrarr[0]), 'CAMERAS')
for ihdr=1, n_elements(hdrarr)-1 do begin
  if (sxpar(*(hdrarr[ihdr]), 'CAMERAS') EQ cameras0) then $
      sxcombinepar, hdrarr[ihdr], cardname, bighdr, func='total'
endfor

; Add keycards
sxaddpar, bighdr, 'NAXIS1', n_elements(bestflux)
sxaddpar, bighdr, 'NAXIS2', nfiber

spawn, 'uname -n', uname
sxaddpar, bighdr, 'UNAME', uname[0]

; Compute the fraction of bad pixels in total, and on each spectrograph.
; Bad pixels are any with SKYMASK(INVVAR)=0, excluding those where
; the NODATA bit is set in the pixel mask.
ifib1 = where(slitmap.spectrographid EQ 1, nfib1)
ifib2 = where(slitmap.spectrographid EQ 2, nfib2)
qbadpix = skymask(finalivar, finalandmask, finalormask) EQ 0 $
    AND (finalandmask AND pixelmask_bits('NODATA')) EQ 0
if (nfib1 GT 0) then $
    fbadpix1 = total(qbadpix[*,ifib1]) / (nfib1 * nfinalpix) $
else $
    fbadpix1 = 0
if (nfib2 GT 0) then $
    fbadpix2 = total(qbadpix[*,ifib2]) / (nfib2 * nfinalpix) $
else $
    fbadpix2 = 0
if (nfib1 GT 0 AND nfib2 GT 0) then $
    fbadpix = total(qbadpix[*,[ifib1,ifib2]]) / ((nfib1+nfib2) * nfinalpix) $
else if (nfib1 GT 0) then $
    fbadpix = fbadpix1 $
else if (nfib2 GT 0) then $
    fbadpix = fbadpix1 $
else $
    fbadpix = 0

sxaddpar, bighdr, 'FBADPIX', fbadpix, ' Fraction of bad pixels'
sxaddpar, bighdr, 'FBADPIX1', fbadpix1, ' Fraction of bad pixels on spectro-1'
sxaddpar, bighdr, 'FBADPIX2', fbadpix2, ' Fraction of bad pixels on spectro-2'

; Add wavelength and flux information keywords
sxaddpar, bighdr, 'BSCALE', 1., ' Flux unit scaling'
sxaddpar, bighdr, 'BZERO', 0., ' Flux zeropoint'
sxaddpar, bighdr, 'BUNIT', '1E-17 erg/s/cm^2/Ang/fiber', ' Flux units are per fiber'

; Add keywords for IRAF-compatability
if (keyword_set(linwave)) then begin
    sxaddpar, bighdr, 'CTYPE1','WAVE-WAV'
    sxaddpar, bighdr, 'CRPIX1',1,' Starting pixel (1-indexed)'
    sxaddpar, bighdr, 'CRVAL1',min(10.^finalwave),' Central wavelength of first pixel'
    sxaddpar, bighdr, 'CD1_1',abs(10.^finalwave[1]-10.^finalwave[0]),' Linear dispersion per pixel'
endif else begin
    sxaddpar, bighdr, 'DC-FLAG',1,' Log-linear flag'
    sxaddpar, bighdr, 'CTYPE1','WAVE-LOG'
    sxaddpar, bighdr, 'CRPIX1',1,' Starting pixel (1-indexed)'
    sxaddpar, bighdr, 'CRVAL1',min(finalwave),' Central wavelength (log10) of first pixel'
    sxaddpar, bighdr, 'CD1_1',abs(finalwave[1]-finalwave[0]),' Log10 dispersion per pixel'
endelse

; Define the output filename from the input file
tempname=fileandpath(filenames[0],path=path)
tempname=(strsplit(tempname,'.',/extract))[0]
; This is the exposure number
tempname=(strsplit(tempname,'-',/extract))[2]
if (keyword_set(linwave)) then outname=path+'mgCFrame-'+tempname+'-LIN.fits' $
else outname=path+'mgCFrame-'+tempname+'-LOG.fits'

; If desired output is linear sampling, convert the dispersion
; extension to per-Angstrom instead of per-1e-4 log10(wave)
if (keyword_set(linwave)) then begin
  delta=1e4*(alog10(10.^finalwave+0.5)-alog10(10.^finalwave-0.5))
  finaldispersion=finaldispersion/rebin(reform(delta,[n_elements(finalwave),1]),[n_elements(finalwave),nfiber])
endif

; Write sky-subtracted spectra to disk following the MaNGA data model
ml_mwrfits, dummyext, outname,hdr=bighdr, /create ; Blank ext 0 with full header
ml_mwrfits, finalflux, outname, extname='FLUX' ; Flux
ml_mwrfits, finalivar, outname, extname='IVAR' ; Inverse variance
ml_mwrfits, mangamask, outname, extname='MASK' ; Pixel mask
ml_mwrfits, 10.^finalwave, outname, extname='WAVE' ; Wavelength solution
ml_mwrfits, finaldispersion, outname, extname='DISP' ; Dispersion
ml_mwrfits, slitmap, outname, hdr=slithdr, extname='SLITMAP'
ml_mwrfits, finalsky, outname, extname='SKY' ; Sky image

; gzip output file
spawn, ['gzip', '-f', outname], /noshell

return
end
