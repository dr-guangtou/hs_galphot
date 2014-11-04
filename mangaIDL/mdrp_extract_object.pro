;+
; NAME:
;   mdrp_extract_object
;
; PURPOSE:
;   Performs all object extraction tasks
;      0) Locate bright fibers, and test image background
;      1) 4 Step Optimal extraction
;
; CALLING SEQUENCE:
;   mdrp_extract_object, outname, objhdr, image, invvar, plugsort, wset, $
;    xarc, lambda, xtrace, fflat, fibermask, proftype=, color=, $
;    [ widthset=, dispset=, skylinefile=, plottitle=, superflatset=, $
;    /splitsky, ccdmask= , slithdr=]
;
; INPUTS:
;   outname    - Name of outputs FITS file
;   objhdr     - Header cards from object image
;   image      - Object image [nx,ny]
;   invvar     - Inverse Variance of object image [nx,ny]
;   plugsort   - Plugmap structure for [ntrace] spectra
;   wset       - Wavelength solution from arc line spectra
;   xarc       - centroids measured in arc line spectra
;   lambda     - air wavelengths corresponding to xarc
;   xtrace     - spatial traces from flat field
;   fflat      - 1d flat field vectors
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER]
;   proftype   - Which type of profile should we use, (default=1 Gaussian)
;   superflatset- If present, then divide by median superflat! ???
;   slithdr    - Yanny header from the slitmap
;   color      - ???
;   widthset   - ???
;   dispset    - ???
;
; REQUIRED KEYWORDS:
;   color      - camera color (red or blue)
;
; OPTIONAL KEYWORDS:
;   plottitle  - Prefix for titles in QA plots.
;   ccdmask    - If set, then use this to set some pixel values in pixmask
;
; OUTPUTS:
;   A fits file is output in outname, which contains
;      FLOAT flux [NX,NTRACE]
;      FLOAT flux_invvar [NX,NTRACE]
;      LONG  manga pixelmask [NX,NTRACE]
;      STRUCT vacuum wavelengths
;      STRUCT wavelength dispersion
;      STRUCT slitmap [NTRACE]
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   mdrp_divideflat
;   djs_median()
;   djs_oplot
;   djs_plot
;   extract_boxcar()
;   extract_image
;   fibermask_bits()
;   fitsn()
;   fitvacset()
;   get_tai
;   heliocentric()
;   mwrfits
;   pixelmask_bits()
;   qaplot_scatlight
;   qaplot_skydev
;   spadd_guiderinfo
;   splog
;   sxaddpar
;   sxpar()
;   traceset2xy
;   find_whopping()
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   24-Jan-2000  Written by S. Burles, Chicago
;   26-Jul-2001  Also pass proftype and superflatset
;   29-Mar-2011  Switching to pure IDL bundle-wise extraction (A. Bolton, U. Utah)
;   15-Aug-2011  Enabling split sky model across halves of CCD. (A. Bolton, U. Utah)
;      Nov-2012  Modified for MaNGA survey (B. Cherinka, Toronto)
;   15-Jul-2013  Separated out from Boss (renamed mdrp_), removed all sky subtraction stuff
;   15-May-2014  Added quality control bitmask logic (David Law)
;-
;------------------------------------------------------------------------------


pro mdrp_extract_object, outname, objhdr, image, invvar, slitmap, wset, $
 xarc, lambda, xtrace, fflat, fibermask, color=color, proftype=proftype, skylinefile=skylinefile, $
 widthset=widthset, dispset=dispset, plotname=plotname, plottitle=plottitle, superflatset=superflatset, $
 ccdmask=ccdmask, slithdr=slithdr, VISUAL=VISUAL, SURVEY=survey, fiberparam=fiberparam, writefiles=writefiles
 
   on_error,0
  compile_opt idl2
  
  if (keyword_set(visual)) then begin
    set_plot,'x'
    device,decomposed=0
    loadct,0,/silent
    !p.multi=[0,1,2]
    !p.thick=0 & !x.thick=0 & !y.thick=0 & !p.charthick=0
  endif else if (keyword_set(plotname)) then begin
    set_plot,'ps'
    origp=!p
    origx=!x
    origy=!y
    !p.multi=[0,1,2]
    !p.thick=0 & !x.thick=0 & !y.thick=0 & !p.charthick=0
    dfpsplot, plotname, /color
  endif

   configuration=obj_new('configuration', sxpar(objhdr, 'MJD'))

   objname = strtrim(sxpar(objhdr,'OBJFILE'),2) 
   flavor  = strtrim(sxpar(objhdr,'FLAVOR'),2) 
   camera  = strtrim(sxpar(objhdr,'CAMERAS'),2) 
   spectrographid=fix(strmid(camera,1,1))
   partialslit=slitmap[where(slitmap.spectrographid eq spectrographid)]
   nfiber = fiberparam.nfiber
   radius = fiberparam.radius
   totfiber = fiberparam.totfiber
   nbundle = (long(total(fiberparam.bundlegap NE 0)))[0]
   bundleid = fiberparam.bundleid
 
   ;------------------
   ; Identify very bright objects
   ; Do a boxcar extraction, and look for fibers where the median counts are 10000 ADU per row.

   fextract = ml_extract_boxcar(image*(invvar GT 0), xtrace, nfiber=nfiber, radius=radius)
   scrunch = djs_median(fextract, 1) ; Find median counts/row in each fiber
   whopping = find_whopping(scrunch, 10000.0, whopct)
   scrunch_sort = sort(scrunch)
   i5 = n_elements(scrunch)/20
   i95 = i5 * 19

   splog, 'Whopping fibers: ', whopping
   splog, 'Median counts in all fibers = ', djs_median(scrunch)
   splog, 'Number of bright fibers = ', whopct
   
   if (whopct GT 20) then begin
      splog, 'WARNING: Disable whopping terms ' + objname
      whopping = -1
      whopct = 0
   endif   

   ;------------------------------------------------------------
   ;  Check for bad pixels within 3 pixels of trace
   badcheck = ml_extract_boxcar((invvar LE 0), xtrace, radius=2.5, nfiber=nfiber)
   badplace = where(badcheck GT 0)

   nx = (size(fextract,/dim))[0] 
   ny = (size(fextract,/dim))[1] 
   pixelmask = lonarr(nx,ny)

   badcolumns = where(total(badcheck GT 0,1) GT 0.45 * nx)
   if (badplace[0] NE -1) then pixelmask[badplace] = pixelmask[badplace] OR pixelmask_bits('NEARBADPIXEL')
   if (badcolumns[0] NE -1) then fibermask[badcolumns] = fibermask[badcolumns] OR pixelmask_bits('MANYBADCOLUMNS')

   if (whopping[0] NE -1) then begin
      ; Set the mask bit for whopping fibers themselves
      fibermask[whopping] = fibermask[whopping] OR pixelmask_bits('WHOPPER')

      ; Set the mask bit for fibers near whopping fibers, excluding the whopping fibers themselves.  Note that a fiber could still have both
      ; WHOPPER and NEARWHOPPER set if it is both bright and near another bright fiber.
      wp = [whopping - 2 , whopping -1, whopping+1 , whopping+2]
      wp = wp[ where(wp GE 0 AND wp LT ny) ]
      fibermask[wp] = fibermask[wp] OR pixelmask_bits('NEARWHOPPER')
   endif

   ;-----
   ; Inherit any mask bits from the ccdmask, by setting the pixmask bits for anything that would be hit in a boxcar extraction
   if (keyword_set(ccdmask)) then begin
      for ibit=0, 31 do begin
         thischeck = ml_extract_boxcar((ccdmask AND 2L^ibit) NE 0, xtrace, radius=2.5, nfiber=nfiber)
         pixelmask = pixelmask OR (2L^ibit * (thischeck GT 0))
      endfor
   endif

   ;-----------------------------------------------------------------------
   ;  This is a kludge to fix first and last column ???
   ;-----------------------------------------------------------------------
  if (configuration->extract_object_fixcolumns()) then begin
     image[0,*] = image[0,*]*0.7
     image[nx-1,*] = image[nx-1,*]*0.7
  endif

   ;
   ;  First we should attempt to shift trace to object flexure
   xnow = mdrp_match_trace(image, invvar, xtrace,radius=radius, nfiber=nfiber)
   bestlag = median(xnow-xtrace)

   splog, 'Shifting traces by match_trace ', bestlag
   if (abs(bestlag) GT 1.0) then begin
      splog, 'WARNING: pixel shift is large!'
   endif

   highrej = 10  ; just for first extraction steps
   lowrej = 10   ; just for first extraction steps
                 ; We need to check npoly with new scattered light backgrounds
   npoly = 16    ; maybe more structure, lots of structure
   nrow = (size(image))[2]
   yrow = lindgen(nrow) 
   nfirst = n_elements(yrow)

   splog, 'Extracting frame '+objname+' with 4 step process'

   traceset2xy, widthset, xx, sigma2
   ntrace = (size(sigma2,/dimens))[1]
   wfixed = [1,1] ; Fit gaussian height + width (fixed center position)
   nterms = n_elements(wfixed)

   ;-----------------------------------------------------------------------
   ;  Now, subtract halo image and do final extraction with all rows
   ;-----------------------------------------------------------------------

  ; Usually the scattered light is taken out well enough by the extract_image
  ; routine.  However, in bright time there can be a higher glow
  ; in the sdssproc'ed image that isn't well subtracted by normal routine.
  ; We must also be careful not to use this routine on calibration lamps,
  ; as these leak light from the brightest regions.  Therefore look for
  ; cases where the total counts are low, but background counts high.

  ; Use fixed breakpoint chosen to sample the BOSS cameras well
  bkpt=[0,150,300,700,1100,1500,1900,2300,2700,3100,3500,3800,3950,4100]
  ; Define a mask of the 150 pixel wide boundary of the detector
   boundmask=intarr(size(image,/dimension))
   boundmask[150:(size(image,/dimension))[0]-150,150:(size(image,/dimension))[1]-150]=1
   theedge=where(boundmask eq 0)
  scat_warning=0
  if ((median(image) lt 300)and(median(image[theedge]) gt 30)) then begin
    splog,'WARNING- bad scattered light, applying secondary scattered light model'
    ; Scat_warning=1 if this routine runs successfully, =2 if it fails
    image=mdrp_scattered(image,invvar,xnow,scatmodel=scatmodel,bkpt=bkpt,nbin=5,nord=4,nreq=400,maskrad=3,upper=5.,lower=5.,badthresh=1.5,errflag=scat_warning)
  endif

   ; (6) Final extraction
   splog, 'Step 6: Final Object extraction'

      ; DRL set these parameters May 22 2014
      ; These give truly flat 30-sec flatfields?, and also good
      ; flat sky-subtracted data for the dark all-sky plate
      ; 7341-56693.
      highrej = 8
      lowrej = 5
      wfixed = [1] ; Fit gaussian amplitude only (not width or centroid)
      nterms = n_elements(wfixed)
      reject=[0.1, 0.6, 0.6]
      npoly = 2L

    extract_image, image, invvar, xnow, sigma2, flux, fluxivar, proftype=proftype, wfixed=wfixed, ansimage=ansimage3, highrej=highrej, lowrej=lowrej, npoly=npoly, $
      chisq=chisq, ymodel=ymodel, pixelmask=pixelmask, reject=reject, /relative, visual=visual, survey=survey

    ;write out post-extracted flux (uncalibrated) and ymodel from science extraction
    if keyword_set(writefiles) then begin
      sciextname = ml_strreplace(outname,'mgFrame','mgFrame-flux')
      mwrfits, flux, sciextname, /create
      scimodname=ml_strreplace(outname,'mgFrame','mgFrame-model')
      mwrfits, ymodel, scimodname, /create  
      spawn, ['gzip', '-f', sciextname], /noshell
      spawn, ['gzip', '-f', scimodname], /noshell           
    endif
    
   ;----------------------------------------------------------------------
   ; Can we find cosmic rays by looking for outlandish ansimage ratios???
   ; a = where(ansimage[lindgen(ntrace)*nterms, *] LT (-2*ansimage[lindgen(ntrace)*nterms+1, *])

   ;------------------
   ; QA chisq plot for fit calculated in extract image (make QAPLOT ???)
   xaxis = lindgen(n_elements(chisq)) + 1
   ymax = 2.*median(chisq)
   djs_plot, xaxis, chisq, xrange=[0,N_elements(chisq)], xstyle=1, yrange=[0,ymax], ystyle=1, xtitle='Row number',  ytitle = '\chi^2', title=plottitle+'Extraction chi^2 for '+objname

   djs_oplot, !x.crange, [1,1]
;x   djs_oplot, yrow, secondchisq[yrow], color='blue'
;x   djs_oplot, yrow, firstchisq[yrow], color='green'

   xyouts, 100, 0.05*!y.crange[0]+0.95*!y.crange[1], 'BLACK = Final chisq extraction'
;x   xyouts, 100, 0.08*!y.crange[0]+0.92*!y.crange[1], 'BLUE = Initial-scattered chisq extraction'
;x   xyouts, 100, 0.08*!y.crange[0]+0.89*!y.crange[1], 'GREEN = Initial chisq extraction'

   ;------------------
   ; Flat-field the extracted object fibers with the global flat
   mdrp_divideflat, flux, fflat, invvar=fluxivar, /quiet
 
   pixelmask = pixelmask OR ((fflat LT 0.5) * pixelmask_bits('LOWFLAT'))

   ;------------------
   ; Tweak up the wavelength solution to agree with the sky lines.
   ; xshet contains polynomial coefficients to shift arc to sky line frame.

   locateskylines, skylinefile, flux, fluxivar, wset, xarc, arcshift=arcshift, xsky=xsky, skywaves=skywaves, skyshift=skyshift, radius=radius, nfiber=nfiber

   if keyword_set(VISUAL) then getwindow,/open 
   qaplot_skyshift, wset, xsky, skywaves, skyshift, title=plottitle+'Sky Line Deviations for '+objname

   if ~keyword_set(arcshift) then splog, 'WARNING: Cannot shift to sky lines'

   ;------------------
   ; Fit for the widths of the sky lines (relative to those of the arcs)

   ; The values in XSKY are noisy measurements, so replace them with
   ; the predicted positions based upon the arc solution.
   ; We should also apply the arcshift here, but I haven't yet ???

   xsky = transpose(traceset2pix(wset, alog10(skywaves)))
   skydispset = skyline_dispersion(flux, fluxivar, xsky, iskies, dispset)
   splog, 'Not applying skyline adjusted line widths'

   ;------------------
   ; Apply heliocentric correction
   ; Correct LAMBDA, which is used to shift to vacuum wavelengths.

   helio=0.0
   ra = sxpar(objhdr, 'RA', count=ct_ra)
   dec = sxpar(objhdr, 'DEC', count=ct_dec)
   if (ct_ra NE 1 OR ct_dec NE 1) then splog, 'WARNING: Missing RA and/or DEC from header'

   ;--------------------------------------------------------
   ; Call standard proc to determine time-stamps

   get_tai, objhdr, tai_beg, tai_mid, tai_end

   ; Set TAI equal to the time half-way through the exposure
   ; If all these keywords are present in the header, they will be either type FLOAT or DOUBLE.  Note that SDSS will put NaN in the header for these values if they are unknown.
   if ( size(ra, /tname) NE 'INT' AND size(dec, /tname) NE 'INT' AND size(tai_mid, /tname) NE 'INT' AND finite(ra) AND finite(dec) AND finite(tai_mid) ) then begin
      helio = heliocentric(ra, dec, tai=tai_mid)
      splog, 'Heliocentric correction = ', helio, ' km/s'
      sxaddpar, objhdr, 'HELIO_RV', helio, ' Heliocentric correction (added to velocities)'
   endif else begin
      splog, 'WARNING: Header info not present to compute heliocentric correction'
   endelse
   if (size(tai_mid, /tname) EQ 'INT' OR finite(tai_mid) EQ 0) then begin
      splog, 'WARNING: Header info not present to compute airmass correction to sky level'
      tai_mid = 0
   endif

   ;------------------
   ; Shift to skylines and fit to vacuum wavelengths
   vacset = fitvacset(xarc, lambda, wset, arcshift, helio=helio, airset=airset)
   ; No longer make the following QA plot ???
   ;qaplot_skydev, flux, fluxivar, vacset, plugsort, color, title=plottitle+objname
   sxaddpar, objhdr, 'VACUUM', 'T', ' Wavelengths are in vacuum'

   ;------------------
   ;  If present, reconstruct superflat and normalize
   sxaddpar, objhdr, 'SFLATTEN', 'F', ' Superflat has not been applied'
   superfit = 0

   if keyword_set(superflatset) AND keyword_set(airset) then begin
     superfit = float(smooth_superflat(superflatset, airset, plottitle=plottitle+'Smooth superflat for '+objname))
     if keyword_set(superfit) then begin
       mdrp_divideflat, flux, superfit, invvar=fluxivar, /quiet
       sxaddpar, objhdr, 'SFLATTEN', 'T', ' Superflat has been applied'
     endif
   endif  

   ; Save current pixelmask for later MaNGA use
   mangamask=pixelmask

   ;----------
   ; Add keywords to object header
   objhdr=mdrp_makeheader(head=objhdr,/drp2d)
   if (keyword_set(osigma)) then sxaddpar, objhdr, 'OSIGMA',  sigma, ' Original guess at spatial sigma in pix'
   sxaddpar, objhdr, 'PREJECT', reject[1], ' Profile area rejection threshold'
   sxaddpar, objhdr, 'LOWREJ', lowrej, ' Extraction: low rejection'
   sxaddpar, objhdr, 'HIGHREJ', highrej, ' Extraction: high rejection'
   sxaddpar, objhdr, 'SCATPOLY', npoly, ' Extraction: Order of scattered light polynomial'
   sxaddpar, objhdr, 'PROFTYPE', proftype, ' Extraction profile: 1=Gaussian'
   sxaddpar, objhdr, 'NFITPOLY', nterms, ' Extraction: Number of parameters in each profile'
   sxaddpar, objhdr, 'XCHI2', mean(chisq), ' Extraction: Mean chi^2'
   sxaddpar, objhdr, 'DESIGNID', yanny_par(slithdr,'designid'), ' Design ID'
   sxaddpar, objhdr, 'NAXIS1', (size(flux))[1], ' X Size of flux ext'
   sxaddpar, objhdr, 'NAXIS2', (size(flux))[2], ' Y Size of flux ext'

   sxaddpar,objhdr,'EQUINOX',2000.0,after='DEC'
   sxaddpar,objhdr,'RADECSYS', 'FK5', after='EQUINOX'
   sxaddpar,objhdr,'AIRMASS', float(tai2airmass(ra, dec, tai=tai_mid)) * (ct_ra EQ 1) * (ct_dec EQ 1) * (tai_mid GT 0), after='ALT'
   sxaddpar, objhdr, 'EXTEND', 'T', after='NAXIS2'

  ; Set ivar to zero where NODATA, FULLREJECT, or LOWFLAT maskbits set
  index=where((pixelmask and sdss_flagval('SPPIXMASK','NODATA')) $
    or (pixelmask and sdss_flagval('SPPIXMASK','FULLREJECT')) $
    or (pixelmask and sdss_flagval('SPPIXMASK','LOWFLAT')) ne 0,nindex)
  if (nindex ne 0) then fluxivar[index]=0.

  ; Interpolate over masked pixels for cosmetic purposes
  flux=djs_maskinterp(flux,fluxivar EQ 0, /const, iaxis=0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Quality control
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; Quality control: add information to the DRPQUAL bitmask
  ; File is valid
  drpqual=sdss_flagval('MANGA_DRPQUALFLAG','VALIDFILE')
  ; Look for scattered light warning flags
  if (scat_warning ne 0) then drpqual=drpqual OR sdss_flagval('MANGA_DRPQUALFLAG','HIGHSCAT')
  if (scat_warning eq 2) then drpqual=drpqual OR sdss_flagval('MANGA_DRPQUALFLAG','SCATFAIL')
  ; Look for cases where less than 30% of the pixels have good value
  index=where(fluxivar ne 0,nindex)
  if (nindex/float(n_elements(fluxivar)) lt 0.3) then $
    drpqual=drpqual OR sdss_flagval('MANGA_DRPQUALFLAG','EXTRACTBAD')
  ; Flag case where extracted fluxes are huge
  ; (e.g., perhaps this is a twilight flat)
  if (nindex gt 0) then begin
    if (median(flux[index]) gt 1e4) then $
      drpqual=drpqual OR sdss_flagval('MANGA_DRPQUALFLAG','EXTRACTBRIGHT')
  endif
  ; Flag cases where the exposure time was abnormally short
  ; (less than 5 minutes)
  if (float(sxpar(objhdr,'EXPTIME')) lt 300.) then $
    drpqual=drpqual OR sdss_flagval('MANGA_DRPQUALFLAG','LOWEXPTIME')

  ; Loop over harnesses and look for entire missing IFUs
  ; (e.g., that might have fallen out)
  badifu=''
  harnames=partialslit[uniq(partialslit.harname)].harname
  nhar=n_elements(harnames)
  for i=0,nhar-1 do begin
    ; Identify the IFU fibers in this block
    thesefibers=where((partialslit.harname eq harnames[i])and(partialslit.fibertype eq 'IFU'),nthese)
    ; Sum the fluxes (silly way to do it, but it works)
    theflux=total(flux[*,thesefibers])
    if (theflux eq 0.) then begin
      splog, strcompress('WARNING! Missing all IFU fibers in harness '+harnames[i])
      badifu+=harnames[i]+' '
      drpqual=drpqual OR sdss_flagval('MANGA_DRPQUALFLAG','BADIFU')
    endif
  endfor

  ; Check the dither position information
  ; If MGDPOS doesn't match MGDRA and MGDDEC there's a problem
  mgdpos=strtrim(sxpar(objhdr,'MGDPOS'),2)
  mgdra=float(sxpar(objhdr,'MGDRA'))
  mgddec=float(sxpar(objhdr,'MGDDEC'))
  result=mdrp_checkdither(mgdpos,mgdra,mgddec,platescale=60.)
  if (result eq 1) then drpqual=drpqual OR sdss_flagval('MANGA_DRPQUALFLAG','BADDITHER')

  ; Add DRPQUAL and MISSIFU to the FITS header
  sxaddpar, objhdr, 'DRPQUAL', drpqual, 'DRP quality bitmask'
  sxaddpar, objhdr, 'BADIFU', badifu, 'Harness names of missing/bad ifus'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Write extracted, lambda-calibrated spectra to disk following the MaNGA data model
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

   ml_mwrfits,dummyext,outname,hdr=objhdr, /create ; Blank ext 0 with full header
   ml_mwrfits, flux, outname, extname='FLUX' ;extracted flux
   ml_mwrfits, fluxivar, outname, extname='IVAR'    ; extracted inverse variance
   ml_mwrfits, mangamask, outname, extname='MASK'  ;pixel mask
   ml_mwrfits, vacset, outname, extname='WSET'     ;trace-set for wavelength sol
   ml_mwrfits, dispset, outname, extname='DISPSET'    ;trace-set for dispersion sol
   ml_mwrfits, slitmap, outname, hdr=slithdr, extname='SLITMAP'   ;slitmap
   ml_mwrfits, xnow, outname, extname='XPOS'       ;x pos of traces on CCD
   ml_mwrfits, superfit, outname, extname='SUPERFLAT'   ;superflat vector from quartz lamps

   spawn, ['gzip', '-f', outname], /noshell

  if ~keyword_set(visual) then begin
    dfpsclose
    ; Convert to PDF and remove old PS file
    spawn, strcompress('ps2pdf '+plotname+' '+ml_strreplace(plotname,'.ps','.pdf'))
    spawn, strcompress('rm -f '+plotname)
    !p.multi=0
    !p=origp
    !x=origx
    !y=origy
  endif

   heap_gc
   obj_destroy,configuration

   return   
 
 end
