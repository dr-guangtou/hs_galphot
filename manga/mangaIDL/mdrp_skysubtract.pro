;+
; NAME:
;   mdrp_skysubtract
;
; PURPOSE:
;   This is the high-level sky subtraction routine, analagous to what was
;   done in BOSS extract_object.pro.  It calls mdrp_calcsky to actually
;   figure out the sky values.
;
; CALLING SEQUENCE:
;   mdrp_skysubtract, innames, [outnames=, $
;      /visual, ssmode=, /scalesky, /pimg]
;
; INPUTS:
;   innames - Array of names of the input mgFrame files
;
; OPTIONAL INPUTS:
;   /visual - Send plots to display instead of to a file
;   plotname - File for the plot (overrides default)
;   ssmode - '1d', '2d', or 'hybrid' sky subtraction mode.  Default is 'hybrid'
;   /scalesky - Scale sky fiber spectra to a mean value first
;   /pimg - Spit out poisson ratio images
;   /writesn2 - writes outs sn2 and fibermag values into a fits file
;
; OPTIONAL OUTPUTS:
;   outnames - Array of names of the output mgSFrame files
;
; PROCEDURES CALLED:
;   mgframe_read
;   mlstrreplace()
;   ml_ensurenogz()
;   ml_psfilename()
;   mdrp_calcsky()
;   djs_median()
;   qaplot_skysub
;   bspline_iterfit()
;   mlmeanclip()
;   djs_maskinterp()
;   bspline_valu()
;   sdss_flagval()
;   ml_mwrfits
;
; REVISION HISTORY:
;   v1.0: 28-Jan-2013  David R. Law (Dunlap Institute; drlaw@di.utoronto.ca)
;       Imported code from BOSS extract_object.pro, start integrating
;       with MaNGA algorithms and calling
;   v2.0: 27-Sep-2013  D. Law
;       Major overhaul for MaNGA production run
;   07-Feb-2014 Updates based on real data (D. Law)
;   19-Feb-2014 Updates for improved performance (D. Law)
;   02-May-2014 Added in SN2 calculations, QA plots (D. Law)
;   15-May-2014 Added quality control bitmask logic, tweaked sky fiber
;     selection to avoid masked fibers (D. Law)
;   30-May-2014 Retool processing for scaled sky subtraction (D. Law)
;   14-Jul-2014 Added writesn2 keyword for writing out the sn2, and fibermag, values for each exposure
;-

pro mdrp_skysubtract, innames, visual=visual, plotname=plotname, $
 outnames=outnames, ssmode=ssmode, pimg=pimg, scalesky=scalesky, writesn2=writesn2

  on_error,0
  compile_opt idl2

; How many frames are there to process?
nframes=n_elements(innames)

outnames=strarr(nframes)

; Default subtraction mode is 'hybrid'
if (~keyword_set(ssmode)) then ssmode='hybrid'
; Valid modes are '1d','2d','hybrid'
if ((ssmode ne '1d')and(ssmode ne '2d')and(ssmode ne 'hybrid')) then begin
  splog,'WARNING: ssmode '+string(ssmode)+' not recognized, defaulting to hybrid.'
  ssmode='hybrid'
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LOOP over exposure number
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

for i=0,nframes-1 do begin

; Set the output name
outnames[i]=ml_strreplace(innames[i],'mgFrame','mgSFrame')
; Make sure outname name doesn't end in '.gz'
outnames[i]=ml_ensurenogz(outnames[i])

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Read in the file, set key information
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

ml_mgframeread,innames[i],objflux=flux,objivar=fluxivar,mask=pixelmask,wset=wset,loglam=loglam, $
 dispset=dispset,dispimg=dispimg,slitmap=slitmapFULL,ximg=ximg,superflat=superflat,hdr=hdr,slithdr=slithdr

camera=strtrim(fxpar(hdr,'CAMERAS'),2)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Quality control checks.  If fails go to next frame.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Quality check- check that input file read properly
; If not, go to next frame in the loop
if (~keyword_set(hdr)) then begin
  splog,strcompress('File '+innames[i]+' not found, skipping!')
  continue
endif
; Check that no DRP status flags set that would prohibit processing
failout=0
drpqual=long(fxpar(hdr,'DRPQUAL'))
; If the VALIDFILE bit wasn't set, skip frame and go to next frame in loop
if ((drpqual and sdss_flagval('MANGA_DRPQUALFLAG','VALIDFILE')) eq 0) then failout=1
; If the EXTRACTBAD bit was set, skip frame and go to next frame in loop
if ((drpqual and sdss_flagval('MANGA_DRPQUALFLAG','EXTRACTBAD')) ne 0) then failout=1
; If the EXTRACTBRIGHT bit was set, skip frame and go to next frame in loop
if ((drpqual and sdss_flagval('MANGA_DRPQUALFLAG','EXTRACTBRIGHT')) ne 0) then failout=1
; If the SCATFAIL bit was set, skip frame (HIGHSCAT is ok, SCATFAIL is not)
if ((drpqual and sdss_flagval('MANGA_DRPQUALFLAG','SCATFAIL')) ne 0) then failout=1
; If the LOWEXPTIME bit was set, skip frame
if ((drpqual and sdss_flagval('MANGA_DRPQUALFLAG','LOWEXPTIME')) ne 0) then failout=1

if (failout) then begin
  splog,strcompress('File '+innames[i]+' has quality problems ('+sdss_flagname('MANGA_DRPQUALFLAG',drpqual,/concat)+'), skipping!')
  continue
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End quality control checks
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Reading slithdr from the mgFrame is buggy and doesn't work so well
; when trying to write back out again.  Read it from mangacore
plugname=yanny_par(slithdr,'plugmap')
slitname=ml_strreplace(plugname,'plPlugMapM','slitmap')
temp=ml_readslitmap(slitname,hdr=slithdr)

; Which spectrograph, 1 or 2?
specid=0
if ((camera eq 'b1')or(camera eq 'r1')) then specid=1
if ((camera eq 'b2')or(camera eq 'r2')) then specid=2
; Fail out if spectrograph wasn't set to either 1 or 2
if ((specid ne 1)and(specid ne 2)) then $
  message,'FAIL- unknown camera '+camera+' in sky subtraction!'
; Which color camera, b or r?
color='r'
if ((camera eq 'b1')or(camera eq 'b2')) then color='b'
if ((camera eq 'r1')or(camera eq 'r2')) then color='r'

; Cull the slitmap lines pertinent to this detector
slitmap=slitmapFULL[where(slitmapFULL.spectrographid eq specid)]

; Set up display if specified, otherwise set plot file
if (keyword_set(visual)) then begin
  set_plot,'x'
  device,decomposed=0
  loadct,0,/silent
  !p.multi=[0,1,2]
  !p.thick=0 & !x.thick=0 & !y.thick=0 & !p.charthick=0
endif else begin
  set_plot,'ps'
  origp=!p
  origx=!x
  origy=!y
  !p.multi=[0,1,2]
  !p.thick=0 & !x.thick=0 & !y.thick=0 & !p.charthick=0
  if (keyword_set(plotname)) then plotname=plotname $
  else begin ; Set default plotname
    plotname=fileandpath(ml_psfilename(outnames[i]),path=path)
    path=path+'qa/'
    if file_test(path,/directory) eq 0 then spawn, '\mkdir -p '+path
    plotname=path+plotname
  endelse
  splog,'Printing to ',plotname
  dfpsplot, plotname, /color
endelse

nx = (size(flux,/dim))[0] 
ny = (size(flux,/dim))[1] 

; Set up the blue/red frame parameters
if (color eq 'b') then nbkpt = 3*nx/4
if (color eq 'r') then nbkpt = nx

; Find which fibers are masked out in the input data
; (e.g., if they fell out after mapping)
fibermask=bytarr(ny)
for k=0,ny-1 do begin
  index=where((pixelmask[*,k] and sdss_flagval('SPPIXMASK','NODATA')) $
    or (pixelmask[*,k] and sdss_flagval('SPPIXMASK','FULLREJECT')) $
    or (pixelmask[*,k] and sdss_flagval('SPPIXMASK','LOWFLAT')) ne 0,nindex)
  ; If over 50% of the fiber pixels are bad, the fiber is not  a good sky fiber
  ; DRL changed this to 50% on June 2 2014- previous 90% wasn't
  ; catching some major failures with 72% fail points
  if (nindex/float(nx) gt 0.5) then fibermask[k]=1
endfor

; Define where the sky fibers are using the slitmap
; and fibermask array
iskies=where((slitmap.fibertype eq 'SKY') AND (slitmap.plugstatus eq 0) AND (fibermask eq 0), nskies)

; Define observation times
tai_start=fxpar(hdr,'TAI-BEG'); MDJ(TAI) in seconds at exposure start
exptime=fxpar(hdr,'EXPTIME')
tai_mid=tai_start+exptime/2.

  ; Read in the skyline flagging file
  skyfile=concat_dir(ml_getenv('MANGADRP_DIR'),'etc/skylines.par')
  skylines=yanny_readone(skyfile,'SKYLINE')
  nlines=n_elements(skylines)
  if (nlines eq 0) then $
     message,'WARNING: skyline file cannot be found or is empty!'

  ; Make a reference mask where the skylines are for this input frame
  skypad=3.; Pad around input values by 3 Angstroms
  skylinemask=intarr(nx,ny)
  ; DRL- this loop takes a long time, how can it be improved?
  splog,'Computing skyline mask'
  for j=0,nlines-1 do begin
    index=where((10.^loglam gt skylines[j].skywave-skypad)and(10.^loglam lt skylines[j].skywave+skypad),nindex)
    if (nindex gt 0) then skylinemask[index]=1
  endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; First call to mdrp_calcsky
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

skystruct = mdrp_calcsky(flux, fluxivar, loglam, slitmap, $
  skysub, skysubivar, iskies=iskies, pixelmask=pixelmask, skylinemask=skylinemask, $
  fibermask=fibermask, upper=3.0, lower=3.0, tai=tai_mid, nbkpt=nbkpt, scalesky=scalesky)

if (NOT keyword_set(skystruct)) then $
  message,'Problem with skysubtract- quit!'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Second call to mdrp_calcsky
; If any of the sky-fibers are bad, then re-do sky-subtraction.
; Make quality plot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Bad fibers are where median sky subtracted flux > sqrt(2) * sigma
ibadfib = where(djs_median(skysub[*,iskies]^2 * $
  skysubivar[*,iskies], 1) GT 2.0)               
if (ibadfib[0] NE -1) then begin
  fibermask[iskies[ibadfib]] = fibermask[iskies[ibadfib]] OR fibermask_bits('BADSKYFIBER')

  splog, 'Calling skysubtract again; masked skyfibers',string(iskies[ibadfib])
  skystruct = mdrp_calcsky(flux, fluxivar, loglam, slitmap, $
    skysub, skysubivar, iskies=iskies, pixelmask=pixelmask, skylinemask=skylinemask, $
     fibermask=fibermask, upper=10.0, lower=10.0, tai=tai_mid, nbkpt=nbkpt, scalesky=scalesky)

  if (NOT keyword_set(skystruct)) then $
    message,'Problem with skysubtract- quit!'
endif

; QA plots for chi^2 from 1D sky-subtraction.
;qaplot_skysub, flux, fluxivar, skysub, skysubivar, vacset, iskies, title=' 1D Sky-subtraction'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Third call to mdrp_calcsky
; This time incorporate 2d information along the slit.
; Only call if ssmode='2d' or 'hybrid'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if ((ssmode eq '2d')or(ssmode eq 'hybrid')) then begin
  ; Two d fit along the slit direction
  ; DRL: Found on 7-Feb-2014 that using slit order for the sort dimension
  ; gives *slightly* better results than using dispersion values
  twodaxis=replicate(1.0,nx) # findgen(ny)

  ; DRL 7-Feb-2014: 3rd order better than 2nd, but didn't see any
  ; gains going to 4th order
  nskypoly = 3L
  skystruct2d = mdrp_calcsky(flux, fluxivar, loglam, slitmap, $
    skysub2d, skysubivar2d, iskies=iskies, pixelmask=pixelmask, $
    fibermask=fibermask, upper=10.0, lower=10.0, tai=tai_mid, $
    twodaxis=twodaxis, skylinemask=skylinemask, $
    npoly=nskypoly, nbkpt=nbkpt, $
    relchi2set=relchi2set, newmask=newmask, scalesky=scalesky)
 ; DRL Feb 2014: Don't use newmask until we've gone through it in depth
 ; pixelmask = newmask

  if (NOT keyword_set(skystruct2d)) then $
    message,'Problem with 2d skysubtract- quit!'

; QA plots for chi^2 from 2D sky-subtraction.
;qaplot_skysub, flux, fluxivar, skysub, skysubivar, $
 ; vacset, iskies, title=' 2D Sky-subtraction'

endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; If ssmode='hybrid' was selected, use the 2d fit in the regions of
; skylines and the ordinary fit elsewhere.  This optimizes the sky
; model quality in the continuum.
; Otherwise use either the 1d or 2d fit if ssmode='1d' or '2d' respectively
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if (ssmode eq '1d') then begin
  ; Nothing to do
  splog,'Using 1d sky fitting only'
endif

; If ssmode is '2d' then assign the 2d results to skysub
if (ssmode eq '2d') then begin
  splog,'Using 2d sky fitting'
  skysub=skysub2d
  skysubivar=skysubivar2d
endif

; If ssmode is 'hybrid' then combine fits
if (ssmode eq 'hybrid') then begin
  splog,'Using hybrid 1d+2d sky fitting'

  ; DRL- Feb 19 2014.  Hack for low fiberid on r2, which
  ; looks better without twod fit on some lines??
  ; How general is this???
  if (camera eq 'r2') then begin
    for j=0,141 do begin
      index=where((10.^loglam[*,j] gt 5934.)and(10.^loglam[*,j] lt 9238.),nindex)
      if (nindex gt 0) then skylinemask[index,j]=0
    endfor
  endif

  index=where(skylinemask eq 1,nindex)
  if (nindex gt 0) then begin
    skysub[index]=skysub2d[index]
    skysubivar[index]=skysubivar2d[index]
  endif
endif


; Save the sky-subtracted flux values as is, and now modify flambda.
flambda = skysub ; Sky-subtracted flux
flambdaivar = skysubivar ; Sky-subtracted ivar
skyimg = flux - flambda ; Sky flux

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Copy input to output header
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Update version information in FITS headers if needed
objhdr=mdrp_makeheader(head=hdr,/drp2d)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Calculate SN2 based on PlateMags file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

sn2 = mdrp_calcsn2(flux, flambda, flambdaivar, loglam, slitmap, hdr, slithdr, color=color, totfiber=ny, exptime=float(exptime))
; Put the info in the output header
sxaddpar,objhdr,'FRAMESN2', sn2[0].sn2, "(S/N)^2 at fidicial magnitude"

if keyword_set(writesn2) then begin
  sn2 = jjadd_tag(sn2, 'name', camera+'-'+(strsplit(outnames[i],'-.',/extract))[2])
  if i eq 0 then begin 
    ml_mwrfits, 0, 'sn2vals.fits', /create
    ml_mwrfits, sn2,'sn2vals.fits', extname=(camera+'-'+(strsplit(outnames[i],'-.',/extract))[2])
  endif else ml_mwrfits, sn2,'sn2vals.fits', extname=(camera+'-'+(strsplit(outnames[i],'-.',/extract))[2])  
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Clean up results
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;----------
; Combine FIBERMASK and PIXELMASK to FINALMASK
finalmask = pixelmask
ntrace=ny ;???
for itrace=0, ntrace-1 do $
  finalmask[*,itrace] = finalmask[*,itrace] OR fibermask[itrace]

;----------
; Disable some mask bits in regions where 'NODATA' is set
q_nodata = (finalmask AND sdss_flagval('SPPIXMASK','NODATA')) NE 0
discards = ['NEARBADPIXEL','LOWFLAT','SCATTEREDLIGHT','NOSKY']
for j=0, n_elements(discards)-1 do $
  finalmask = finalmask - q_nodata $
  * (finalmask AND sdss_flagval('SPPIXMASK',discards[j]))

;----------
; Get an estimate of the relative chi^2 at each pixel.
; Do this with a simple linear interpolation.
if (keyword_set(relchi2set)) then begin
  xx = 0
;   rchi2img = interpol(relchi2struct.chi2, relchi2struct.wave, loglam)
  rchi2img = bspline_valu(loglam, relchi2set)
   ; Compute the mean relative chi2 of sky-subtraction, after masking
   ; bad regions of the CCD
  fval = sdss_flagval('SPPIXMASK','NOPLUG') $
    + sdss_flagval('SPPIXMASK','BADTRACE') $
    + sdss_flagval('SPPIXMASK','BADFLAT') $
    + sdss_flagval('SPPIXMASK','BADARC') $
    + sdss_flagval('SPPIXMASK','LOWFLAT') $
    + sdss_flagval('SPPIXMASK','NOSKY') $
    + sdss_flagval('SPPIXMASK','NODATA') $
    + sdss_flagval('SPPIXMASK','BADFLUXFACTOR')
  indx = where((finalmask AND fval) EQ 0 AND flambdaivar NE 0, ct)
  if (ct EQ 0) then skychi2 = 0. $
  else skychi2 = mean(rchi2img[indx])
endif else begin
  rchi2img = 0 * flambda + 1.
  skychi2 = 0.
endelse
sxaddpar, objhdr, 'SKYCHI2', skychi2, ' Mean chi^2 of sky-subtraction'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Cosmetic Cleanup (DRL Feb 2014):
; Wherever flambdaivar=0 set flambda to 0.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

index=where(flambdaivar eq 0.,nindex)
if (nindex gt 0) then flambda[index]=0.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Poisson metric calculations
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Calculate official pstat
pstat=-1.; Default value

; Calculate the Poisson ratio image for the entire frame
pratio=flambda*superflat/sqrt(skyimg*superflat)
; Give zero value to anything with ivar=0
index=where(flambdaivar eq 0,nindex)
if (nindex ne 0) then pratio[index]=0.
; If /pimg keyword set, write it out
if (keyword_set(pimg)) then begin
  pfile=ml_strreplace(outnames[i],'Frame','pimg')
  writefits,pfile,pratio
endif

wave=ml_setwcalib(); Global wavelength solution
nxnew=n_elements(wave)
;Wavelength at which to calculate Pstat
; (chosen to be a bright skyline)
if ((camera eq 'r1')or(camera eq 'r2')) then lampeak=9378.
if ((camera eq 'b1')or(camera eq 'b2')) then lampeak=5462.
; Identify x location closest to this wavelength on the global solution
temp=abs(wave-lampeak)
ipeak=where(temp eq min(temp))
ipeak=ipeak[0]
splog,'Calculating Poisson ratio statistic at pixel '+strcompress(string(ipeak),/remove_all)

; Identify only plugged sky fibers
index=where((strcompress(slitmap.fibertype,/remove_all) eq 'SKY')and(slitmap.plugstatus eq 0),nindex)
; Crop the Poisson ratio image if sky fibers successfully identified
if (nindex gt 0) then begin
  pratio=pratio[*,index]
  lam=10.^loglam[*,index]
  pivar=flambdaivar[*,index]
  ; Wave rectify.  Propagate ivar as well since we'll use it
  ; to reject bad pixels in the calculation.
  prationew=fltarr(nxnew,nindex)
  pivarnew=fltarr(nxnew,nindex)
  for j=0,nindex-1 do begin
    prationew[*,j]=interpol(pratio[*,j],lam[*,j],wave)
    pivarnew[*,j]=interpol(pivar[*,j],lam[*,j],wave)
  endfor
  ; If /pimg keyword set, write it out
  if (keyword_set(pimg)) then begin
    pfile=ml_strreplace(outnames[i],'Frame','pimgsky')
    writefits,pfile,prationew
  endif
  ; Crop out subregion around the line of interest
  subregion=prationew[ipeak-2:ipeak+2,*]
  subivar=pivarnew[ipeak-2:ipeak+2,*]
  ; Do the computation on pixels that aren't bad
  index2=where(subivar ne 0,nindex2)
  if (nindex2 gt 0) then begin
    temp=ml_meanclip(subregion[index2],mean,rms,clipsig=5.)
    pstat=rms
  endif

endif

splog,'P='+strcompress(string(pstat),/remove_all)

; Add Poisson ratio statistic to FITS header
sxaddpar, objhdr, 'PSTAT',pstat,' Poisson ratio statistic'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Quality control
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ; Quality control: add information to the DRPQUAL bitmask
  drpqual=long(fxpar(objhdr,'DRPQUAL'))
  ; Flag cases where the Pstat was greater than 1e5
  ; (example case)
  if (pstat gt 1e5) then $
    drpqual=drpqual OR sdss_flagval('MANGA_DRPQUALFLAG','SKYSUBBAD')
  ; Add DRPQUAL to the FITS header
  sxaddpar, objhdr, 'DRPQUAL', drpqual, 'DRP quality bitmask'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Write out results
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Write sky-subtracted spectra to disk following the MaNGA data model
ml_mwrfits, dummyext, outnames[i],hdr=objhdr, /create ; Blank ext 0 with full header
ml_mwrfits, flambda, outnames[i], extname='FLUX' ;sky subtracted flux
ml_mwrfits, flambdaivar, outnames[i], extname='IVAR'   ; sky subtracted inverse variance
ml_mwrfits, finalmask, outnames[i], extname='MASK' ; final pixel mask
ml_mwrfits, wset, outnames[i], extname='WSET'    ;trace-set for wavelength sol; wset
ml_mwrfits, dispset, outnames[i], extname='DISPSET'   ;trace-set for dispersion sol
ml_mwrfits, slitmapFULL, outnames[i], hdr=slithdr, extname='SLITMAP'  ;slitmap used
ml_mwrfits, ximg, outnames[i], extname='XPOS' ;x pos on CCD
ml_mwrfits, superflat, outnames[i], extname='SUPERFLAT'  ;superflat vector from quartz lamps
ml_mwrfits, skyimg, outnames[i], extname='SKY' ;sky flux

; gzip output file
spawn, ['gzip', '-f', outnames[i]], /noshell

; Close out plots
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; End LOOP over exposure number
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

endfor

; Take out the trash
heap_gc

return
end
