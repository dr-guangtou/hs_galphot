;+
; NAME:
;   mdrp_reduceoneifu
;
; PURPOSE:
;   Reduce a single IFU through the 3d MaNGA pipeline.
;   This program is called once by mdrp_reduce3d for each IFU.
;
; CALLING SEQUENCE:
;   mdrp_reduceoneifu
;
; INPUTS:
;   framenames: String array containing full filepaths to the input frames
;   plan: mgPlan3d structure describing the exposures to combine
;   ifuDesign: The one ifu to combine data for (e.g., 12701).
;
; OPTIONAL INPUTS:
;   outdir: Output directory, by default $MANGA_SPECTRO_REDUX/DRPVER/PLATE/stack/
;   outnameroot: Root output name, by default manga-$plateid-$ifudesign-
;         Program will scan this string looking for '$' symbols.  If none found,
;         it simply uses this string as the root.  If found, it replaces these
;         elements with proper variable values.  Requires variables to be '-' delimited.
;   cubescale: Output spaxel size in cube in arcsec (default 0.5 arcsec/spaxel)
;   extast_log: Output logfile for the extended astrometry module
;
; OPTIONAL KEYWORDS:
;   /astweak: Set this to turn on extended astrometry module
;     (broadband registration)
;
; OUTPUT:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   This routine spawns the Unix command 'mkdir'.
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
; REVISION HISTORY:
;   October-2013  Written by David Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   19-Mar-2014 Updated for production run (D. Law)
;   02-Oct-2014 Make default spaxel scale 0.5 arcsec (D. Law)
;-
;------------------------------------------------------------------------------

pro mdrp_reduceoneifu, framenames, plan, ifuDesign, outdir=outdir, outnameroot=outnameroot1, astweak=astweak, logonly=logonly, linonly=linonly, dominimal=dominimal, cubescale=cubescale, extast_log=extast_log

; Default error behavior is to stop immediately
on_error,0
; Use 32-bit integers by default, and strict array usage
compile_opt idl2

; Initialize drp3qual to 0
drp3qual=0L

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Setup- define number of exposures, output directory, and
; partial slitmap for the first exposure (i.e., slitmap cropped
; to only the ifuDesign value)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

nframes=n_elements(framenames)

; Default output spaxel size is 0.5 arcsec
if (keyword_set(cubescale)) then cubescale=cubescale $
else cubescale=0.5

; Check that all frames exist
if (total(file_test(framenames)) ne nframes) then $
  message,'Frames not found!'

; Check that the number of exposures matches the plan file
if (nframes ne n_elements(plan)) then $
  message, 'Error: Number of exposures doesnt match plan file!'

splog,'nframes = '+strcompress(string(nframes),/remove_all)

; Check that ifuDesign has either a single value or 1 value for each exposure
ndesign=n_elements(ifuDesign)
if ((ndesign ne 1)and(ndesign ne nframes)) then $
  message, 'ERROR: Too many or too few values of ifuDesign'

; Read in slitmap for first exposure
ml_mgframeread,framenames[0],slitmap=slitmap
; partialslit is the slitmap entry for the 1st exposure for this IFU only
partialslit=slitmap[where(slitmap.ifuDesign eq ifuDesign[0])]

; Set output directory using provided filepath, otherwise construct it
; based on MaNGA data model
if (keyword_set(outdir)) then outdir=outdir $
else begin
  outdir=concat_dir(ml_getenv('MANGA_SPECTRO_REDUX'),mangadrp_version(/simple))
  outdir=concat_dir(outdir,strcompress(string(plan[0].plate),/remove_all))
  outdir=concat_dir(outdir,'stack/')
endelse

splog,'Output directory set to: '+outdir
; If output directory doesn't exist, create it and any necessary parent directories
spawn, 'mkdir -p ' + outdir

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Setup- define output filename root
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
; Set the output filename root (manga-$plateid-$ifudesign-)
if (keyword_set(outnameroot1)) then outnameroot=outnameroot1 $
else outnameroot='manga-$plateid-$ifudesign-'

; If outnameroot has any variable names in it (indicated by $) replace them
; now with actual values.  NOTE!  This needs to be explicitly done for
; each possible variable.  Requires variables to be '-' delimited.
pieces=strsplit(outnameroot,'-',/extract)
npiece=n_elements(pieces)
for i=0,npiece-1 do begin
  ; If this piece starts with a $ symbol, remove the symbol to get variable name
  if (strpos(pieces[i],'$') eq 0) then begin
    varname=strmid(pieces[i],1)
    ; Depending on varname, change to the correct variable value
    if (varname eq 'plateid') then pieces[i]=strcompress(string(plan[0].plate),/remove_all)
    if (varname eq 'ifudesign') then pieces[i]=strcompress(string(ifuDesign[0]),/remove_all)
    if (varname eq 'mangaid') then pieces[i]=strcompress(string(partialslit[0].mangaID),/remove_all)
  endif
endfor
; Reassemble the pieces into the final output filename root
outnameroot=''
for i=0,npiece-1 do outnameroot+=pieces[i]+'-'

splog,'Outname root set to: '+outnameroot

; If not provided with an extended astrometry filename, set it here
if (~keyword_set(extast_log)) then extast_log=outdir+'tweakastrometry.out'

; Set directory for qa plots
qaplotdir=concat_dir(outdir,'qa/')
; Create directory if needed
spawn, 'mkdir -p ' + qaplotdir
; Define qa filename
qaplotfile=qaplotdir+outnameroot+'qa.ps'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Setup- define total number of fibers, preimaging directory
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Determine total number of fibers across all exposures by
; reading in each of the slitmaps.  In principle, there's no
; reason why this can't handle different physical ifus, with a
; different number of fibers, for each exposure.
nftotal=0
; Loop over exposures
for i=0,nframes-1 do begin
  ; Allow for multiple values of ifuDesign
  if (ndesign eq 1) then thisdesign=ifuDesign $
  else thisdesign=ifuDesign[i]

  ; Read in the slitmap for this exposure
  ml_mgframeread,framenames[i],slitmap=slitmap
  ; Crop out lines appropriate for this IFU
  partialslit=slitmap[where(slitmap.ifuDesign eq thisdesign)]
  nfiber=n_elements(partialslit)

  splog, 'Number of fibers in exposure '+strcompress(string(i+1),/remove_all)+' is '+strcompress(string(nfiber),/remove_all)

  if nfiber eq 0 then $
    message, 'ERROR- ifuDesign= '+strcompress(string(ifuDesign),/remove_all)+' not found in '+framenames[i]

  nftotal+=nfiber
endfor

splog,strcompress(string(nftotal),/remove_all)+' fibers total.'

; Define the preimaging directory.
preim_dir=ml_getenv('MANGA_PREIMAGING')
; If it's not defined, then turn off astrometric tweaking
if (preim_dir eq '') then begin
  splog,'MANGA_PREIMAGING not defined, turning off astrometric tweaking.'
endif

; Define the blank observing info structure.  This is the .par structure
; that consolidates key information about what went into each cube
; in a manner that can be read by the data simulator.
; We'll add one line to this structure later for each exposure
; combined.
obsinfo_line=ml_createobsstruc()

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Now we've got multiple sets of data to process.
; We have nframes worth of data for 3 data flavors:
; LINEAR: Flux calibrated data on a linear wave grid
; LOG: Flux calibrated data on a log wave grid
; MINIMAL: Flux calibrated data per camera on original wave grids
;
; Do the LOG loop first.  We'll work out astrometry here and
; apply it to the LINEAR and MINIMAL loops later.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Data has the MaNGA log wavelength solution
wave=ml_setwcalib()
nwave=n_elements(wave)
; Figure out indices where the wavelength solution best matches
; SDSS u,g,r,i,z broadband effective lambda
bbwave=[3551.,4686.,6166.,7480.,8932.]
bbindex=intarr(n_elements(bbwave))
for j=0,n_elements(bbwave)-1 do begin
  temp=min(abs(wave-bbwave[j]),minat)
  bbindex[j]=minat
endfor

; Set up arrays for the data
flux_arr=fltarr(nwave,nftotal)
ivar_arr=fltarr(nwave,nftotal)
mask_arr=lonarr(nwave,nftotal)
disp_arr=fltarr(nwave,nftotal)
x_arr=fltarr(nwave,nftotal)
y_arr=fltarr(nwave,nftotal)

; Set up arrays for the psf reference
psf_flux=fltarr(n_elements(bbwave),nftotal)
psf_xarr=fltarr(n_elements(bbwave),nftotal)
psf_yarr=fltarr(n_elements(bbwave),nftotal)

; Running index for which fiber in nftotal a frame starts on
fiberstart=0L

; Loop over individual exposures
for i=0,nframes-1 do begin
  splog,''
  splog,strcompress('Processing exposure number '+string(i+1))
  splog,'File = '+framenames[i]

  ; Allow for multiple values of ifuDesign
  if (ndesign eq 1) then thisdesign=ifuDesign $
  else thisdesign=ifuDesign[i]

  ; Read slitmap from file first
  ml_mgframeread,framenames[i],slitmap=slitmap,slithdr=slithdr

  ; Reading slithdr from the mgFrame is buggy and doesn't work so well
  ; when trying to write back out again.  Read it from mangacore
  plugname=yanny_par(slithdr,'plugmap')
  slitname=ml_strreplace(plugname,'plPlugMapM','slitmap')
  temp=ml_readslitmap(slitname,hdr=slithdr)

  ; index is the lines from slitmap appropriate for this IFU
  index=where(slitmap.ifuDesign eq thisdesign)

  ; Read data from file for only these ifu lines only
  ml_mgframeread,framenames[i],index,objflux=objflux,objivar=objivar, mask=objmask, $
    dispimg=dispimg, slitmap=partialslit, hdr=hdr

  ; If this is the first frame, save the FITS header
  if (i eq 0) then hdrsave=hdr

  ; Number of fibers read in from this exposure
  nfiber=n_elements(partialslit)

  ; Check that loglam (the wavelength solution from the data file)
  ; has nwave elements
  if ((size(objflux))[1] ne nwave) then $
    message,'ERROR- Incorrect number of wavelength samples.'
  ; DRL- check that it matches default wave solution too???
  ; Should this be generalized?

  ; Figure out the simplistic on-sky PSF fwhm for this IFU at u,g,r,i,z bands
  ; taking into account focal offset as a function of wavelength and plate location
  ; First get the basic psf from the guider images
  ; Directory for the cogimg guider images
  temp=fileandpath(framenames[i],path=cogdir)
  psfpar=psfparams(long(fxpar(hdr,'PLATEID')),long(fxpar(hdr,'MJD')),long(fxpar(hdr,'EXPOSURE')),reduxdir=cogdir)
  ; Compute the effective fwhm at these wavelengths, for this plate position
  fcpsf1d,psfpar.params,partialslit[0].xfocal,partialslit[0].yfocal,psfplane,wavearr=bbwave,xarr=xpsfplane,fwhmout=fwhm_eff,/dofwhm

  ; Populate this line of obsinfo with basic relevant information
  ml_addobsinfo,obsinfo_line,partialslit,slithdr,hdr,pf_fwhm=fwhm_eff
  ; Add this line to the obsinfo file
  if (~keyword_set(obsinfo)) then obsinfo=obsinfo_line $
  else obsinfo=[obsinfo,obsinfo_line]
  obsinfo[i].set=plan[i].set

  ; Define a new mask that will combine previous SPPIXMASK values
  ; into a MANGA_DRPPIXFLAG bitmask
  newmask=replicate(0L,size(objmask,/dimen))

  ; DRL- look for where the SPPIXMASK bits 'FULLREJECT', 'COMBINEREJ',
  ; 'LOWFLAT', or 'NODATA' are set.  Set these pixels to zero inverse variance
  ; and activate the BADDATA bit in MANGA_DRPPIXFLAG
  index=where((objmask and sdss_flagval('SPPIXMASK','COMBINEREJ')) ne 0,nindex)
  if (nindex ne 0) then begin
    objivar[index]=0.
    newmask[index]=sdss_flagval('MANGA_DRPPIXFLAG','BADDATA')
  endif
  index=where((objmask and sdss_flagval('SPPIXMASK','FULLREJECT')) ne 0,nindex)
  if (nindex ne 0) then begin
    objivar[index]=0.
    newmask[index]=sdss_flagval('MANGA_DRPPIXFLAG','BADDATA')
  endif
  index=where((objmask and sdss_flagval('SPPIXMASK','LOWFLAT')) ne 0,nindex)
  if (nindex ne 0) then begin
    objivar[index]=0.
    newmask[index]=sdss_flagval('MANGA_DRPPIXFLAG','BADDATA')
  endif
  index=where((objmask and sdss_flagval('SPPIXMASK','NODATA')) ne 0,nindex)
  if (nindex ne 0) then begin
    objivar[index]=0.
    newmask[index]=sdss_flagval('MANGA_DRPPIXFLAG','BADDATA')
  endif

  ; Figure out where the broken fibers were using the slitmap.gbu flag
  ; and add this information to the MANGA_DRPPIXFLAG bitmask
  index=where(partialslit.gbu ne 0,nindex)
  if (nindex ne 0) then begin
    tempmask=replicate(0L,nfiber)
    tempmask[index]=sdss_flagval('MANGA_DRPPIXFLAG','DEADFIBER')
    newmask = newmask OR rebin(reform(tempmask,1,nfiber),nwave,nfiber)
  endif

  ; Look for where the inverse variance is zero, and zero out the
  ; flux at these locations too so that they don't confuse anyone
  index=where(objivar eq 0.,nindex)
  if (nindex gt 0) then objflux[index]=0.

  ; Stuff the flux, inverse variance, and new bitmask into the big arrays
  flux_arr[*,fiberstart:fiberstart+nfiber-1]=objflux[*,*]
  ivar_arr[*,fiberstart:fiberstart+nfiber-1]=objivar[*,*]
  mask_arr[*,fiberstart:fiberstart+nfiber-1]=newmask[*,*]
  disp_arr[*,fiberstart:fiberstart+nfiber-1]=dispimg[*,*]

  ; Work out the basic astrometry for these fiber positions and stick them in
  ; the vectors xfiber and yfiber on-sky.  If this is the first exposure,
  ; use the /setbase keyword to mark this as the zeropoint against which
  ; DC shifts should be measured
  if (i eq 0) then begin
    mdrp_astrometry,xfiber,yfiber,partialslit,slithdr,hdr,wave
  endif else begin
    mdrp_astrometry,xfiber,yfiber,partialslit,slithdr,hdr,wave
  endelse

  ; Run the extended astrometry module to tweak up positions, so long
  ; as the preimaging directory is defined
  if ((keyword_set(astweak))and(keyword_set(preim_dir))) then begin
    ; Convert from x,y fiber positions in arcsec into RA,DEC
    rafiber=partialslit[0].ra-xfiber/3600.D/cos(partialslit[0].dec*!DPI/180.)
    decfiber=partialslit[0].dec+yfiber/3600.D

    ; Construct filepath to the reference image
    refimagepath=concat_dir(preim_dir,ml_plategroup(fxpar(hdr,'DESIGNID'),/design))+'/'+strcompress(string(fxpar(hdr,'DESIGNID')),/remove_all)+'/'
    refimagefile=lookforgzip(refimagepath+'preimage-'+strcompress(string(partialslit[0].mangaID),/remove_all)+'.fits')

    ; Can we find this file?  If so, do the astrometric tweaking
    if (file_test(refimagefile)) then begin
      ; Set up the output qa plot path
      frametemp=fileandpath(framenames[i],path=pathtemp)
      qadir=concat_dir(pathtemp,'qa/')
      ; Create qa directory if it doesn't exist
      if file_test(qadir,/directory) eq 0 then spawn, '\mkdir -p '+qadir
      ; Construct full path to output qa file
      outplot=concat_dir(qadir,(strsplit(frametemp, '.fits', /regex, /extract))[0]+'_'+ifuDesign+'_astweak.pdf')

      ; Run the module.
      ; rafiber, decfiber, objflux, objivar, objmask all have dimension [nwave, nfiber]
      ; wave is a vector of size [nwave]

       mdrp_tweakastrometry, rafiber, decfiber, objflux, objivar, objmask, wave, obsinfo[i].pf_fwhm[1:4], refimagefile, raout=raout, decout=decout, fluxout=fluxout, errorout=errorout, plan=plan[i], ifudesign=thisdesign, outputfile=extast_log, /plot, outplot=outplot, outpar=outpar;, /norot;, /debug

      ; Save derived values for dRA, dDEC, dTHETA (mean combined across bands)
      obsinfo[i].EAMfit=outpar[0:2]
      obsinfo[i].EAMfiterr=outpar[3:5]

      ; Convert from RA,DEC back to x,y fiber offsets in arcsec relative to fiducial pointing
      xfiber=-(raout-partialslit[0].ra)*3600.D*cos(partialslit[0].dec*!DPI/180.)
      yfiber=(decout-partialslit[0].dec)*3600.D
    endif
  endif

  ; Stuff the x,y fiber relative position information into the big arrays
  x_arr[*,fiberstart:fiberstart+nfiber-1]=xfiber[*,*]
  y_arr[*,fiberstart:fiberstart+nfiber-1]=yfiber[*,*]

  ; Work out the psf information.  Using the actual xfiber,yfiber info,
  ; figure out the flux that would have gone into each from a point
  ; source convolved with true on-sky profile (pre-fiber).
  ; Do this for 5 wavelengths, u,g,r,i,z
  for j=0,n_elements(bbindex)-1 do begin
    psf_xarr[j,fiberstart:fiberstart+nfiber-1]=xfiber[bbindex[j],*]
    psf_yarr[j,fiberstart:fiberstart+nfiber-1]=yfiber[bbindex[j],*]

    roff=reform(sqrt(xfiber[bbindex[j],*]^2+yfiber[bbindex[j],*]^2),nfiber)
    ftemp=replicate(0.,nfiber)
    ; Only calculate flux in fibers within 6'' of center
    index=where(roff lt 6.,nindex)
    for k=0,nindex-1 do begin
      junk=min(abs(roff[index[k]]-xpsfplane),minat)
      ftemp[index[k]]=psfplane[minat,j]
    endfor
    psf_flux[j,fiberstart:fiberstart+nfiber-1]=ftemp[*]
  endfor

  ; Increment the fiber running index
  fiberstart+=nfiber
endfor


; DRL- loop over sets to calculate Omega


; Quality control plots based on the obsinfo structure
ml_qaobsinfo,obsinfo,filename=qaplotfile

; DRL- add a second loop through EAM here



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Create data cube in log wavelength grid
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Cube size in arcsec based on relative fiber positions
; across all exposures
cube_zsize=nwave
; Require square format output
cube_xysize=fix(max(x_arr)-min(x_arr)) > fix(max(y_arr)-min(y_arr))

; For consistency make Odd integers true
; Convert to output pixel scale and add some padding
cube_xysize = fix(cube_xysize/cubescale) +10
; For consistency make the number of pixels even
; (Note that odd integers count as TRUE)
if (cube_xysize) then cube_xysize+=1
; Assign to both x and y
cube_xsize = cube_xysize
cube_ysize = cube_xysize

; Create data cube and inverse variance structure
cube_flux=fltarr(cube_xsize,cube_ysize,cube_zsize)
cube_ivar=fltarr(cube_xsize,cube_ysize,cube_zsize)
cube_mask=lonarr(cube_xsize,cube_ysize,cube_zsize)

; Define 1d vectors of x,y,flux,ivar,mask to be used
; at each wavelength slice
xvec=fltarr(nftotal)
yvec=fltarr(nftotal)
fvec=fltarr(nftotal)
ivec=fltarr(nftotal)
mvec=lonarr(nftotal)

; Define output psf cube
psfcube_xsize=cube_xsize*1.
psfcube_ysize=cube_ysize*1.
psfcube_scale=cubescale/1.
psfcube=fltarr(psfcube_xsize,psfcube_ysize,n_elements(bbwave))
; FWHM of single-gaussian fit to resulting cube
cubefwhm=fltarr(n_elements(bbwave))

; Figure out reconstructed psf in each of the 5 wavebands
for i=0,n_elements(bbwave)-1 do begin
  xvec=psf_xarr[i,*]/psfcube_scale+psfcube_xsize/2.
  yvec=psf_yarr[i,*]/psfcube_scale+psfcube_ysize/2.
  fvec=psf_flux[i,*]
  scale=psfcube_scale*psfcube_scale/!PI
  psfcube[*,*,i]=ml_griddata(xvec,yvec,fvec,[psfcube_xsize,psfcube_ysize],1.6/psfcube_scale,0.7/psfcube_scale,scale=scale)
  temp=gauss2dfit(psfcube[*,*,i],coef)
  cubefwhm[i]=mean(coef[2:3])*2.355*psfcube_scale
endfor

; Loop over all wavelengths figuring out the cube at that wavelength
for i=0,cube_zsize-1 do begin
  ; Print a message to splog every 10% of the loop
  if (i mod fix(cube_zsize/10) eq 0) then $
    splog,'Constructing cube: '+strcompress(string(round(i*100./cube_zsize)),/remove_all)+'% complete'

  ; Create input vectors for this wavelength slice
  xvec[*]=x_arr[i,*]/cubescale+cube_xsize/2.
  yvec[*]=y_arr[i,*]/cubescale+cube_ysize/2.
  fvec[*]=flux_arr[i,*]
  ivec[*]=ivar_arr[i,*]
  mvec[*]=mask_arr[i,*]

  ; Scale correction factor is the ratio of area between a PI area fiber
  ; (in arcsec^2) and the output pixel size in arcsec^2
  ; The result means that the cube will be in calibrated units/pixel
  scale=cubescale*cubescale/!PI

;if (i eq 2545) then stop
;if (i eq 1312) then stop

  cube_flux[*,*,i]=ml_griddata(xvec,yvec,fvec,[cube_xsize,cube_ysize],1.6/cubescale,0.7/cubescale,ivar=ivec,invarimg=invarimg,maskvec=mvec,maskimg=maskimg,scale=scale)
  ; Populate inverse variance and mask arrays
  cube_ivar[*,*,i]=invarimg
  cube_mask[*,*,i]=maskimg
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Work out the median spectral resolving power
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Median across all fibers in this IFU
; Max variability in 12705 is about XX
disp_med=median(disp_arr,dimension=2)

; Figure out variability at each lambda between
; fibers and between exposures (i.e., variability
; along the slit, and from exp to exposure)
specvar=fltarr(nwave)
for i=0,nwave-1 do begin
  ; Sigma clip to take out lone bad values
  temp=ml_meanclip(disp_arr[i,*],themean,thesig,ignore=0.)
  ; Fractional RMS variability of the LSF
  specvar[i]=thesig/themean
endfor

; QA flag if the LSF varies by more than 5% across 5% or more of the
; wavelength range
temp=where(specvar gt 0.05,nhighvar)
if (float(nhighvar)/nwave gt 0.05) then $
  drp3qual=drp3qual OR sdss_flagval('MANGA_DRP3QUAL','VARIABLELSF')

; Take a b-spline of the median to ensure it is smooth
; and to interpolate over any bad values
; Weight anything with disp_med=0 with 0 weight
tmpwght=replicate(1.,nwave)
index=where(disp_med eq 0.,nindex)
if (nindex gt 0) then tmpwght[index]=0.
nord=3
everyn=100
bkpt=0

sset=bspline_iterfit(wave,disp_med,nord=nord,everyn=everyn,bkpt=bkpt,maxrej=0,invvar=tmpwght,yfit=yfit)
; Convert from units of per log-sampled wavelength pixel to units of Angstroms
dlam=fltarr(nwave)
for i=0,nwave-2 do dlam[i]=wave[i+1]-wave[i]
dlam[nwave-1]=dlam[nwave-2]
; Spectral resolution R=lambda/(FWHM)
specres=wave/(yfit*dlam*2.35)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make the basic cube and RSS file headers
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

mkhdr,rsshdr,flux_arr
mkhdr,cubehdr,cube_flux
; Add MaNGA keycards
mdrp_drp3dhead,drp3qual=drp3qual,cubehdr=cubehdr,rsshdr=rsshdr,obsinfo=obsinfo,mgcf_hdr=hdrsave,wave=wave,cubescale=cubescale,cubefwhm=cubefwhm

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Make some broadband images from the data cube
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

u_img=ml_mangatosdssimg(cube_flux,wave,cubehdr,'u',imghdr=u_hdr)
g_img=ml_mangatosdssimg(cube_flux,wave,cubehdr,'g',imghdr=g_hdr)
r_img=ml_mangatosdssimg(cube_flux,wave,cubehdr,'r',imghdr=r_hdr)
i_img=ml_mangatosdssimg(cube_flux,wave,cubehdr,'i',imghdr=i_hdr)
z_img=ml_mangatosdssimg(cube_flux,wave,cubehdr,'z',imghdr=z_hdr)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Write out data files
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Write out the RSS-format file
rssname=outdir+outnameroot+'LOGRSS.fits'
ml_mwrfits, dummyext, rssname, hdr=rsshdr, /create
ml_mwrfits, flux_arr, rssname, extname='FLUX'
ml_mwrfits, ivar_arr, rssname, extname='IVAR'
ml_mwrfits, mask_arr, rssname, extname='MASK'
ml_mwrfits, disp_arr, rssname, extname='DISP'
ml_mwrfits, wave, rssname, extname='WAVE'
ml_mwrfits, specres, rssname, extname='SPECRES'
ml_mwrfits, obsinfo, rssname, extname='OBSINFO'
ml_mwrfits, x_arr, rssname, extname='XPOS'
ml_mwrfits, y_arr, rssname, extname='YPOS'

; Write out the cube
cubename=outdir+outnameroot+'LOGCUBE.fits'
ml_mwrfits, dummyext, cubename, hdr=cubehdr, /create
; Need to stick the full header in the FLUX ext to in order to make
; some cube viewers (ml2) happy.
ml_mwrfits, cube_flux, cubename, hdr=cubehdr, extname='FLUX'
ml_mwrfits, cube_ivar, cubename, extname='IVAR'
ml_mwrfits, cube_mask, cubename, extname='MASK'
ml_mwrfits, wave, cubename, extname='WAVE'
ml_mwrfits, specres, cubename, extname='SPECRES'
ml_mwrfits, obsinfo, cubename, extname='OBSINFO'
ml_mwrfits, x_arr, cubename, extname='XPOS'
ml_mwrfits, y_arr, cubename, extname='YPOS'
ml_mwrfits, u_img, cubename, hdr=u_hdr, extname='UIMG'
ml_mwrfits, g_img, cubename, hdr=g_hdr, extname='GIMG'
ml_mwrfits, r_img, cubename, hdr=r_hdr, extname='RIMG'
ml_mwrfits, i_img, cubename, hdr=i_hdr, extname='IIMG'
ml_mwrfits, z_img, cubename, hdr=z_hdr, extname='ZIMG'
ml_mwrfits, psfcube[*,*,0], cubename, hdr=u_hdr, extname='UPSF'
ml_mwrfits, psfcube[*,*,1], cubename, hdr=g_hdr, extname='GPSF'
ml_mwrfits, psfcube[*,*,2], cubename, hdr=r_hdr, extname='RPSF'
ml_mwrfits, psfcube[*,*,3], cubename, hdr=i_hdr, extname='IPSF'
ml_mwrfits, psfcube[*,*,4], cubename, hdr=z_hdr, extname='ZPSF'

; gzip output files
spawn, ['gzip', '-f', rssname], /noshell
spawn, ['gzip', '-f', cubename], /noshell

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; DRL- Need to do the LINEAR and MINIMAL stuff still
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

return
end
