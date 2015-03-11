;+
; NAME:
;   mdrp_calib
;
; PURPOSE:
;   Extract calibration frames.
;
; CALLING SEQUENCE:
;   mdrp_calib, flatname, arcname, slitmap=, fibermask=, cartid=, $
;    lampfile=, indir=, timesep=, ecalibfile=, plottitle=, $
;    minflat=, maxflat=, arcinfoname=, flatinfoname=, $
;    arcstruct=, flatstruct=, writeflatmodel=, /bbspec]
;
; INPUTS:
;   flatname   - Name(s) of flat-field SDSS image(s)
;   arcname    - Name(s) of arc SDSS image(s)
;   cartid     - Cartridge ID from plugmap
;
; OPTIONAL KEYWORDS:
;   slitmap    - Full slitmap for the set.
;   fibermask  - Fiber status bits, set nonzero for bad status [NFIBER].
;                Note this is not modified, but modified copies appear
;                in the returned structures ARCSTRUCT and FLATSTRUCT.
;   lampfile   - Name of file describing arc lamp lines, which would
;                over-ride the default file read by FITARCIMAGE.
;   indir      - Input directory for FLATNAME, ARCNAME, OBJNAME;
;                default to '.'
;   timesep    - Maximum time separation between flats and arcs to pair them;
;                set to zero to disable this test; default to 7200 sec.
;   ecalibfile - opECalib file to pass to SDSSPROC
;   plottitle  - Prefix for titles in QA plots.
;   minflat    - Parameter for SDSSPROC for pixel flats; default to 0.8
;   maxflat    - Parameter for SDSSPROC for pixel flats; default to 1.2
;   arcinfoname- File name (with path) to output arc extraction and fitting
;                information
;   flatinfoname-File name (with path) to output flat field extraction and
;                fitting information
; writeflatmodel-Set this keyword to write flat data image, ivar, and
;                final extraction model image to a file.  Will only
;                work if "flatinfoname" is present also (ASB).
; writearcmodel- Set this keyword to write arc data image, ivar, and
;                final extraction model image to a file.  Will only
;                work if "arcinfoname" is present also (ASB).
;   bbspec         - use bbspec extraction code
;
; OUTPUTS:
;   arcstruct  - Structure array with extracted arc calibration information
;   flatstruct - Structure array with extracted flat calibration information
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Always pair arcs to the nearest good flat, and flats to the nearest good arc
;   (nearest in time, as defined by the TAI-BEG keyword in the FITS headers).
;
;   Also store SUPERFLATSET from fiberflat, since we need this to remove
;   small scale features present in all spectra

; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   extract_image
;   fiberflat()
;   fitarcimage
;   fitdispersion
;   fitflatwidth()
;   get_tai
;   reject_arc()
;   reject_flat()
;   sdssproc
;   splog
;   trace320crude()
;   traceset2xy
;   xy2traceset
;
; INTERNAL SUPPORT ROUTINES:
;   create_arcstruct()
;   create_flatstruct()
;
; REVISION HISTORY:
;   24-Jan-2000  Written by D. Schlegel, Princeton
;   27-Nov-2000  Changed to proftype 3, added minflat, maxflat keywords
;    8-Jan-2001  And now back to proftype 1, more robust against bad columns
;   26-Jan-2001  And now let's check both 1&3, and use the better fit
;      Apr-2010  Added "write[flat,arc]model" option (A. Bolton, Utah)
;   25-Jan-2011  Added "twophase" test and switching, A. Bolton, Utah
;   29-Mar-2011  Switched to bundle-wise pure IDL extraction, A. Bolton, Utah
;
;-
;------------------------------------------------------------------------------
function create_arcstruct, narc

   on_error, 0
   compile_opt idl2
   compile_opt idl2, hidden
   
  ftemp = create_struct( name='ARC_STRUCT', $
    'NAME', '', $
    'TAI', 0D, $
    'TSEP', 0D, $
    'QBAD', 0B, $
    'IFLAT', -1, $
    'BESTCORR', 0.0, $
    'NMATCH', 0L, $
    'MEDWIDTH', fltarr(4), $
    'LAMBDA', ptr_new(), $
    'REJLINE', ptr_new(), $
    'XPEAK', ptr_new(), $
    'XDIF_TSET', ptr_new(), $
    'WSET', ptr_new(), $
    'DISPSET', ptr_new(), $
    'FIBERMASK', ptr_new() )

  arcstruct = replicate(ftemp, narc)
  
  return, arcstruct
end
;------------------------------------------------------------------------------
function create_flatstruct, nflat

   on_error, 0
   compile_opt idl2
   compile_opt idl2, hidden
   
  ftemp = create_struct( name='FLAT_STRUCT', $
    'NAME', '', $
    'TAI', 0D, $
    'TSEP', 0D, $
    'QBAD', 0, $
    'IARC', -1, $
    'PROFTYPE', 0, $
    'MEDWIDTH', fltarr(4), $
    'FIBERMASK', ptr_new(), $
    'TSET', ptr_new(), $
    'XSOL', ptr_new(), $
    'WIDTHSET', ptr_new(), $
    'FFLAT', ptr_new(), $
    'SUPERFLATSET', ptr_new() )
    
  flatstruct = replicate(ftemp, nflat)
  
  return, flatstruct
end
;------------------------------------------------------------------------------

pro mdrp_calib, flatname, arcname, slitmap=slitmap, $
    spectrographid=spectrographid, fibermask=fibermask, cartid=cartid, $
    lampfile=lampfile, indir=indir, timesep=timesep, $
    ecalibfile=ecalibfile, plottitle=plottitle, $
    arcinfoname=arcinfoname, flatinfoname=flatinfoname, $
    arcstruct=arcstruct, flatstruct=flatstruct, simulation=simulation, $
    minflat=minflat, maxflat=maxflat, $
    writeflatmodel=writeflatmodel, writearcmodel=writearcmodel, $
    visual=visual, fiberparam=fiberparam, writefiles=writefiles

   on_error, 0
   compile_opt idl2
    
  if (~keyword_set(indir)) then indir = '.'
  if (~keyword_set(timesep)) then timesep = 7200
  if (~keyword_set(minflat)) then minflat = 0.8
  if (~keyword_set(maxflat)) then maxflat = 1.2
  
  stime1 = systime(1)

;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
; STAGE 1 CALIBRATION - FIBER TRACING AND INITIAL FLUX EXTRACTION FOR LSF WIDTH DETERMINATION
  
  ;---------------------------------------------------------------------------
  ; Determine spectrograph ID and color from first flat file
  ;---------------------------------------------------------------------------
  
  sdssproc, flatname[0], indir=indir, spectrographid=spectrographid, color=color

  nflat = n_elements(flatname)
    
  ;---------------------------------------------------------------------------
  ; LOOP THROUGH FLATS + TRACE
  ;---------------------------------------------------------------------------
  splog, 'LOOP THROUGH FLATS-----------------------------'  
  
  flatstruct = create_flatstruct(nflat)
  
  for iflat=0, nflat-1 do begin
  
    splog, iflat+1, nflat, format='("Tracing flat #",I3," of",I3)'
    
    ;---------------------------------------------------------------------
    ; Read flat-field image
    ;---------------------------------------------------------------------

    if ~keyword_set(SIMULATION) then begin     
        splog, 'Reading flat ', flatname[iflat]
        sdssproc, flatname[iflat], flatimg, flativar, indir=indir, hdr=flathdr, /applybias, /applypixflat, nsatrow=nsatrow, fbadpix=fbadpix,$
          ecalibfile=ecalibfile, minflat=minflat, maxflat=maxflat,/applycrosstalk
          
        configuration=obj_new('configuration', sxpar(flathdr, 'MJD'))
    
        ; Decide if this flat is bad
        qbadflat = reject_flat(flatimg, flathdr, nsatrow=nsatrow, fbadpix=fbadpix, percent80thresh=configuration->spcalib_reject_calib_percent80thresh())
        ; When using the whitelights, we need a hack to turn
        ; off the usual flatfield QA
;        qbadflat=0; HACK for whitelights

    endif else begin
        flatimg = mrdfits(indir+flatname[iflat]+'.gz',0,flathdr)
        flativar = mrdfits(indir+flatname[iflat]+'.gz',2)
        qbadflat=0
        configuration=obj_new('configuration', sxpar(flathdr, 'MJD'))
    endelse

    flatinfofile = string(format='(a,i8.8,a)',flatinfoname, sxpar(flathdr, 'EXPOSURE'), '.fits')
    flattmp=fileandpath(flatinfofile,path=outdir)
    
    ;if set then write out the 2d bias-subtraced, pixel-flatted flat field
    if keyword_set(writefiles) then begin
      rawflatname=ml_strreplace(flatinfofile,'mgFlat','mgFlat-2d')
      mwrfits, flatimg, rawflatname, /create
      splog, 'Writing 2d flat image to file: ', rawflatname 
    endif
 
    if (~keyword_set(fibermask)) then tmp_fibmask = 0 else tmp_fibmask = fibermask

    if (~qbadflat) then begin
      ;------------------------------------------------------------------
      ; Create spatial tracing from flat-field image
      ;------------------------------------------------------------------
      splog, 'Tracing fibers in ', flatname[iflat]
      xsol = mdrp_tracefibers(flatimg, flativar, yset=ycen, maxdev=1.0, fibermask=tmp_fibmask, cartid=cartid, xerr=xerr, flathdr=flathdr, miss=miss, slitmap=slitmap, $
       padding=configuration->spcalib_tracefiber_padding(), plottitle=plottitle+' Traces '+flatname[iflat],visual=visual,fiberparam=fiberparam, simulation=simulation,outdir=outdir)

      ;retrieve bundle and fiber info
      nbundle = fiberparam.nbundle
      bundleid = fiberparam.bundleid
      nfiber=fiberparam.nfiber
      radius = fiberparam.radius
        
      splog, 'Fitting traces in ', flatname[iflat]
      ntrace = (size(xsol, /dimens))[1]
      outmask = 0
      ; Ignore values whose central point falls on a bad pixel
      ; ASB: New recipe for inmask, just masking fully useless rows, since traceFibers has already done clever fill-ins:
      inmask = (total(flativar gt 0., 1) gt 0.) # replicate(1B, ntrace)
      xy2traceset, ycen, xsol, tset, ncoeff=configuration->spcalib_xy2traceset_ncoeff(), maxdev=0.5, outmask=outmask, /double, xerr=xerr, inmask=inmask

      junk = where(outmask EQ 0, totalreject)
      if (totalreject GT configuration->spcalib_rejecttheshold()) then begin
        splog, 'Reject flat ' + flatname[iflat] + ': ' + string(format='(i8)', totalreject) + ' rejected pixels'
        qbadflat = 1
      endif
      
      traceset2xy, tset, ycen, xsol
      flatstruct[iflat].tset = ptr_new(tset)
      flatstruct[iflat].xsol = ptr_new(xsol)
      flatstruct[iflat].fibermask = ptr_new(tmp_fibmask)
    endif else begin
      ;bad flat
      xsol = 0
      flatstruct[iflat].qbad = 1
    endelse
    
    ;----------
    ; Verify that traces are separated by > 3 pixels   ;modify this for 2 pixel radius for sims -> to varying radius
    if (qbadflat EQ 0) then begin
      ;print, 'CHECKING IF FIBERS SEPARATED BY > 3 pixels - PRODUCES WARNING if TOO CLOSE'
      sep = xsol[*,1:ntrace-1] - xsol[*,0:ntrace-2]
      tooclose = where(sep LT 3.0)
      if (tooclose[0] NE -1) then begin
        splog, 'Reject flat ' + flatname[iflat] + ': Traces not separated by more than '+strtrim(3.0,2)+' pixels'
      endif
    endif
    
    print, 'QBADFLAT: ',qbadflat
    if (~qbadflat) then begin
      ;---------------------------------------------------------------------
      ; Extract the flat-field image to obtain width and flux
      ;---------------------------------------------------------------------

      print, 'EXTRACTING FLAT FIELD----------------------'
      sigma = configuration->spcalib_sigmaguess() ; Initial guess for gaussian width
      highrej = 15
      lowrej = 15
      npoly = 10 ; Fit 1 terms to background
      wfixed = [1,1] ; Fit the first gaussian term + gaussian width

      ;BOSS better with proftype=1
      proftype = configuration->spcalib_extract_image_proftype()  
      splog, 'Extracting flat with proftype=', proftype

      extract_image, flatimg, flativar, xsol, sigma, flux, fluxivar, proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
                      npoly=npoly, relative=1, ansimage=ansimage, reject=[0.1, 0.6, 0.6], chisq=chisq3, visual=visual, surve=survey,ymodel=ymodel

      cam=strmid(color,0,1)+strtrim(spectrographid,2)
      qawidthfile = outdir+'qa/mgFlat-'+cam+'-00'+strtrim(sxpar(flathdr, 'EXPOSURE'),2)+'-flatwidths.ps'
      widthset3 = mdrp_fitflatwidth(flux, fluxivar, ansimage, tmp_fibmask, ncoeff=configuration->spcalib_fitflatwidth_ncoeff(), sigma=sigma, $
                            medwidth=medwidth, mask=configuration->spcalib_fitflatwidth_mask(flux,fluxivar), plotname=qawidthfile,$
                            inmask=configuration->spcalib_fitflatwidth_inmask(flux,fluxivar,ntrace), /double, nbundle=nbundle,nfiber=nfiber)
      widthset = widthset3
      
      ;relieve some memory
      ansimage = 0
      widthset3 = 0

      junk = where(flux GT 1.0e5, nbright)  ;select out bright pixels     
      splog, 'Using proftype=', proftype
      splog, 'Found ', nbright, ' bright pixels in extracted flat ', flatname[iflat], format='(a,i7,a,a)'
        
      flatstruct[iflat].proftype  = proftype
      flatstruct[iflat].fibermask = ptr_new(tmp_fibmask)
      flatstruct[iflat].widthset = ptr_new(widthset)
      flatstruct[iflat].medwidth  = medwidth
      
    endif
    
    flatstruct[iflat].name = flatname[iflat]
    get_tai, flathdr, tai_beg, tai_mid, tai_end
    flatstruct[iflat].tai = tai_mid
    flatstruct[iflat].qbad = qbadflat
    obj_destroy,configuration
    
    if flatstruct[iflat].qbad eq 1 then splog, 'WARNING: QBAD is 1 for FLAT, FLAT REJECT' else splog, 'GOOD FLAT!, '+flatname[iflat]    
    
    ;write out the flat field extraction after the 1st pass -- extraction and determination of profile widths
    if keyword_set(writefiles) then begin
      flat1extname = ml_strreplace(flatinfofile,'mgFlat','mgFlat-1flux')
      mwrfits, flux, flat1extname, /create
      flat1modname=ml_strreplace(flatinfofile,'mgFlat','mgFlat-1model')
      mwrfits, ymodel, flat1modname, /create  
      spawn, ['gzip', '-f', flat1extname], /noshell
      spawn, ['gzip', '-f', flat1modname], /noshell     
    endif

  endfor

;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
; STAGE 2 CALIBRATION - WAVELENGTH CALIBRATION USING ARCS
 
  ;--------REDUCE THE ARCS
    narc = n_elements(arcname)
  
    ;---------------------------------------------------------------------------
    ; LOOP THROUGH ARCS + FIND WAVELENGTH SOLUTIONS
    ;---------------------------------------------------------------------------
    splog, 'LOOP THROUGH ARCS-------------------------'
    
    arcstruct = create_arcstruct(narc)
    
    for iarc=0, narc-1 do begin
    
      splog, iarc+1, narc, format='("Extracting arc #",I3," of",I3)'
      
      ;---------------------------------------------------------------------
      ; Read the arc
      ;---------------------------------------------------------------------
  
      if ~keyword_set(SIMULATION) then begin    
          splog, 'Reading arc ', arcname[iarc]
          
          sdssproc, arcname[iarc], arcimg, arcivar, indir=indir, hdr=archdr, /applybias, /applypixflat, nsatrow=nsatrow, fbadpix=fbadpix, $
            ecalibfile=ecalibfile, minflat=minflat, maxflat=maxflat,/applycrosstalk
          ny = (size(arcimg,/dimens))[1]
            
          configuration=obj_new('configuration', sxpar(archdr, 'MJD'))
          splog, 'Fraction of bad pixels in arc = ', fbadpix
          
          ;----------
          ; Decide if this arc is bad
          qbadarc = reject_arc(arcimg, archdr, nsatrow=nsatrow, fbadpix=fbadpix)
          
      endif else begin
          arcimg = mrdfits(indir+arcname[iarc]+'.gz',0,archdr)
          arcivar = mrdfits(indir+arcname[iarc]+'.gz',2)
          qbadarc=0  
          ny = (size(arcimg,/dimens))[1]
          fbadpix=0
          configuration=obj_new('configuration', sxpar(archdr, 'MJD')) 
      endelse

     ;----------
     ; Identify the nearest flat-field for this arc, which must be within TIMESEP seconds and be a good flat.      
      get_tai, archdr, tai_beg, tai_mid, tai_end
      tai = tai_mid
      
      iflat = -1
      igood = where(flatstruct.qbad EQ 0)
      if (igood[0] NE -1) then begin
        tsep = min( abs(tai - flatstruct[igood].tai), ii )
        if (tsep LE timesep AND timesep NE 0) then iflat = igood[ii]
      endif
      
      if (iflat GE 0) then begin
        print, 'FOUND FLAT TO PAIR WITH ARC!'
        splog, 'Arc ' + arcname[iarc] + ' paired with flat ' + flatname[iflat]
      endif else begin
        print, 'NO FLAT TO PAIR WITH ARC, SETTING BAD ARC to 1'
        splog, 'Arc ' + arcname[iarc] + ' paired with no flat'
        qbadarc = 1
      endelse
      
      if (~qbadarc) then begin
        xsol = *(flatstruct[iflat].xsol)
        widthset = *(flatstruct[iflat].widthset)
        tmp_fibmask = *(flatstruct[iflat].fibermask)
        proftype = flatstruct[iflat].proftype
        
        ;----------
        ; Calculate possible shift between arc and flat
        xcor = mdrp_match_trace(arcimg, arcivar, xsol, radius=radius, nfiber=nfiber)
        
        bestlag = median(xcor-xsol)
        if (abs(bestlag) GT 2.0) then begin
          qbadarc = 1
          splog, 'Reject arc: pixel shift is larger than 2 pixel'
          splog, 'Reject arc ' + arcname[iarc] + ': Pixel shift = ', bestlag
        endif
      endif
      
      if (~qbadarc) then begin
        splog, 'Shifting traces with match_trace', bestlag
     
        ;---------------------------------------------------------------------
        ; Extract the arc image
        ;---------------------------------------------------------------------

        traceset2xy, widthset, xx, traced_sigma
        
        highrej = 15
        lowrej = 15
        wfixed = [1,0] ; ASB: Don't fit for width terms.
  
        splog, 'Extracting arc'
        mdrp_extract_bundle_image, arcimg, arcivar, xcor, traced_sigma, flux, fluxivar, proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, npoly=2L, relative=1, $
          reject=[0.1, 0.6, 0.6], ymodel=ymodel, nperbun=20L, buffsize=8L, visual=visual, survey=survey, nbundle=nbundle, nfiber=nfiber
          
        ;flag to determine whether or not to do 2-phase arc solution:
        twophase = sxpar(archdr, 'TWOPHASE')
        if keyword_set(twophase) then splog, 'Setting 2-phase readout flag'
  
        if keyword_set(VISUAL) then begin
          getwindow,/open
          h=bytscl(flux,min=0,max=255)
          ml_tvimage, h, /axis, axkeywords={charsize:2,xtitle:'Row',ytitle:'Fiber Number',title:'Extracted Arc Image - '+string(f='(a1,i1)',color, spectrographid)}
        endif
  
        ;---------------------------------------------------------------------
        ; Compute correlation coefficient for this arc image
        ;---------------------------------------------------------------------
  
        splog, 'Searching for wavelength solution'
        aset = 0

        mdrp_fitarcimage, flux, fluxivar, aset=aset, color=color, lampfile=lampfile, fibermask=tmp_fibmask, bestcorr=bestcorr, $
                    acoeff=configuration->spcalib_arcfitguess_acoeff(color), dcoeff=configuration->spcalib_arcfitguess_dcoeff(color), $
                    wrange=configuration->spcalib_fitarcimage_wrange(color), twophase=twophase, visual=visual, $
                    nbundle=nbundle, nfiber=nfiber, radius=radius
  
        arcstruct[iarc].bestcorr = bestcorr
        
        if ((color EQ 'blue' AND bestcorr LT 0.5) OR (color EQ 'red'  AND bestcorr LT 0.5) ) then begin
          qbadarc = 1
          splog, 'Reject arc ' + arcname[iarc] + ': correlation is only = ' + string(format='(i4)', bestcorr)
        endif
      endif
      
      ;if we don't have a bad arc then continue
      if (~qbadarc) then begin
      
        ; Check for unplugged IFU's / missing entire blocks before doing wavelength calibration
        ml_checkunpluggedifus, miss, nfiber, newfiber=newfiber, newblock=newblock, fibers=unplugfibers, blockid=unplugblockid
      
        ;---------------------------------------------------------------------
        ; Compute wavelength calibration
        ;---------------------------------------------------------------------
        arccoeff = configuration->spcalib_arccoeff()
        splog, 'Searching for wavelength solution'
        mdrp_fitarcimage, flux[*,newfiber], fluxivar[*,newfiber], xpeak, ypeak, wset, ncoeff=arccoeff, aset=aset, color=color, lampfile=lampfile, $
          fibermask=tmp_fibmask[newfiber], lambda=lambda, rejline=rejline, xdif_tset=xdif_tset, acoeff=configuration->spcalib_arcfitguess_acoeff(color), $
          dcoeff=configuration->spcalib_arcfitguess_dcoeff(color), wrange=configuration->spcalib_fitarcimage_wrange(color), twophase=twophase,visual=visual, $
          nbundle=n_elements(newblock), nfiber=nfiber[newblock], radius=radius[newblock]
  
        splog, 'Number of Rejected lamp lines: ', n_elements(where(rejline eq 'Reject-offset')), ' out of 45 lamp lines'

        ; In case of unplugged IFUs, we need to expand xpeak, ypeak, and the wavelength calibration trace, back to the full complement of NTRACE fibers
        if n_elements(newfiber) ne ntrace then begin
          splog, 'Fixing wavelength calibration from '+string(n_elements(newfiber),' to ',ntrace,f='(I3,A,I3)')+' fibers, due to unplugged IFUs'
          ml_expandtraces, xpeak, ypeak, wset, ntrace, infiber=newfiber, outfiber=unplugfibers, other=xdif_tset
        endif
  
        ;check for no wset
        if (~keyword_set(wset)) then begin
          splog, 'Warning: No wset! Setting qbadarc to 1.  Wavelength solution failed'
          qbadarc = 1
        endif else begin
  
          nfitcoeff = configuration->spcalib_ncoeff(color)
          ilamp = where(rejline EQ '')
          dispset = mdrp_fitdispersion(flux, fluxivar, xpeak[*,ilamp], sigma=configuration->spcalib_sigmaguess(), ncoeff=nfitcoeff, $
                                   xmin=0.0, xmax=ny-1, nfiber=nfiber, medwidth=wsigarr, numbundles=nbundle)
  
          arcstruct[iarc].dispset = ptr_new(dispset)
          arcstruct[iarc].wset = ptr_new(wset)
          arcstruct[iarc].nmatch = N_elements(lambda)
          arcstruct[iarc].lambda = ptr_new(lambda)
          arcstruct[iarc].rejline = ptr_new(rejline)
          arcstruct[iarc].tsep = tsep
          arcstruct[iarc].xpeak = ptr_new(xpeak)
          arcstruct[iarc].xdif_tset = ptr_new(xdif_tset)
          arcstruct[iarc].fibermask = ptr_new(tmp_fibmask)
          arcstruct[iarc].medwidth = wsigarr
          
          ;------------------------------------------------------------------
          ; Write information on arc lamp processing
          
          if (keyword_set(arcinfoname)) then begin
            sxaddpar, archdr, 'FBADPIX', fbadpix, 'Fraction of bad pixels in raw image'
            sxaddpar, archdr, 'BESTCORR', bestcorr, 'Best Correlation coefficient'
              
            arcinfofile = string(format='(a,i8.8,a)',arcinfoname, sxpar(archdr, 'EXPOSURE'), '.fits')
            
            print, 'WRITING ARC: ',arcinfofile
              
            ml_mwrfits, dummyext, arcinfofile, hdr=archdr, /create; Blank ext 0 with full header
            ml_mwrfits, flux, arcinfofile, extname='FLUX'
            ml_mwrfits, [transpose(lambda), xpeak], arcinfofile, extname='LXPEAK'
            ml_mwrfits, *arcstruct[iarc].wset, arcinfofile, extname='WSET'
            ml_mwrfits, *arcstruct[iarc].fibermask, arcinfofile, extname='MASK'
            ml_mwrfits, *arcstruct[iarc].dispset, arcinfofile, extname='DISPSET'
            
            spawn, ['gzip', '-f', arcinfofile], /noshell

           ;write arc image model info if requested:
            ; DRL- MODELIMG needs to be added to data model
            if keyword_set(writearcmodel) then begin
               arcmodelfile = string(format='(a,i8.8,a)',arcinfoname + 'MODELIMG-', sxpar(archdr, 'EXPOSURE'), '.fits')
               ml_mwrfits, dummyext, arcmodelfile, /create
               ml_mwrfits, arcimg, arcmodelfile, extname='ARCIMG'
               ml_mwrfits, arcivar, arcmodelfile, extname='ARCIVAR'
               ml_mwrfits, ymodel, arcmodelfile, extname='YMODEL'
               spawn, ['gzip', '-f', arcmodelfile], /noshell
            endif
            ymodel = 0
          endif  ;end of write arc file
        endelse ;end of no wset if statement
      endif  ;end of qbadarc if 
      
      arcstruct[iarc].name = arcname[iarc]
      arcstruct[iarc].tai = tai
      arcstruct[iarc].iflat = iflat
      arcstruct[iarc].qbad = qbadarc
      
      if arcstruct[iarc].qbad eq 1 then splog, 'WARNING: QBAD FOR ARC = 1, ARC REJECT' else splog, 'GOOD ARC!, '+arcname[iarc]
      
      obj_destroy,configuration
    endfor
    
    arcimg = 0
    arcivar = 0
    
;-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
; STAGE 3 CALIBRATION - FINAL FLUX EXTRACTION FROM FLATS + FIBER-TO-FIBER VARIATION REMOVAL + SUPERFLAT REMOVAL

  ;---------------------------------------------------------------------------
  ; LOOP THROUGH FLATS + CREATE FIBERFLATS
  ;---------------------------------------------------------------------------
  splog, 'CREATING FIBERFLATS----------------------'
  for iflat=0, nflat-1 do begin
  
    splog, iflat+1, nflat, format='("Create fiberflats for flat #",I3," of",I3)'
      
    ;----------
    ; Identify the nearest arc for each flat-field, which must be within TIMESEP seconds and be good.
    iarc = -1
    igood = where(arcstruct.qbad EQ 0)
    if (igood[0] NE -1) then begin
      tsep = min( abs(flatstruct[iflat].tai - arcstruct[igood].tai), ii )
      if (tsep LE timesep AND timesep NE 0) then iarc = igood[ii]
      flatstruct[iflat].tsep = tsep
    endif
    
    if (iarc GE 0) then splog, 'Flat ' + flatname[iflat] + ' paired with arc ' + arcname[iarc] $
    else begin
      splog, 'Flat ' + flatname[iflat] + ' paired with no arc'
      flatstruct[iflat].qbad = 1 ; Flat is bad if no companion arc exists
    endelse
    
    flatstruct[iflat].iarc = iarc
    
    if (~flatstruct[iflat].qbad) then begin
    
      widthset = *(flatstruct[iflat].widthset)
      wset = *(arcstruct[iarc].wset)
      xsol = *(flatstruct[iflat].xsol)
      tmp_fibmask = *(flatstruct[iflat].fibermask)
      proftype = flatstruct[iflat].proftype
      
      ;---------------------------------------------------------------------
      ; Read flat-field image (again)
      ;---------------------------------------------------------------------
      
      ; If there is only 1 flat image, then it's still in memory
      if (nflat GT 1) then begin
        splog, 'Reading flat ', flatname[iflat]
        sdssproc, flatname[iflat], flatimg, flativar, indir=indir, hdr=flathdr, /applybias, /applypixflat, ecalibfile=ecalibfile, minflat=minflat, maxflat=maxflat,/applycrosstalk
      endif
      configuration=obj_new('configuration',sxpar(flathdr, 'MJD'))

      ;---------------------------------------------------------------------
      ; Extract the flat-field image
      ;---------------------------------------------------------------------
      
      traceset2xy, widthset, xx, sigma2   ; sigma2 is real width

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

       extract_image, flatimg, flativar, xsol, sigma2, flux, fluxivar, proftype=proftype, wfixed=wfixed, highrej=highrej, lowrej=lowrej, $
         npoly=npoly, relative=1, ansimage=ansimage, reject=reject, chisq=chisq3, visual=visual, survey=survey,ymodel=ymodel

      ;write out the flat field extraction after the 2nd pass -- actual flux extraction over blocks or rows 
      if keyword_set(writefiles) then begin
        flat2extname=ml_strreplace(flatinfofile,'mgFlat','mgFlat-2flux')
        mwrfits, flux, flat2extname, /create
        flat2modname=ml_strreplace(flatinfofile,'mgFlat','mgFlat-2model')
        mwrfits, ymodel, flat2modname, /create  
        spawn, ['gzip', '-f', flat2extname], /noshell
        spawn, ['gzip', '-f', flat2modname], /noshell 
      endif
        
      ;---------------------------------------------------------------------
      ; Compute fiber-to-fiber flat-field variations
      ;---------------------------------------------------------------------
      sigma2 = 0
      xsol = 0

      qasflatfile = outdir+'qa/mgFlat-'+cam+'-00'+strtrim(sxpar(flathdr, 'EXPOSURE'),2)+'-superflat.ps'

      fflat = mdrp_fiberflat(flux, fluxivar, wset, fibermask=tmp_fibmask, /dospline, pixspace=5, plotname=qasflatfile, plottitle=plottitle+' Superflat '+flatstruct[iflat].name, $
        badflatfracthresh=configuration->spcalib_fiberflat_badflatfracthresh(), minval=configuration->spcalib_fiberflat_minval(flux), visual=visual, $
        fiberparam=fiberparam, superflatset=superflatset, medval=medval)
        
      if (n_elements(fflat) EQ 1) then begin
        flatstruct[iflat].qbad  = 1
        splog, 'Reject flat ' + flatname[iflat] + ': No good traces'
      endif
      
      flatstruct[iflat].fflat = ptr_new(fflat)
      flatstruct[iflat].superflatset = ptr_new(superflatset)
      flatstruct[iflat].fibermask = ptr_new(tmp_fibmask)
      
      ;------------------------------------------------------------------
      ; Write information on flat field processing
      if (keyword_set(flatinfoname)) then begin
        sxaddpar, flathdr, 'NBRIGHT', nbright, 'Number of bright pixels (>10^5) in extracted flat-field'
        flatinfofile = string(format='(a,i8.8,a)',flatinfoname, sxpar(flathdr, 'EXPOSURE'), '.fits')
        print, 'WRITING FLAT: ',flatinfofile

        ml_mwrfits, dummyext, flatinfofile, hdr=flathdr, /create
        ml_mwrfits, *flatstruct[iflat].fflat, flatinfofile, extname='FLUX'
        ml_mwrfits, *flatstruct[iflat].tset, flatinfofile, extname='TSET'
        ml_mwrfits, *flatstruct[iflat].fibermask, flatinfofile, extname='MASK'
        ml_mwrfits, *flatstruct[iflat].widthset, flatinfofile, extname='DISPSET'
        ml_mwrfits, *flatstruct[iflat].superflatset, flatinfofile, extname='SUPERFLATSET'
        
        spawn, ['gzip', '-f', flatinfofile], /noshell

        ; Make QA plots for the flat per ifu
        ; Plot file
        qafile=outdir+'qa/mgFlat-'+cam+'-00'+strtrim(sxpar(flathdr, 'EXPOSURE'),2)+'-ifuflat.ps'
        ; Yanny .par file (of the format ifuflat-02-b1-56741-177343.par)
        thecart=fix(sxpar(flathdr, 'CARTID'))
        cartstr=strtrim(string(thecart),2)
        if (thecart lt 10) then cartstr='0'+cartstr
        theexposure=strtrim(sxpar(flathdr, 'EXPOSURE'),2)
        themjd=strtrim(sxpar(flathdr, 'MJD'),2)
        qapar=outdir+'qa/ifuflat-'+cartstr+'-'+cam+'-'+themjd+'-'+theexposure+'.par'
	ml_qaifuflat,*flatstruct[iflat].fflat,*flatstruct[iflat].superflatset,*arcstruct[0].wset,slitmap,cam,flathdr,plotname=qafile,parfile=qapar

        ;write the flat image model if requested
        ; DRL- needs to be added to data model
        if keyword_set(writeflatmodel) then begin
           flatmodelfile = string(format='(a,i8.8,a)',flatinfoname + 'MODELIMG-', sxpar(flathdr, 'EXPOSURE'), '.fits')
           ml_mwrfits, dummyext, flatmodelfile, /create
           ml_mwrfits, flatimg, flatmodelfile, extname='FLATIMG'
           ml_mwrfits, flativar, flatmodelfile, extname='FLATIVAR'
           ml_mwrfits, ymodel, flatmodelfile, extname='YMODEL'
           spawn, ['gzip', '-f', flatmodelfile], /noshell
        endif
        ymodel = 0
      endif
      
      obj_destroy,configuration
    endif
  endfor

  ;write out calibration structures for the case where you want to skip the reduction of them next time
  tmp=fileandpath(flatinfoname,path=outdir)
  save, flatstruct, arcstruct, fiberparam, filename=outdir+'calibstruct-'+string(format='(a1,i1)',color,spectrographid)+'.sav', description='structures containing the flat and arc parameters'
  splog, 'Writing calibration structure ', outdir+'calibstruct-'+string(format='(a1,i1)',color,spectrographid)+'.sav'
  
  splog, 'Elapsed time = ', systime(1)-stime1, ' seconds', format='(a,f6.0,a)'
  return
end
;------------------------------------------------------------------------------
