;+
; NAME:
;   mdrp_fluxcal
;
; PURPOSE:
;    Flux calibrate sky-subtracted MaNGA frames (mgSFrame) on a 
;    per-camera basis.  This routine needs to be run once
;    per spectrograph; innames should be
;    an N_exposure x N_camera array of names pointing to the 
;    corresponding mgSFrame files. N_camera=2
;
;   b1calib is essentially like extension 0 of the
;   BOSS files spFluxcalib, which is calibimg = float( bspline_valu(loglam1, thisset, x2=x2) )
; CALLING SEQUENCE:
;   ...
;
; INPUTS:
;   innames - filenames for one or many exposures of a plate in one of the 
;            spectrographs. 
;            This routine needs to be run per spectrograph. innames should be
;            an N_exposure x N_camera array of names pointing to the 
;            corresponding mgSFrame files. N_camera=2
;
; OPTIONAL INPUTS:
;   plotname - File for the output plot (overrides default)
;   /visual  - Plots go to X-display instead of file
;
; OPTIONAL OUTPUTS:
;   outnames - Array of names of the output mgFFrame files
;   bcalib   - 
;   rcalib   -
;
; KEYWORDS
;   writefile - when specified, it outputs the stellar model, mratio, mrativar, etc.
;               for each star in a file named fluxcal-sp[1,2]-expno.fits
;
; PROCEDURES CALLED:
;   ...
;
; REVISION HISTORY:
;   25-Mar-2014 First written based on BOSS modules (Renbin Yan)
;   06-May-2014 Tweaked to write out mgFFrame files (David Law)
;   15-May-2014 Added quality control bitmask logic (David Law)
;-

function spflux_mratio_flatten, loglam1, mratio1, mrativar1, pres=pres

   ;--------
   ; Re-form the input data arrays from multi-dimensional to N x M

   ndim = size(loglam1, /n_dimen)
   dims = size(loglam1, /dimens)
   npix = dims[0]
   nobj = n_elements(loglam1) / npix
   loglam = reform(loglam1, npix, nobj)
   mratio = reform(mratio1, npix, nobj)
   mrativar = reform(mrativar1, npix, nobj)

   ;--------
   ; Re-bin the spectra to the same spacing

   minlog1 = min(loglam, max=maxlog1)
   newloglam = wavevector(minlog1, maxlog1)
   nnewpix = n_elements(newloglam)

   newratio = fltarr(nnewpix, nobj)
   newivar = fltarr(nnewpix, nobj)

   for iobj=0L, nobj-1 do begin
      isort = sort(loglam[*,iobj])
      combine1fiber, loglam[isort,iobj], mratio[isort,iobj], $
       mrativar[isort,iobj], $
       newloglam=newloglam, newflux=newratio1, newivar=newivar1
      newratio[*,iobj] = newratio1
      newivar[*,iobj] = newivar1
   endfor

   ;--------
   ; Compute the straight weighted mean at each wavelength
   ; (Avoid divide-by- zeros.)

   if (ndim EQ 1) then begin
      meanratio = (newratio * newivar) / (newivar + (newivar EQ 0))
   endif else begin
      denom = total(newivar, 2)
      meanratio = total(newratio * newivar, 2) / (denom + (denom EQ 0))
   endelse

   qbadpix = meanratio LE 0
   ibadpix = where(qbadpix, nbadpix)
   if (nbadpix GT 0) then newivar[ibadpix,*] = 0

   ;--------
   ; Actually take this "mean" and turn it into something more like
   ; a median, to protect us against standard stars that have bad
   ; magnitudes from the imaging.

; Comment-out ???
;   igoodpix = where(qbadpix EQ 0)
;   if (ndim EQ 1) then medratio = newratio $
;    else medratio = djs_median(newratio, 2)
;   rescale = median( medratio[igoodpix] / meanratio[igoodpix] )
;   if (rescale LE 0) then begin
;      splog, 'Warning: RESCALE = ', rescale
;   endif else begin
;      meanratio = rescale * meanratio
;      splog, 'Rescale factor median/mean = ', rescale
;   endelse

   ;--------
   ; Now for each object, compute the polynomial fit of it relative to the mean

   npoly = 4 ; ???
   flatarr = fltarr(npix, nobj)
   pres = fltarr(npoly, nobj)
   for iobj=0L, nobj-1 do begin
      ii = where(newivar[*,iobj] GT 0, ct)
      if (ct GT npoly+1) then begin ; At least NPOLY+1 pixels for a fit...
         thisloglam = newloglam[ii]
         thisratio = newratio[ii,iobj] / meanratio[ii]
         thisivar = newivar[ii,iobj] * meanratio[ii]^2

         ; This fit requires no rejection, because this function falls
         ; within an iteration loop that rejects points.

         ; The following is a weighted fit...
;         pres1 = poly_fit(thisloglam-3.5d0, thisratio, npoly-1, $
;          measure_errors=1./sqrt(thisivar))

         ; The following would be an unweighted fit...
;         pres1 = poly_fit(thisloglam-3.5d0, thisratio, npoly-1)

         ; The following is an unweighted fit but with outlier-rejection...
         poly_iter, thisloglam-3.5d0, thisratio, npoly-1, 3., coeff=pres1

         flatarr[*,iobj] = poly(loglam[*,iobj]-3.5d0, pres1)
         pres[*,iobj] = reform(pres1, npoly)
       endif else begin
         flatarr[*,iobj] = 1
         pres[*,iobj] = 0
         pres[0,iobj] = 1
       endelse
   endfor

   if (ndim GT 1) then $
    pres = reform(pres, [npoly, dims[1:ndim-1]])
   return, reform(flatarr, dims)
end

;-------------------------------------------------------------------------
function spflux_masklines, loglam, hwidth=hwidth, stellar=stellar, $
 telluric=telluric

   if (NOT keyword_set(hwidth)) then $
    hwidth = 5.7e-4 ; Default is to mask +/- 5.7 pix = 400 km/sec

   mask = bytarr(size(loglam,/dimens))

   if (keyword_set(stellar)) then begin
      starwave = [ $
       3830.0 , $ ; ? (H-7 is at 3835 Ang)
       3889.0 , $ ; H-6
       3933.7 , $ ; Ca_k
       3968.5 , $ ; Ca_H (and H-5 at 3970. Ang)
       4101.7 , $ ; H-delta
       4300.  , $ ; G-band
       4305.  , $ ; G-band
       4310.  , $ ; more G-band
       4340.5 , $ ; H-gamma
       4861.3 , $ ; H-beta
       5893.0 , $ ; Mg
       6562.8 , $ ; H-alpha
       8500.8 , $
       8544.6 , $
       8665.0 , $
       8753.3 , $
       8866.1 , $
       9017.5 , $
       9232.0 ]
      airtovac, starwave

      for i=0L, n_elements(starwave)-1 do begin
         mask = mask OR (loglam GT alog10(starwave[i])-hwidth $
          AND loglam LT alog10(starwave[i])+hwidth)
      endfor
   endif

   if (keyword_set(telluric)) then begin
      tellwave1 = [6850., 7150., 7560., 8105., 8930.]
      tellwave2 = [6960., 7350., 7720., 8240., 9030.]
      for i=0L, n_elements(tellwave1)-1 do begin
         mask = mask OR (loglam GT alog10(tellwave1[i]) $
          AND loglam LT alog10(tellwave2[i]))
      endfor
   endif

   return, mask
end
;------------------------------------------------------------------------------
; Divide the spectrum by a median-filtered spectrum.
; The median-filtered version is computed ignoring stellar absorp. features.

function spflux_medianfilt, loglam, objflux, objivar, width=width, $
 newivar=newivar, _EXTRA=KeywordsForMedian

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   npix = dims[0]
   if (ndim EQ 1) then nspec = 1 $
    else nspec = dims[1]

   ;----------
   ; Loop over each spectrum

   medflux = 0 * objflux
   if (arg_present(objivar)) then newivar = 0 * objivar
   for ispec=0L, nspec-1 do begin

      ; For the median-filter, ignore points near stellar absorp. features,
      ; but keep points near telluric bands.
      qgood = 1 - spflux_masklines(loglam[*,ispec], /stellar)

      ; Median-filter, but skipping masked points
      igood = where(qgood, ngood)
      thisback = fltarr(dims[0])
      if (ngood GT 1) then begin
         thisback[igood] = djs_median(objflux[igood,ispec], width=width, $
          _EXTRA=KeywordsForMedian)
      endif
      thisback = djs_maskinterp(thisback, (qgood EQ 0), /const)

      ; Force the ends of the background to be the same as the spectrum,
      ; which will force the ratio of the two to be unity.
      hwidth = ceil((width-1)/2.)
      thisback[0:hwidth] = objflux[0:hwidth,ispec]
      thisback[npix-1-hwidth:npix-1] = objflux[npix-1-hwidth:npix-1,ispec]
      czero2 = where(thisback eq 0., count2)
      if count2 gt 0 then thisback[czero2] = 1.
      medflux[*,ispec] = objflux[*,ispec] / thisback
      if (arg_present(objivar)) then $
      newivar[*,ispec] = objivar[*,ispec] * thisback^2
   endfor

   return, medflux
end
;------------------------------------------------------------------------------
function spflux_bestmodel, loglam, objflux, objivar, dispimg, kindx=kindx1, $
 plottitle=plottitle

   filtsz = 99 ; ???
   cspeed = 2.99792458e5

   ndim = size(objflux, /n_dimen)
   dims = size(objflux, /dimens)
   npix = dims[0]
   if (ndim EQ 1) then nspec = 1 $
    else nspec = dims[1]

   ;----------
   ; Median-filter the object fluxes

   medflux = spflux_medianfilt(loglam, objflux, objivar, $
    width=filtsz, /reflect, newivar=medivar)
   sqivar = sqrt(medivar)

   ;----------
   ; Mask out the telluric bands

   sqivar = sqivar * (1 - spflux_masklines(loglam, /telluric))

   ;----------
   ; Load the Kurucz models into memory

   junk = spflux_read_kurucz(kindx=kindx)
   nmodel = n_elements(kindx)

   ;----------
   ; Fit the redshift just by using a canonical model

   ifud = where(kindx.teff EQ 6000 AND kindx.g EQ 4 AND kindx.feh EQ -1.5)
   if (ifud[0] EQ -1) then $
    message, 'Could not find fiducial model!'
   nshift = 20
   logshift = (-nshift/2. + findgen(nshift)) * 1.d-4
   chivec = fltarr(nshift)
   for ishift=0L, nshift-1 do begin
      modflux = spflux_read_kurucz(loglam-logshift[ishift], $
       dispimg, iselect=ifud)
      ; Median-filter this model
      medmodel = spflux_medianfilt(loglam, modflux, $
       width=filtsz, /reflect)
      for ispec=0L, nspec-1 do begin
         chivec[ishift] = chivec[ishift] + computechi2(medflux[*,ispec], $
          sqivar[*,ispec], medmodel[*,ispec])
      endfor
   endfor
   zshift = (10.d^logshift - 1) ; Convert log-lambda shift to redshift
   zpeak = find_nminima(chivec, zshift, errcode=errcode)
   splog, 'Best-fit velocity for std star = ', zpeak * cspeed, ' km/s'
   if (errcode NE 0) then $
    splog, 'Warning: Error code ', errcode, ' fitting std star'

   ;----------
   ; Generate the Kurucz models at the specified wavelengths + dispersions,
   ; using the best-fit redshift

   modflux = spflux_read_kurucz(loglam-alog10(1.+zpeak), dispimg)

   ;----------
   ; Loop through each model, computing the best chi^2
   ; as the sum of the best-fit chi^2 to each of the several spectra
   ; for this same object.
   ; We do this after a median-filtering of both the spectra + the models.

   chiarr = fltarr(nmodel,nspec)
   chivec = fltarr(nmodel)
   for imodel=0L, nmodel-1 do begin
      ; Median-filter this model
      medmodel = spflux_medianfilt(loglam, modflux[*,*,imodel], $
       width=filtsz, /reflect)

      for ispec=0L, nspec-1 do begin
         chiarr[imodel,ispec] = computechi2(medflux[*,ispec], $
          sqivar[*,ispec], medmodel[*,ispec])
      endfor
      chivec[imodel] = total(chiarr[imodel,*])
   endfor

   ;----------
   ; Return the best-fit model

   minchi2 = min(chivec, ibest)
   dof = total(sqivar NE 0)
   splog, 'Best-fit total chi2/DOF = ', minchi2/(dof>1)
   bestflux = modflux[*,*,ibest]

   ;----------
   ; Compute the chi^2 just around the stellar absorp. lines
   ; for the best-fit model star

   mlines = spflux_masklines(loglam, hwidth=12e-4, /stellar)
   linesqivar = sqivar * mlines
   linechi2 = 0.
   for ispec=0L, nspec-1 do begin
      thismodel = spflux_medianfilt(loglam, modflux[*,ispec,ibest], $
       width=filtsz, /reflect)
      linechi2 = linechi2 + computechi2(medflux[*,ispec], $
       linesqivar[*,ispec], thismodel)
   endfor
   linedof = total(linesqivar NE 0)
   splog, 'Best-fit line chi2/DOF = ', linechi2/(linedof>1)

   ;----------
   ; Compute the median S/N for all the spectra of this object,
   ; and for those data just near the absorp. lines

   sn_median = median(objflux * sqrt(objivar))
   indx = where(mlines, ct)
   if (ct GT 1) then $
    linesn_median = median(objflux[indx] * sqrt(objivar[indx])) $
   else $
    linesn_median = 0.
   splog, 'Full median S/N = ', sn_median
   splog, 'Line median S/N = ', linesn_median

   kindx1 = create_struct(kindx[ibest], $
    'IMODEL', ibest, $
    'Z', zpeak, $
    'SN_MEDIAN', sn_median, $
    'CHI2', minchi2, $
    'DOF', dof, $
    'LINESN_MEDIAN', linesn_median, $
    'LINECHI2', linechi2, $
    'LINEDOF', linedof)

   ;----------
   ; Plot the filtered object spectrum, overplotting the best-fit Kurucz model

   ; Select the observation to plot that has the highest S/N,
   ; and one that goes blueward of 4000 Ang.
   snvec = total(objflux * sqrt(objivar), 1) $
    * (10.^loglam[0,*] LT 4000 OR 10.^loglam[npix-1,*] LT 4000)
   junk = max(snvec, iplot) ; Best blue exposure

   snvec = total(objflux * sqrt(objivar), 1) $
    * (10.^loglam[0,*] GT 8600 OR 10.^loglam[npix-1,*] GT 8600)
   junk = max(snvec, jplot) ; Best red exposure

   csize = 0.85
   djs_plot, [3840., 4120.], [0.0, 1.4], /xstyle, /ystyle, /nodata, $
    xtitle='Wavelength [Ang]', ytitle='Normalized Flux', $
    title=plottitle
   if (iplot[0] NE -1) then begin
      djs_oplot, 10^loglam[*,iplot], medflux[*,iplot]
      djs_oplot, 10^loglam[*,iplot], medmodel[*,iplot], color='red'
   endif
   xyouts, 3860, 1.25, kindx1.model, charsize=csize
   djs_xyouts, 4000, 0.3, charsize=csize, $
    string(minchi2/(dof>1), format='("Total \chi^2/DOF=",f5.2)')
   djs_xyouts, 4000, 0.2, charsize=csize, $
    string(linechi2/(linedof>1), format='("Lines \chi^2/DOF=",f5.2)')
   djs_xyouts, 3860, 0.1, string(kindx1.feh, kindx1.teff, kindx1.g, $
    zpeak*cspeed, $
    format='("Fe/H=", f4.1, "  T_{eff}=", f6.0, "  g=", f3.1, "  cz=",f5.0)'), $
    charsize=csize

   djs_plot, [8440., 9160.], [0.0, 1.4], /xstyle, /ystyle, /nodata, $
    xtitle='Wavelength [Ang]', ytitle='Normalized Flux'
   if (jplot[0] NE -1) then begin
      djs_oplot, 10^loglam[*,jplot], medflux[*,jplot]
      djs_oplot, 10^loglam[*,jplot], medmodel[*,jplot], color='red'
   endif

;   outmodflux = spflux_read_kurucz(outloglam-alog10(1.+zpeak), outdispimg,iselect=ibest)

   return,bestflux
end
;------------------------------------------------------------------------------
function spflux_goodfiber, pixmask
   qgood = ((pixmask AND pixelmask_bits('NOPLUG')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('BADTRACE')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('BADFLAT')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('BADARC')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('MANYBADCOLUMNS')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('NEARWHOPPER')) EQ 0) $
       AND ((pixmask AND pixelmask_bits('MANYREJECTED')) EQ 0)
   return, qgood
end

;------------------------------------------------------------------------------
function spflux_bspline, loglam, mratio, mrativar, outmask=outmask, $
 everyn=everyn, airmass=airmass

   isort = sort(loglam)
   nord = 3

   ; Choose the break points using the EVERYN option, but masking
   ; out more pixels near stellar features just when selecting them.
   mask1 = 1 - spflux_masklines(loglam, hwidth=12.e-4, /stellar)
   ii = where(mrativar[isort] GT 0 AND mask1[isort] EQ 1)
   bkpt = 0
   fullbkpt = bspline_bkpts(loglam[isort[ii]], everyn=everyn, $
    bkpt=bkpt, nord=nord)

   outmask1 = 0
   if (keyword_set(airmass)) then begin
      x2 = airmass[isort]
   endif
   sset = bspline_iterfit(loglam[isort], mratio[isort], $
    invvar=mrativar[isort], lower=3, upper=3, fullbkpt=fullbkpt, $
    maxrej=ceil(0.05*n_elements(indx)), outmask=outmask1, nord=nord, $
    x2=x2, npoly=2*keyword_set(airmass), requiren=(everyn-1)>1)
   if (max(sset.coeff) EQ 0) then $
    message, 'B-spline fit failed!!'
   if (arg_present(outmask)) then begin
      outmask = bytarr(size(loglam,/dimens))
      outmask[isort] = outmask1
   endif

   return, sset
end
;------------------------------------------------------------------------------

function typingmodule, objflux,loglam,objivar,dispimg,sfd_ebv,psfmag,kindx=kindx,modflux=modflux,plottitle=plottitle


   npix = (size(objflux,/dimen))[0]
   ;----------
   ; For each star, find the best-fit model.

;   !p.multi = [0,1,2]
   modflux = 0 * objflux
   ; Find the best-fit model -- evaluated for each exposure [NPIX,NEXP]
   thismodel = spflux_bestmodel(loglam, objflux, $
       objivar, dispimg, kindx=kindx, plottitle=plottitle)

   ; Also evaluate this model over a big wavelength range [3006,10960] Ang.
   tmploglam = 3.4780d0 + lindgen(5620) * 1.d-4
   tmpdispimg = 0 * tmploglam + 1.0 ; arbitrarily select this resolution
   tmpdispimg = interpol(dispimg,loglam,tmploglam)
   tmpflux = spflux_read_kurucz(tmploglam-alog10(1+kindx.z), tmpdispimg, $
       iselect=kindx.imodel)

   ; The returned models are redshifted, but not fluxed or
   ; reddened.  Do that now...  we compare data vs. model reddened.
   extcurve1 = ext_odonnell(10.^loglam, 3.1)
   thismodel = thismodel * 10.^(-extcurve1 * 3.1 * sfd_ebv / 2.5)
   extcurve2 = ext_odonnell(10.^tmploglam, 3.1)
   tmpflux = tmpflux * 10.^(-extcurve2 * 3.1 * sfd_ebv / 2.5)

   ; Now integrate the apparent magnitude for this spectrum,
   ; The units of FTHRU are such that m = -2.5*alog10(FTHRU) + (48.6-2.5*17)
   ; Note that these computed magnitudes, THISMAG, should be equivalent
   ; to THISINDX.MAG in the case of no reddening.
   wavevec = 10.d0^tmploglam
   flambda2fnu = wavevec^2 / 2.99792e18
   fthru = filter_thru(tmpflux * flambda2fnu, waveimg=wavevec, /toair)
   thismag = -2.5 * alog10(fthru) - (48.6-2.5*17)

   ; Compute SCALEFAC = (plugmap flux) / (flux of the model spectrum)

   scalefac = 10.^((thismag[2]-psfmag[2])/2.5)
   thismodel = thismodel * scalefac
   kindx.mag = reform(thismag)+psfmag[2]-thismag[2]

   modflux = thismodel
   splog, prelog=''
   !p.multi = 0
   return,0
end
;------------------------------------------------------------------------------

function matchbundleratio2, bundleratio, ratioivar, reffiber, distarr, cenwave=cenwave,dwave=dwave,predictratio=predictratio,usemask=usemask,weight=weight,ivarweight=ivarweight,lambda=lambda
;common com_fluxbundle, cube,xarr,yarr,wavearr
common com_fluxbundle2, psfplane,xarr,wavearr
    
    if n_elements(cenwave) ne n_elements(dwave) then message,'Error: cenwave and dwave must have the same number of elements.'
    if n_elements(weight) ne n_elements(lambda) then message,'Error: lambda and weight must have the same number of elements.'
    if (size(weight,/dimen))[0] ne (size(ivarweight,/dimen))[0] then message,"Error: ivarweight's first dimension must match weight."
;    stop
    xpos = interpol(findgen(n_elements(xarr)),xarr,distarr)
    nfiber = (size(distarr,/dimen))[0]
    if nfiber ne (size(ivarweight,/dimen))[2] then message,"Error: ivarweight's third dimension must match nfiber."

    nwave = n_elements(cenwave)
    n_wavearr = n_elements(wavearr)

    siz = size(bundleratio,/dimen)
    nfile = siz[1]
    if siz[2] ne nfiber then message,'Error: the 3rd dimension of bundleratio does not match nfiber.'
    if siz[0] ne nwave then message,'Error: the 1st dimension of bundleratio does not match nwave.'

    covfn=fltarr(nwave,nfile,nfiber)
    for i=0,nfiber-1 do begin
      covfn_wavearr = interpolate(psfplane,xpos[i,*],indgen(n_wavearr)) 
      covfn_allwave = interpol(covfn_wavearr,wavearr,lambda[*,*,i])
      weightedflux = covfn_allwave*weight[*,*,i]*ivarweight[*,*,i]
      for ifile=0,nfile-1 do begin
        for j=0,nwave-1 do begin
           lam1=cenwave[j]-dwave[j]/2.
	   lam2=cenwave[j]+dwave[j]/2.
	   ind = where(lambda[*,ifile,i] ge lam1 and lambda[*,ifile,i] lt lam2 and ivarweight[*,ifile,i] gt 0, ngood)
	   if ngood eq 0 then continue
           covfn[j,ifile,i] = total(weightedflux[ind,ifile])/total(ivarweight[ind,ifile,i])
	endfor
      endfor
    endfor
    refpredict = rebin(covfn[*,*,reffiber],nwave,nfile,nfiber)
    predictratio=covfn/(refpredict + (refpredict eq 0))
    useid = where(indgen(nfiber) ne reffiber and usemask,nuse)
    diff = (bundleratio[*,*,useid]-predictratio[*,*,useid])*sqrt(ratioivar[*,*,useid] < (1/0.003)^2)
    ind_valid = where(ratioivar[*,*,useid] gt 0,nvalid)
    chi2 = total(diff[ind_valid]*diff[ind_valid])/(nvalid-1)
    return,chi2
end

function mcmcbundleratio,fluxratio,fluxratioivar,offsetx0=offsetx0,offsety0=offsety0,fiberoffx=fiberoffx,fiberoffy=fiberoffy,nfiber=nfiber,n_wavearr=n_wavearr,reffiber=reffiber,centerwave=centerwave,dwave=dwave,innermask=innermask,weight=weight,ivarweight=ivarweight,lambda=lambda,scale=scale,rotation=rotation,xbest=xbest,ybest=ybest,nstep=nstep,debug=debug,errorreport=errorreport
      stepsize = [0.06,2,0.08] ; [positional move radius, rotation in degree, scale of ADR]
      if NOT keyword_set(nstep) then nstep =1000L
      random=randomu(seed,nstep*4)
      rr = random[0:nstep-1]
      theta = random[nstep:2*nstep-1]*!pi*2
      xstep = rr*cos(theta)
      ystep = rr*sin(theta)
      anglestep = random[2*nstep:3*nstep-1]-0.5;*stepsize[1]*!dtor
      scalestep = random[3*nstep:4*nstep-1]-0.5;*stepsize[2])
      offsetx = offsetx0
      offsety = offsety0
      dxarr = fiberoffx#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsetx
      dyarr = fiberoffy#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsety
      distarr = sqrt(dxarr*dxarr+dyarr*dyarr)
      oldchi2 = matchbundleratio2(fluxratio,fluxratioivar,reffiber,distarr,cenwave=centerwave,dwave=dwave,usemask=innermask,weight=weight,ivarweight=ivarweight,lambda=lambda)
      xposarr = fltarr(nstep+1)
      yposarr = fltarr(nstep+1)
      scalearr = fltarr(nstep+1)
      rotationarr = fltarr(nstep+1)
      chi2arr = fltarr(nstep+1)
      stepnum = intarr(nstep+1)
      chi2arr[0] = oldchi2
      newscale = 1.0
      newrotation = 0.0
      newx = fiberoffx
      newy = fiberoffy 
      npoint=0
      badstep_ct=0
      bestcounter=0
      bestchi2 = oldchi2
      for i=0,nstep-1 do begin  
;         if i mod 1000 eq 0 then print,i
         tmpscale = (newscale*10^(scalestep[i]*stepsize[2]) < 1.2) > 0.8
         tmprotation = (newrotation + anglestep[i]*stepsize[1]*!dtor < 0.17 ) > (-0.17)
         rotmatrix = [[cos(tmprotation),-sin(tmprotation)],[sin(tmprotation),cos(tmprotation)]]
         offsetx=offsetx0*tmpscale
         offsety=offsety0*tmpscale
         xy2= [[offsetx],[offsety]]#rotmatrix
         offsetx = xy2[*,0]
         offsety = xy2[*,1]
         dx = xstep[i]*stepsize[0]
         dy = ystep[i]*stepsize[0]
         dxarr = (newx+dx)#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsetx
         dyarr = (newy+dy)#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsety
         distarr = sqrt(dxarr*dxarr+dyarr*dyarr)
         
         newchi2 = matchbundleratio2(fluxratio,fluxratioivar,reffiber,distarr,cenwave=centerwave,dwave=dwave,usemask=innermask,weight=weight,ivarweight=ivarweight,lambda=lambda)
         if newchi2 lt oldchi2 or randomu(seed,1) lt exp(oldchi2-newchi2) then begin
             npoint += 1
       	     xposarr[npoint] = xposarr[npoint-1]+dx
  	     yposarr[npoint] = yposarr[npoint-1]+dy
             scalearr[npoint] = tmpscale
             rotationarr[npoint] = tmprotation
	     chi2arr[npoint] = newchi2
	     oldchi2 = newchi2
             newx = newx+dx
	     newy = newy+dy 
             newscale = tmpscale
             newrotation = tmprotation
	     if newchi2 lt bestchi2*0.95 then begin
	       oldbestchi2 = bestchi2
	       bestchi2 = newchi2
	       bestcounter = 0
	     endif else bestcounter += 1
	     stepnum[npoint]=i
         endif else begin
             badstep_ct +=1
             if badstep_ct gt 300 then begin
                stepsize = stepsize*0.6
                if stepsize[0] lt 0.003 and stepsize[1] lt 0.25 and stepsize[2] lt 0.004 then break
                badstep_ct=0
             endif 
         endelse

	 if bestcounter gt 300 then break
      endfor
;      print,'Iteration stopped on step ',i
;      print,'Number of bad steps:',badstep_ct
      chi2arr = chi2arr[0:npoint]
      xposarr = xposarr[0:npoint]
      yposarr = yposarr[0:npoint]
      scalearr = scalearr[0:npoint]
      rotationarr = rotationarr[0:npoint]
      tmp = min(chi2arr,ind)
      xbest = xposarr[ind]
      ybest = yposarr[ind]
      scale = scalearr[ind]
      rotation = rotationarr[ind]
      fchi2 = chi2arr[ind]
      if keyword_set(errorreport) then begin
        q=where(chi2arr lt 2*min(chi2arr),nq)
        xerr = stddev(xposarr[q])
        yerr = stddev(yposarr[q])
        print,'MCMC X uncertainty:',xerr
        print,'MCMC Y uncertainty:',yerr
      endif
      erase
      multiplot,[0,1,5,0,0]
      plot,stepnum,xposarr,ytitle='X POS'
      multiplot
      plot,stepnum,yposarr,ytitle='Y POS'
      multiplot
      plot,stepnum,scalearr,ytitle='Scale'
      multiplot
      plot,stepnum,rotationarr,ytitle='Rotation'
      multiplot
      plot,stepnum,chi2arr,ytitle='Chi Square'
      multiplot,[1,1],/default
      if keyword_set(debug) then stop
      return,fchi2
end


function modelthrupt, distarr, wavelength
common com_fluxbundle2, psfplane,xarr,wavearr

   nfiber = (size(distarr,/dimen))[0]
   if (size(distarr,/dimen))[1] ne n_elements(wavelength) then message,"Error:distarr's second dimension must agree with the size of wavelength."
   xpos = interpol(findgen(n_elements(xarr)),xarr,distarr)
   wavepos = interpol(findgen(n_elements(wavearr)),wavearr,wavelength)
   nwave = n_elements(wavepos)
   covfn=fltarr(nwave,nfiber)
   for i=0,nfiber-1 do begin
     covfn[*,i] = interpolate(psfplane,xpos[i,*],wavepos)
   endfor
   return,covfn
end

; innames should be an nexp x ncamera array of names
; pointing to the corresponding mgSFrame files.
pro mdrp_fluxcal, innames, outnames=outnames, bcalib=bcalib,rcalib=rcalib,docams=docams, visual=visual, plotname=plotname, centerwave=centerwave, dwave=dwave,staroffx=staroffx,staroffy=staroffy, combinedir=combinedir,writefile=writefile,writename=writename

common com_fluxbundle2, psfplane,xarr,wavearr
common plotcolors, black, red, green, blue, cyan, magenta, yellow, white,grey

    manga_redux_dir = getenv('MANGA_SPECTRO_REDUX')
    drpver = ~keyword_set(drpver) ? mangadrp_version(/simple) : strtrim(drpver,2)


  ; Set the output names
  outnames=strarr(n_elements(innames))
  for i=0,n_elements(innames)-1 do begin
    outnames[i]=ml_strreplace(innames[i],'mgSFrame','mgFFrame')
    ; Make sure outname name doesn't end in '.gz'
    outnames[i]=ml_ensurenogz(outnames[i])
  endfor      

  nfile = n_elements(innames)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Quality control checks.  If fails return from function
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  for ifile=0, nfile-1 do begin
    infile=lookforgzip(innames[ifile])
    ml_mgframeread, infile, hdr=hdr
    ; Check that input file read properly
    if (~keyword_set(hdr)) then begin
      splog,strcompress('File '+infile+' not found, skipping!')
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
      plotname=fileandpath(ml_psfilename(outnames[0]),path=path)
      path=path+'qa/'
      plotname=path+plotname
      if file_test(path,/directory) eq 0 then spawn, '\mkdir -p '+path
      plotname=ml_strreplace(plotname,'mgFFrame-b','mgFFrame-spec')
    endelse
    splog,'Printing to ',plotname
    dfpsplot, plotname, /color
  endelse

    lat = 32+46/60.+49/3600.d
    lon = -(105.+49/60.+13/3600.d)

    plateid = lonarr(nfile)
    mjd = lonarr(nfile)
    camname = strarr(nfile)
    camid = strarr(nfile)
    expnum = lonarr(nfile)
    spectroid = lonarr(nfile)
    npixarr = lonarr(nfile)
    dcnra = fltarr(nfile)
    dcndec = fltarr(nfile)
    jd = dblarr(nfile)
    ha = fltarr(nfile)

    for ifile=0, nfile-1 do begin
       infile=lookforgzip(innames[ifile])
       ml_mgframeread, infile, hdr=hdr
       plateid[ifile] = strtrim(sxpar(hdr, 'PLATEID'),2)
       mjd[ifile] = strtrim(sxpar(hdr, 'MJD'),2)
       camname[ifile] = strtrim(sxpar(hdr, 'CAMERAS'),2)
       spectroid[ifile] = strmid(camname[ifile],1,1)
       camid[ifile] = strmid(camname[ifile],0,1)
       expnum[ifile] = sxpar(hdr, 'EXPOSURE')
      ; hdr1 = headfits(innames[ifile],exten='FLUX')
       npixarr[ifile] = sxpar(hdr, 'NAXIS1')
       dcnra[ifile] = sxpar(hdr,'MGDRA')
       dcndec[ifile] = sxpar(hdr,'MGDDEC')
       jd[ifile] = (sxpar(hdr,'TAI-BEG')+sxpar(hdr,'EXPTIME')/2.)/86400.+2400000.5
    endfor
    maxmjd = max(mjd)

    uexp = uniq(expnum, sort(expnum))
    exposures = expnum[uexp]
    nexp = n_elements(uexp)

    for iexp=0,nexp-1 do begin
       ind = where(expnum eq exposures[iexp])
;       platestr = string(plateid[0],format='(i0.0)')
;       mjdstr = string(mjd[0],format='(i0.0)')
;       cogimg_dir = djs_filepath('', root_dir=ml_getenv('MANGA_SPECTRO_REDUX'), subdir=[drpver,platestr,mjdstr]) 
;       cogimg_file= file_search(cogimg_dir,'cogimg-'+mjdstr+'-'+string(exposures[iexp],format='(i8.8)')+'.fits*',count=ct)
;       if ct eq 0 then $
       coaddgimg,plateid[ind[0]],mjd[ind[0]],exposures[iexp]
    endfor

    ml_mgframeread,innames[0],slitmap=slitmap

    plugged=where(slitmap.plugstatus eq 0)
    eq2hor,slitmap[plugged[0]].cenra*(dblarr(nfile)+1),slitmap[plugged[0]].cendec*(dblarr(nfile)+1),jd,alt,az,ha,lat=lat,lon=lon,altitude=2788.

    usp = uniq(spectroid,sort(spectroid))
    if n_elements(usp) ne 1 then message, 'Run mdrp_fluxcal per spectrograph'

   ;Truncate the slitmap to one of the spectrograph
    slitmap = slitmap[where(slitmap.spectrographid eq spectroid[0])]

    iphoto = where(strtrim(slitmap.targettype,2) eq 'standard' and slitmap.plugstatus eq 0,nphoto)

    if (nphoto EQ 0) then begin 
        splog, 'WARNING: NO standard stars in this spectrograph for flux calibration'
        return
    endif

    npix = max(npixarr)
    loglam = dblarr(npix, nfile, nphoto)
    objflux = fltarr(npix, nfile, nphoto)
    objivar = fltarr(npix, nfile, nphoto)
    dispimg = fltarr(npix, nfile, nphoto)
;    airmass = fltarr(npix, nfile, nphoto1+nphoto2)

    for ifile=0,nfile-1 do begin
       ml_mgframeread,innames[ifile],iphoto,objflux=objflux1,objivar=objivar1,mask=mask1,loglam=loglam1,$
         wset=wset1, dispset=dispset1,dispimg=dispimg1

      ; Make a map of the size of each pixel in delta-(log10-Angstroms).
      ; Re-normalize the flux to ADU/(dloglam).
      ; Re-normalize the dispersion from /(raw pixel) to /(new pixel).
       correct_dlam,objflux1,objivar1,wset1,dlam=dloglam
       correct_dlam,dispimg1,0,wset1,dlam=dloglam,/inverse

       objivar1 = objivar1 * spflux_goodfiber(mask1)

       loglam[0:npixarr[ifile]-1,ifile,*] = loglam1
       ;it wont do having a tail of zeros in the wavelength so add some dummy values
       if (npix GT npixarr[ifile]) then begin
         dllam=loglam1[npixarr[ifile]-1,*]-loglam1[npixarr[ifile]-2,*]
         for j=0, nphoto-1 do $
            loglam[npixarr[ifile]:*,ifile,j] = loglam1[npixarr[ifile]-1,j]+dllam[0,j]*(1+findgen(npix-npixarr[ifile]))
       endif
       objflux[0:npixarr[ifile]-1,ifile,*] = objflux1
       ;hopefully the inverse variance of 0 of non-filled objects will indicate the uselessness
       ; of the extra
       objivar[0:npixarr[ifile]-1,ifile,*] = skymask(objivar1, mask1, mask1)
       dispimg[0:npixarr[ifile]-1,ifile,*] = dispimg1
    endfor 

    thru=mrdfits(getenv('MANGADRP_DIR')+'/etc/meanthrupt.fits',1)
    init_thru=interpol(thru.thrupt,thru.loglam,loglam)

    frlplug = slitmap[iphoto].frlplug
    uuplug = uniq(frlplug,sort(frlplug))
    starplug = frlplug[uuplug]
    nstar = n_elements(starplug)

    euler,slitmap[iphoto[uuplug]].ra,slitmap[iphoto[uuplug]].dec,l,b,1
    sfd_ebv = dust_getval(l,b,/interp)
    starttime=systime(/sec)
    for istar=0,nstar-1 do begin
       starttime_istar=systime(/sec)
       print,'Processing standard star in Ferrule '+string(starplug[istar],format='(i0.0)')
       bind = where(slitmap[iphoto].frlplug eq starplug[istar],nfiber)
       ss = sort(slitmap[iphoto[bind]].fnum)
       bind = bind[ss]   
    ; slitmap[iphoto[bind]] gives the fibers in the bundle sorted in order of fnum
    ; bind also correspond to the row numbers in objflux and objivar.
    ; as those correspond to part of the original frame as read in with indices iphoto.

       platescale = 3.62730/60. ; (mm/arcsec) as given in Gunn et al. 2006
       fiberx = slitmap[iphoto[bind]].xpmm/platescale
       fibery = slitmap[iphoto[bind]].ypmm/platescale
;        innermask = (sqrt(fiberx*fiberx+fibery*fibery) lt 3) and bmap.gbu ge 0

;  Now identify which fibers belong to the central 7 fiber group in the bundle. This is for cases when we have standard star in a science IFU bundle.
       xmedian = median(slitmap[iphoto[bind]].xpmm)/platescale
       ymedian = median(slitmap[iphoto[bind]].ypmm)/platescale

       dist = sqrt((fiberx-xmedian)^2 + (fibery-ymedian)^2)
       tmp = min(dist,center)
       innermask = (sqrt((fiberx-fiberx[center])^2+(fibery-fibery[center])^2) lt 3) and (slitmap[iphoto[bind]].gbu eq 0)
       innerfiber = where(innermask,ninner)

       bind = bind[innerfiber]
       nfiber = ninner

       lam0 = 3500.

       if NOT (keyword_set(centerwave) and keyword_set(dwave)) then begin
          wave1 = [3500.,4000.,4500.,5000.,5500.,6500.,7500.,9000.]
          wave2 = [4000.,4500.,5000.,5500.,6500.,7500.,9000.,10500.]
          nwave = n_elements(wave1)
          centerwave = (wave1+wave2)/2.
          dwave = wave2-wave1
       endif else begin
          if n_elements(centerwave) ne n_elements(dwave) then message,'centerwave and dwave must have the same number of elements.'
          nwave = n_elements(centerwave)
          wave1 = centerwave-dwave/2.
          wave2 = centerwave+dwave/2.
       endelse

       totnumer = total(total(objflux[*,*,bind]*objivar[*,*,bind],1),1)
       totdenom = total(total(objivar[*,*,bind],1),1)
       fiberflux = totnumer/totdenom
       tmp = max(fiberflux,reffiber)

       thisbundle = starplug[istar]
       thisfiber = iphoto[bind[reffiber]]
       plottitle = 'PLATE=' + string(plateid[0], format='(i4.4)') $
                   + ' MJD=' + string(maxmjd, format='(i5.5)') $
                   + ' Spectro-Photo Star' $
                   + ' Bundle '+ strtrim(thisbundle,2) $
                   + ' Fiber ' + strtrim(thisfiber,2)

       tmp = typingmodule(objflux[*,*,bind[reffiber]],loglam[*,*,bind[reffiber]],objivar[*,*,bind[reffiber]],dispimg[*,*,bind[reffiber]],sfd_ebv[istar],slitmap[iphoto[bind[reffiber]]].psfmag,kindx=thisindx,modflux=modflux,plottitle=plottitle)
;     modflux returend from typingmodule has the dimensions of [npix,nfile]
;     we duplicate it for each fiber in the bundle.
       modflux = rebin(modflux,npix,nfile,nfiber)

; We will generate a correction curve for each camera/spectrograph/exposure. 
; The two cameras from the same spectrograph and the same exposure should share a common positional offset and a common DAR model, but they will provide separate data points and get separate correction vector.

;   Build fluxarr which will contain the binned flux for each fiber and each file
;     blue camera contain 3500-6500A
;     red  camera contain 5800-10500A. 

       fluxarr = fltarr(nwave,nfile,nfiber)
       ivararr = fltarr(nwave,nfile,nfiber)
       for ifile=0,nfile-1 do begin
  	   for ifiber=0,nfiber-1 do begin
              for j=0,nwave-1 do begin
                lam1 = alog10(wave1[j])
                lam2 = alog10(wave2[j])
                indwave=where(loglam[*,ifile,bind[ifiber]] gt lam1 and loglam[*,ifile,bind[ifiber]] lt lam2 and objivar[*,ifile,bind[ifiber]] gt 0,ngood)
                if ngood eq 0 then continue
                fluxarr[j,ifile,ifiber] = total(objflux[indwave,ifile,bind[ifiber]]*objivar[indwave,ifile,bind[ifiber]])/total(objivar[indwave,ifile,bind[ifiber]])
                ivararr[j,ifile,ifiber] = total(objivar[indwave,ifile,bind[ifiber]])
  	      endfor
	   endfor
       endfor

;   Recompute fiberflux, this time do not add different camera files/exposures together.
       fiberflux = total(objflux[*,*,bind]*objivar[*,*,bind],1)/total(objivar[*,*,bind],1)

       refwave = 5300.

       for iexp =0,nexp-1 do begin

          indexp = where(expnum eq exposures[iexp],nexpfile)
          psf=psfparams(plateid[indexp[0]],mjd[indexp[0]],expnum[indexp[0]])

          fluxratio = fltarr(nwave,nexpfile,nfiber)
          fluxratioerr = fltarr(nwave,nexpfile,nfiber)
	  tmp = max(total(fiberflux[indexp,*],1),reffiber)
	  refflux = rebin(fluxarr[*,indexp,reffiber],nwave,nexpfile,nfiber)
	  refivar = rebin(ivararr[*,indexp,reffiber],nwave,nexpfile,nfiber)
	  fluxratio = fluxarr[*,indexp,*]/(refflux + (refflux eq 0))
	  fluxratioerr = abs(fluxratio)*sqrt(1/(fluxarr[*,indexp,*]^2*ivararr[*,indexp,*])+1/(refflux^2*refivar))
	  fluxratioivar = 1/(fluxratioerr*fluxratioerr)
	  q = where(finite(fluxratioerr) eq 0,n_NaN)
          if n_NaN gt 0 then fluxratioivar[q] = 0.0
	  
          if NOT keyword_set(staroffx) then staroffx = 0.0
          if NOT keyword_set(staroffy) then staroffy = 0.0
          fiberoffx = slitmap[iphoto[bind]].xpmm/platescale+dcnra[indexp[0]]-staroffx
          fiberoffy = slitmap[iphoto[bind]].ypmm/platescale+dcndec[indexp[0]]-staroffy

          tmp = min(abs(centerwave-refwave),i_refwave)
	  i_cam = where(camid[indexp] eq 'b')
;	  guess_x = total(fiberoffx*reform(fluxarr[i_refwave,indexp[i_cam],*]))/total(fluxarr[i_refwave,indexp[i_cam],*])
;	  guess_y = total(fiberoffy*reform(fluxarr[i_refwave,indexp[i_cam],*]))/total(fluxarr[i_refwave,indexp[i_cam],*])
;	  fiberoffx = fiberoffx-guess_x
;	  fiberoffy = fiberoffy-guess_y

          nsize = 7
          factorarr=1+(findgen(nsize)-nsize/2)*0.05
          xcen = fltarr(nsize) 
          ycen = fltarr(nsize) 
          scale = fltarr(nsize)
          rotation = fltarr(nsize)
          fchi2 = fltarr(nsize)
          for is=0,nsize-1 do begin
            factor = factorarr[is]
            params=psf.params
            params[2:3] = params[2:3]*factor
            params[0:1] = params[0:1]/factor^2
            fcpsf1d,params,slitmap[iphoto[bind[0]]].xfocal,slitmap[iphoto[bind[0]]].yfocal,psfplane,wavearr=wavearr,xarr=xarr,magnify=2,/resetwave
	    n_wavearr = n_elements(wavearr)

            if is eq 0 then begin 
              offsetx0 = fltarr(n_wavearr)
	      offsety0 = fltarr(n_wavearr)
	      for jj=0,n_wavearr-1 do begin
  	         status=ml_dar((ha[indexp[0]]-(slitmap[iphoto[bind[0]]].ra-slitmap[iphoto[bind[0]]].cenra))/15.,slitmap[iphoto[bind[0]]].dec,wavearr[jj],parangle=parangle,offsetX=offsetx_tmp,offsetY=offsety_tmp,waveref=refwave,raobj=slitmap[iphoto[bind[0]]].ra,racen=slitmap[iphoto[bind[0]]].cenra,deccen=slitmap[iphoto[bind[0]]].cendec,/distort)
	         offsetx0[jj] = offsetx_tmp
		 offsety0[jj] = offsety_tmp
	      endfor
	    endif

            fchi2[is]=mcmcbundleratio(fluxratio,fluxratioivar,$
	              offsetx0=offsetx0, offsety0=offsety0,$
		      fiberoffx=fiberoffx, fiberoffy=fiberoffy,$
		      nfiber=nfiber,n_wavearr=n_wavearr,reffiber=reffiber,$
                      centerwave=centerwave,dwave=dwave,innermask=intarr(nfiber)+1,$
                      weight=modflux[*,indexp,*]*init_thru[*,indexp,bind],ivarweight=objivar[*,indexp,bind],$
                      lambda=10^loglam[*,indexp,bind],xbest=xbest,ybest=ybest,scale=scalebest,$
                      rotation=rotationbest,nstep=200)

            xcen[is] = xbest
            ycen[is] = ybest
            scale[is] = scalebest
            rotation[is] = rotationbest
          endfor

          if nsize ge 5 then begin
            res = svdfit(factorarr,fchi2,3)
            factor = -res[1]/(2*res[2])

            params=psf.params
            params[2:3] = params[2:3]*factor
            params[0:1] = params[0:1]/factor^2
       
            fcpsf1d,params,slitmap[iphoto[bind[0]]].xfocal,slitmap[iphoto[bind[0]]].yfocal,psfplane,wavearr=wavearr,xarr=xarr,magnify=2,/resetwave
            bestchi2=mcmcbundleratio(fluxratio,fluxratioivar,$
               offsetx0=offsetx0,offsety0=offsety0,fiberoffx=fiberoffx,$
               fiberoffy=fiberoffy,nfiber=nfiber,n_wavearr=n_wavearr,$
               reffiber=reffiber,centerwave=centerwave,$
               dwave=dwave,innermask=intarr(nfiber)+1,$
	       weight=modflux[*,indexp,*]*init_thru[*,indexp,bind],$
               ivarweight=objivar[*,indexp,bind],lambda=10^loglam[*,indexp,bind],$
	       xbest=xbest,ybest=ybest,scale=scalebest,rotation=rotationbest,nstep=2000,/errorreport)
            if bestchi2 gt min(fchi2) then begin
              bestchi2 = min(fchi2,ind_s)
              factor = factorarr[ind_s]
              params=psf.params
              params[2:3] = params[2:3]*factor
              params[0:1] = params[0:1]/factor^2
              fcpsf1d,params,slitmap[iphoto[bind[0]]].xfocal,slitmap[iphoto[bind[0]]].yfocal,psfplane,wavearr=wavearr,xarr=xarr,magnify=2,/resetwave
              scalebest = scale[ind_s]
              rotationbest=rotation[ind_s]
              xbest= xcen[ind_s]
              ybest= ycen[ind_s]
            endif 
          endif else begin
            bestchi2 = min(fchi2,ind_s)
            factor = factorarr[ind_s]
            scalebest = scale[ind_s]
            rotationbest=rotation[ind_s]
            xbest= xcen[ind_s]
            ybest= ycen[ind_s]
            params=psf.params
            params[2:3] = params[2:3]*factor
            params[0:1] = params[0:1]/factor^2
            fcpsf1d,params,slitmap[iphoto[bind[0]]].xfocal,slitmap[iphoto[bind[0]]].yfocal,psfplane,wavearr=wavearr,xarr=xarr,magnify=2,/resetwave
          endelse  

          offsetx = offsetx0*scalebest
          offsety = offsety0*scalebest
          rotmatrix = [[cos(rotationbest),-sin(rotationbest)],[sin(rotationbest),cos(rotationbest)]]
          xy2 = [[offsetx],[offsety]]#rotmatrix
          offsetx = xy2[*,0]
          offsety = xy2[*,1]
          dxarr = (fiberoffx+xbest)#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsetx
          dyarr = (fiberoffy+ybest)#(fltarr(n_wavearr)+1)-(fltarr(nfiber)+1)#offsety
          distarr = sqrt(dxarr*dxarr+dyarr*dyarr)

          newchi = matchbundleratio2(fluxratio,fluxratioivar,reffiber,distarr,cenwave=centerwave,dwave=dwave,usemask=intarr(nfiber)+1,weight=modflux[*,indexp,*]*init_thru[*,indexp,bind],ivarweight=objivar[*,indexp,bind],lambda=10^loglam[*,indexp,bind],predict=predict)

          example = [0,3,6]
;          psplot,'fluxratio_example_'+string(starplug[istar],format='(i0.0)')+'.ps',/square,/color
;          psset,/multi
          !p.multi=[0,3,3]
          for iwave=0,n_elements(example)-1 do begin &$
	     if example[iwave] le 3 then icam = where(camid[indexp] eq 'b') else icam = where(camid[indexp] eq 'r')
             plotratiofit,fiberoffx,fiberoffy,string(slitmap[iphoto[bind]].fnum,format='(i0.0)'),string(fluxratio[example[iwave],icam,*]*100,format='(f0.2)'),string(predict[example[iwave],icam,*]*100,format='(f0.2)'),string(fluxratioerr[example[iwave],icam,*]*100,format='(f3.1)'),titles=['Fiber Num','Data','Best-fit Model'],color=1,charsize=1. &$
             oplot,-xbest+offsetx,-ybest+offsety,color=254 &$
             useless = min(abs(wavearr-centerwave[example[iwave]]),id_w) &$
             plots,-xbest+offsetx[id_w],-ybest+offsety[id_w],ps=1,symsize=2,color=254 &$
             legend,'!7k!6='+string(centerwave[example[iwave]],format='(i0.0)')+'A',box=0,charsize=1.5 &$
          endfor
	  !p.multi=0
;          psclose

; We shall use all the fibers that have significant flux to derive the final correction. But we need to figure out how to do this without resampling the native-wavelength grid. (keep these code for now to be reactivated later)
if 0 then begin
          snratio = objflux[*,indexp,bind]*sqrt(objivar[*,indexp,bind])
          snratiobin = rebin(snratio,12,nexpfile,nfiber)
          snratiobinr = snratiobin/rebin(snratiobin[*,*,reffiber],12,nexpfile,nfiber)
          select = (snratiobinr gt (sqrt(2.)-1)) and finite(snratiobinr)
          selectbig = rebin(select,npix,nexpfile,nfiber,/sample)
endif

;          data = total(objflux[*,indexp,bind]*selectbig,3)
;          var = total((1/objivar[*,indexp,bind])*selectbig,3)
;          ivar = 1/var
	  data = objflux[*,indexp,bind[reffiber]]
	  ivar = objivar[*,indexp,bind[reffiber]]
          

;          ivarprep = objivar[*,indexp,bind]
;          ind = where(selectbig eq 0) 
;          ivarprep[ind] = 1.0
;          ivarproduct = product(ivarprep,3)
;          ivar[where(ivarproduct eq 0)] = 0.0
 
          covfn = modelthrupt(distarr,wavearr)
          covfnbig = fltarr(npix,nexpfile,nfiber)
          for i=0,nfiber-1 do covfnbig[*,*,i] = interpol(covfn[*,i],wavearr,10^loglam[*,indexp,bind[i]],/spline)
;          totcovfn = total(covfnbig*selectbig,3)
       
          model = modflux[*,*,reffiber]*covfnbig[*,*,reffiber]
          mratio = data/(model+ (model eq 0))
          mrativar = ivar*model^2
          mrativar = mrativar *(1- spflux_masklines(loglam[*,indexp,bind[reffiber]],/stellar))
          
;	  outloglam= dblarr(npix,nexpfile)
;	  fluxcorr = dblarr(npix,nexpfile)
	  outloglam = reform(loglam[*,indexp,bind[reffiber]])
;	  for ifile=0,nexpfile-1 do begin
;	     CASE camid[indexp[ifile]] OF 
;                'b': BEGIN
;		      bsset=spflux_bspline(loglam[*,indexp[ifile],bind[reffiber]],mratio[*,ifile],mrativar[*,ifile],everyn=10)
;                      fluxcorr[*,ifile]=bspline_valu(loglam[*,indexp[ifile],bind[reffiber]],bsset)
;		      outloglam[*,ifile]=loglam[*,indexp[ifile],bind[reffiber]]
;		     END
;	        'r': BEGIN
;		      rsset=spflux_bspline(loglam[*,indexp[ifile],bind[reffiber]],mratio[*,ifile],mrativar[*,ifile],everyn=1.5)
;		      fluxcorr[*,ifile]=bspline_valu(loglam[*,indexp[ifile],bind[reffiber]],rsset)
;		      outloglam[*,ifile]=loglam[*,indexp[ifile],bind[reffiber]]
;		     END
;             ENDCASE
;	  endfor
	  ss = sort(camid[indexp]) ; sort it so blue camera is first and red camera is second.
	  mratio=mratio[*,ss]
	  mrativar=mrativar[*,ss]

          result={plate:0L,mjd:0L,exposure:0L,dcnra:0.0,dcndec:0.0,frlplug:0L,ra:0.0d,dec:0.0d,$
	          psfmag:fltarr(5),factor:0.0,xbest:0.0,ybest:0.0,$
		  fiberoffx:fiberoffx,fiberoffy:fiberoffy,fnum:slitmap[iphoto[bind]].fnum,$
		  offsetx:offsetx,offsety:offsety,$
                  scalebest:0.0,rotationbest:0.0,fluxratio_chi2:0.0,modelspec:model,modflux:modflux[*,indexp,reffiber],$
                  data:objflux[*,indexp,bind[reffiber]],mratio:mratio,mrativar:mrativar,$
	          loglam:outloglam,fluxratio:fluxratio,fluxratioerr:fluxratioerr,predictratio:predict}
          result = create_struct(result,thisindx)
          if n_elements(totresult) eq 0 then begin
    	     totresult = replicate(result,nexp*nstar)
          endif
	  result.plate = plateid[0]
	  result.mjd = maxmjd
          result.exposure=exposures[iexp]
	  result.dcnra = dcnra[indexp[0]]
	  result.dcndec = dcndec[indexp[0]]
          result.frlplug=starplug[istar]
	  result.ra = slitmap[iphoto[bind[0]]].ra
	  result.dec = slitmap[iphoto[bind[0]]].dec
	  result.psfmag = slitmap[iphoto[bind[0]]].psfmag
          result.factor=factor
          result.xbest=xbest;-guess_x
          result.ybest=ybest;-guess_y
          result.scalebest=scalebest
          result.rotationbest=rotationbest
          result.fluxratio_chi2 =bestchi2
          totresult[iexp*nstar+istar]=result
        endfor
	print,'Time used to process all exposures of star:',systime(/sec)-starttime_istar
     endfor
      print,'Total time used to process all exposures of all stars:',systime(/sec)-starttime
      mjdstr=string(mjd[0],format='(i0.0)')
      if NOT keyword_set(writename) then writename=''
      if keyword_set(writefile) then begin
        for iexp=0,nexp-1 do begin
           filename='fluxcal-sp'+string(spectroid[0],format='(i0.0)')+'-'+string(exposures[iexp],format='(i8.8)')+writename+'.fits'
	   ind = where(totresult.exposure eq exposures[iexp])
           mwrfits,totresult[ind],filename,/create
        endfor
      endif



; Final loop over input frames
; Read in the frames, apply the calibrations, and write
; out the calibrated mgFFrame files.  Note that these will
; have units of flux/pixel rather than flux per some uniform
; wavelength range.
      for iexp=0,nexp-1 do begin
        indexp = where(expnum eq exposures[iexp],nexpfile)
        inds = where(totresult.exposure eq exposures[iexp])        
	; select only good stars --- rejecting those with S/N less than 1/3 of the median and those with Chi2 greater than 3 times the median
	gstar = where(totresult[inds].sn_median gt $
	  median(totresult[inds].sn_median,/even)/3. and $
	  totresult[inds].linechi2/totresult[inds].linedof lt $
          median(totresult[inds].linechi2/totresult[inds].linedof,/even)*3,ngstar) 
        fmratio   = fltarr(npix, 2, ngstar)
	fmrativar = fltarr(npix, 2, ngstar)

        ; take out a lower-order polynomial relative to the mean correction vector
	; This will hopefully average out the systematics introduced in imperfect 
	; model match. 
        flatarrblue = spflux_mratio_flatten(totresult[inds[gstar]].loglam[*,0],$
	             totresult[inds[gstar]].mratio[*,0],totresult[inds[gstar]].mrativar[*,0],pres=pres_b)
	flatarrred = spflux_mratio_flatten(totresult[inds[gstar]].loglam[*,1],$
	             totresult[inds[gstar]].mratio[*,1],totresult[inds[gstar]].mrativar[*,1],pres=pres_r)
        fmratio[*,0,*] = totresult[inds[gstar]].mratio[*,0]/flatarrblue
	fmrativar[*,0,*] = totresult[inds[gstar]].mrativar[*,0]*flatarrblue^2

	fmratio[*,1,*] = totresult[inds[gstar]].mratio[*,1]/flatarrred
	fmrativar[*,1,*] = totresult[inds[gstar]].mrativar[*,1]*flatarrred^2

	sset_b = spflux_bspline(totresult[inds[gstar]].loglam[*,0],$
	            fmratio[*,0,*],fmrativar[*,0,*],everyn=10*ngstar)
        sset_r = spflux_bspline(totresult[inds[gstar]].loglam[*,1],$
	            fmratio[*,1,*],fmrativar[*,1,*],everyn=1.5*ngstar)
        !p.multi=[0,1,2]
        djs_plot,10^totresult[inds[gstar]].loglam[*,0],fmratio[*,0,*],ps=1,symsize=0.1,xran=[3500,6500],xst=1,xtitle='Wavelength (A)',ytitle='Correction Vector'
	djs_oplot,10^totresult[0].loglam[*,0],bspline_valu(totresult[0].loglam[*,0],sset_b),color='red'
        djs_plot,10^totresult[inds[gstar]].loglam[*,1],fmratio[*,1,*],ps=1,symsize=0.1,xran=[6000,10500],xst=1,xtitle='Wavelength (A)',ytitle='Correction Vector'
	djs_oplot,10^totresult[0].loglam[*,1],bspline_valu(totresult[0].loglam[*,1],sset_r),color='red'
	!p.multi=0
        
        for ifile=0,nexpfile-1 do begin
          ml_mgframeread,innames[indexp[ifile]],loglam=loglam1, wset=wset1, objflux=tempflux, $
            objivar=tempivar, mask=tempmask, dispset=dispset1, dispimg=tempdisp, slitmap=slitfull, ximg=tempximg, $
            slithdr=slithdr, superflat=superflat1, skyflux=tempsky, hdr=thishdr, adderr=adderr
          ; Reading slithdr from the mgFrame is buggy and doesn't work so well
          ; when trying to write back out again.  Read it again from mangacore
          plugname=yanny_par(slithdr,'plugmap')
          slitname=ml_strreplace(plugname,'plPlugMapM','slitmap')
          temp=ml_readslitmap(slitname,hdr=slithdr)

          ; Make a map of the size of each pixel in 1e-4 (log10-Angstroms),
          ; and re-normalize the flux, ivar, sky to electrons/(dloglam)
          correct_dlam, tempflux, tempivar, wset1
          correct_dlam, tempsky, 0, wset1

          splog,'Applying flux calibration vectors to mgSFrame inputs.'

          CASE camid[indexp[ifile]] OF 
             'b' : BEGIN
                     bcalib = float(bspline_valu(loglam1,sset_b))
                     minval = 0.05*median(bcalib)
                     ; Apply calibration vector to flux, inverse variance, and sky extensions
                     ; This will put it in units of flux/Angstrom

                     mdrp_divideflat, tempflux, invvar=tempivar, bcalib, minval=minval
                     mdrp_divideflat, tempsky, bcalib, minval=minval
                     ; Flag pixels with bad flux calibration factor
                     tempmask = tempmask OR ((bcalib LE minval) * pixelmask_bits('BADFLUXFACTOR'))
                   END
             'r' : BEGIN
                     rcalib = float(bspline_valu(loglam1,sset_r))
                     minval = 0.05*median(rcalib)
                     ; Apply calibration vector to flux, inverse variance, and sky extensions
                     ; This will put it in units of flux/Angstrom
                     mdrp_divideflat, tempflux, invvar=tempivar, rcalib, minval=minval
                     mdrp_divideflat, tempsky, rcalib, minval=minval
                     ; Flag pixels with bad flux calibration factor
                     tempmask = tempmask OR ((rcalib LE minval) * pixelmask_bits('BADFLUXFACTOR'))
                   END
          ENDCASE

          ; Add flux information keywords
          sxaddpar, thishdr, 'BSCALE', 1., ' Flux unit scaling'
          sxaddpar, thishdr, 'BZERO', 0., ' Flux zeropoint'
          sxaddpar, thishdr, 'BUNIT', '1E-17 erg/s/cm^2/Ang/fiber', ' Flux units are per fiber and per Angstrom'

          ; Add quality control keywords
          drpqual=long(fxpar(thishdr,'DRPQUAL'))
          ; Add new flagging cases here
          ; None yet!
          ; Write back to header
          sxaddpar, thishdr, 'DRPQUAL', drpqual, 'DRP quality bitmask'

          ; Write sky-subtracted spectra to disk following the MaNGA data model
          ml_mwrfits, dummyext, outnames[indexp[ifile]],hdr=thishdr, /create ; Blank ext 0 with full header
          ml_mwrfits, tempflux, outnames[indexp[ifile]], extname='FLUX'       ;sky subtracted flux
          ml_mwrfits, tempivar, outnames[indexp[ifile]], extname='IVAR'   ; sky subtracted inverse variance
          ml_mwrfits, tempmask, outnames[indexp[ifile]], extname='MASK'     ; final pixel mask
          ml_mwrfits, wset1, outnames[indexp[ifile]], extname='WSET'          ;trace-set for wavelength sol; wset
          ml_mwrfits, dispset1, outnames[indexp[ifile]], extname='DISPSET'    ;trace-set for dispersion sol
          ml_mwrfits, slitfull, outnames[indexp[ifile]], hdr=slithdr, extname='SLITMAP' ;slitmap used
          ml_mwrfits, tempximg, outnames[indexp[ifile]], extname='XPOS'                        ;x pos on CCD
          ml_mwrfits, superflat1, outnames[indexp[ifile]], extname='SUPERFLAT'              ;superflat vector from quartz lamps
          ml_mwrfits, tempsky, outnames[indexp[ifile]], extname='SKY'                       ;sky flux

          ; gzip output file
          spawn, ['gzip', '-f', outnames[indexp[ifile]]], /noshell

       endfor
     endfor

  ; Close out plots
  if ~keyword_set(visual) then begin
    dfpsclose
    ; Convert to PDF and remove old PS file
    ;stop
    spawn, strcompress('ps2pdf '+plotname+' '+ml_strreplace(plotname,'.ps','.pdf'))
    spawn, strcompress('rm -f '+plotname)
    !p.multi=0
    !p=origp
    !x=origx
    !y=origy
  endif
        
  ; Take out the trash
  heap_gc



end
