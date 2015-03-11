;+
; NAME:
;   ml_readpreimaging
;
; PURPOSE:
;   Read data from MaNGA preimaging files.  These are multi-extension
;   files giving SDSS broadband imaging data for target galaxies.
;
; CALLING SEQUENCE:
;   ml_readpreimaging, filename, [hdr=, gimg=, givar=, gpsf=, ghdr=, $
;     rimg=, rivar=, rpsf=, rhdr, iimg=, iivar=, ipsf=, ihdr=,$
;     zimg=, zivar=, zpsf=, zhdr= ]
;
; INPUTS:
;   filename   - Input file name
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   hdr        - Global image header
;   gimg       - g band SDSS image
;   givar      - g band inverse variance image
;   gpsf       - g band PSF model
;   ghdr       - g band image header
;   rimg       - r band SDSS image
;   rivar      - r band inverse variance image
;   rpsf       - r band PSF model
;   rhdr       - g band image header
;   iimg       - i band SDSS image
;   iivar      - i band inverse variance image
;   ipsf       - i band PSF model
;   ihdr       - g band image header
;   zimg       - z band SDSS image
;   zivar      - z band inverse variance image
;   zpsf       - z band PSF model
;   zhdr       - g band image header
;
; COMMENTS:
;   The MaNGA preimaging data are called by the HDU names (EXTNAME)
;   rather than by their number.  The extensions are:
;
;   HDU #0 []:          Empty except for global header
;   HDU #1 [g img]:     g band SDSS image
;   HDU #2 [g ivar]:    g band inverse variance image
;   HDU #3 [g psf]:     g band PSF model
;   HDU #4 [r img]:     r band SDSS image
;   HDU #5 [r ivar]:    r band inverse variance image
;   HDU #6 [r psf]:     r band PSF model
;   HDU #7 [i img]:     i band SDSS image
;   HDU #8 [i ivar]:    i band inverse variance image
;   HDU #9 [i psf]:     i band PSF model
;   HDU #10 [z img]:     z band SDSS image
;   HDU #11 [z ivar]:    z band inverse variance image
;   HDU #12 [z psf]:     z band PSF model
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   lookforgzip()
;   headfits()
;   mrdfits()
;
; REVISION HISTORY:
;   06-Oct-2013  Written by D. Law, Dunlap Institute
;-
;------------------------------------------------------------------------------
pro ml_readpreimaging, filename, hdr=hdr, gimg=gimg, givar=givar, gpsf=gpsf, ghdr=ghdr, $
 rimg=rimg, rivar=rivar, rpsf=rpsf, rhdr=rhdr, iimg=iimg, iivar=iivar, ipsf=ipsf, ihdr=ihdr, $
 zimg=zimg, zivar=zivar, zpsf=zpsf, zhdr=zhdr

  ; If no filename specified, error.
  if (NOT keyword_set(filename)) then begin
    splog,'ERROR: Must specify FILENAME'
    message, 'Must specify FILENAME'
  endif

  ; Check to see if the desired file is gzipped, and change name if so
  thisfile = lookforgzip(filename[0])

  ; Get the image header (whatever is in extension 0)
  if (arg_present(hdr)) then hdr = headfits(thisfile[0])

  ; Get the g-band data
  if ((arg_present(gimg))or(arg_present(ghdr))or(arg_present(givar))or(arg_present(gpsf))) then begin
    gimg=mrdfits(thisfile[0], 'g img',ghdr,/silent)
    givar=mrdfits(thisfile[0], 'g ivar',/silent)
    gpsf=mrdfits(thisfile[0], 'g psf',/silent)
  endif

  ; Get the r-band data
  if ((arg_present(rimg))or(arg_present(rhdr))or(arg_present(rivar))or(arg_present(rpsf))) then begin
    rimg=mrdfits(thisfile[0], 'r img',rhdr,/silent)
    rivar=mrdfits(thisfile[0], 'r ivar',/silent)
    rpsf=mrdfits(thisfile[0], 'r psf',/silent)
  endif

  ; Get the i-band data
  if ((arg_present(iimg))or(arg_present(ihdr))or(arg_present(iivar))or(arg_present(ipsf))) then begin
    iimg=mrdfits(thisfile[0], 'i img',ihdr,/silent)
    iivar=mrdfits(thisfile[0], 'i ivar',/silent)
    ipsf=mrdfits(thisfile[0], 'i psf',/silent)
  endif

  ; Get the z-band data
  if ((arg_present(zimg))or(arg_present(zhdr))or(arg_present(zivar))or(arg_present(zpsf))) then begin
    zimg=mrdfits(thisfile[0], 'z img',zhdr,/silent)
    zivar=mrdfits(thisfile[0], 'z ivar',/silent)
    zpsf=mrdfits(thisfile[0], 'z psf',/silent)
  endif

return
end
