pro cgs_galrot, image, xc, yc, r80, pa, img_rot, mask_name = mask_name, $ 
    savefits = savefits, fits_name = fits_name

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; rotate the galaxy image 
;; written by SHUANG, 06.26.2012 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
on_error, 2
compile_opt idl2

if N_params() lt 6 then begin 
    print,  'Syntax - cgs_galrot, image, xc, yc, r80, pa, img_rot, /mask, '
    print,  '         mask_rot=mask_rot, /savefits' 
    return
endif

pix = 0.259 
pix_area = ( pix^2.0 )
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; prepare the input parameters
if ( ( size( image))[0] ne 2 ) then begin 
    message, 'The input image should be a two dimension image ! '
endif else begin 
    temp = size( image )
    naxis1 = ( temp[1] )
    naxis2 = ( temp[2] )
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xc = long( xc ) 
yc = long( yc ) 
if ( ( xc le 1 ) or ( xc ge naxis1 ) or ( yc le 1 ) or $ 
    ( yc ge naxis2 ) ) then begin 
    message, 'Something worng with the location of center !' 
endif 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
r80 = float( r80 ) 
if ( ( r80 le 1 ) or ( r80 ge ( naxis1 / 2.0 ) ) or ( r80 ge ( naxis2 / 2.0 ) $
    ) ) then begin 
    message, 'Something wrong with the R80 !' 
endif 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pa = float( pa ) 
if ( ( pa gt 360.0 ) or ( pa lt -180.0 ) ) then begin 
    message, '!!! Something wrong with the PA !!! ' 
endif
if ( pa lt 0 ) then begin 
    print, '!!! PA should be between 0 and 180 degree !!!' 
    pa = ( pa + 180.0 ) 
    print, ' PA : ' + string( pa - 180.0 ) + ' ===> ' + string( pa ) 
endif else begin 
    if ( pa gt 180.0 ) then begin  
        print, '!!! PA should be between 0 and 180 degree !!!' 
        pa = ( pa - 180.0 ) 
        print, ' PA : ' + string( pa + 180.0 ) + ' ===> ' + string( pa ) 
    endif 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; mask image 
if keyword_set( mask ) then begin 
    temp = size( mask ) 
    if ( ( temp[0] ne 2 ) or ( temp[1] ne naxis1 ) or $ 
        ( temp[2] ne naxis2 ) ) then begin 
        message, 'Something worng with the mask image !!!'
    endif else begin 
        mask[ where( mask gt 0 ) ] = 1 
    endelse 
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; rotate the image 
;angle = pa
angle = ( pa - 90.0 )
;; int = 0 : neares neighbor 
;;       1 : bilinear 
;;       2 : cubic
img_rot = ROT( image, angle, 1.0, xc, yc, /INTERP, MISSING=-999.9 )
if keyword_set( savefits ) then begin 
    if keyword_set( fits_name ) then begin 
        outfits = strcompress( fits_name, /remove_all )
    endif else begin 
        outfits = 'rotate.fits' 
    endelse
    mwrfits, img_rot, outfits, /create 
endif
img_rot = img_rot 

;; if mask is provided, rotate the mask as well 
if keyword_set( mask_name ) then begin 
    msk_rot = ROT( mask_name, angle, 1.0, xc, yc, /INTERP, MISSING=1 ) 
    msk_rot[ where( msk_rot ge 0.99 ) ] = 1
    msk_rot[ where( msk_rot lt 0.99 ) ] = 0
    if keyword_set( savefits ) then begin 
        temp = strsplit( outfits, '.', /extract ) 
        mask_file = strcompress( temp[0], /remove_all ) + '_mask.fit' 
        mwrfits, msk_rot, mask_file, /create 
    endif
    mask_rot = mask_rot
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
