; + 
; 
; NAME:
;              OBJECT_MASK_BUILDING 
;
; PURPOSE:
;              Build a mask image for surface brightness measurements
;
; USAGE:
;    
;
; ARGUMENTS:
;    list_file: 
;
;                                                            
; KEYWORDS:
;
; OUTPUT:
;
; AUTHOR:
;             Song Huang
;
; HISTORY:
;
; TODO LIST: 
;-
; CATEGORY:    HS_SDSS
;------------------------------------------------------------------------------

pro object_mask_building, image, region = region, bg_size = bg_size, $
    thre_l = thre_l, thre_h = thre_h, lt_size = lt_size, ht_size = ht_size, $
    magzpt = magzpt, pixel = pixel, star_gauss = star_gauss, $
    gauss_cut = gauss_cut, star_class = star_class, $
    expand_factor = expand_factor, seeing = seeing 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

on_error, 2
compile_opt idl2

if N_params() lt 1 then begin 
    print,  'Syntax - Object_Mask_Building, Image, [Region], [bg_size],  '
    print,  '        [thre_l], [thre_h], [lt_size], [ht_size], [magzpt], '
    print,  '        [pixel], [star_gauss], [gauss_cut], [star_class],   '
    print,  '        [expand_factor]   '
    return
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Parameter List
if NOT keyword_set( bg_size ) then begin 
    bg_size = string( 1024 )
endif else begin 
    bg_size = string( bg_size )
endelse
if NOT keyword_set( thre_l ) then begin 
    thre_l = string( 2.5 )
endif else begin 
    thre_l = string( thre_l ) 
endelse
if NOT keyword_set( thre_h ) then begin 
    thre_h = string( 2.5 )
endif else begin 
    thre_h = string( thre_h ) 
endelse
if NOT keyword_set( lt_size ) then begin 
    lt_size = string( 128 )
endif else begin 
    lt_size = string( lt_size ) 
endelse
if NOT keyword_set( ht_size ) then begin 
    ht_size = string( 8 )
endif else begin 
    ht_size = string( ht_size ) 
endelse
if NOT keyword_set( magzpt ) then begin 
    magzpt = string( 23.5 )
endif else begin 
    magzpt = string( magzpt ) 
endelse
if NOT keyword_set( pixel ) then begin 
    pixel = string( 0.259 )
endif else begin 
    pixel = string( pixel ) 
endelse
if NOT keyword_set( star_gauss ) then begin 
    star_gauss = 6
endif else begin 
    star_gauss = float( star_gauss )  
endelse
lt_gauss     = star_gauss
ht_gauss     = star_gauss
if NOT keyword_set( gauss_cut ) then begin 
    gauss_cut = 0.6
endif else begin 
    gauss_cut = float( gauss_cut ) 
endelse
lt_gauss_cut = gauss_cut
ht_gauss_cut = gauss_cut
if NOT keyword_set( star_class ) then begin 
    star_class = 0.7
endif else begin 
    star_class = float( star_class )
endelse
non_star_class = ( star_class - 0.2 ) 
if NOT keyword_set( expand_factor ) then begin 
    expand_factor = 1.1
endif else begin 
    expand_factor = float( expand_factor )
endelse
if NOT keyword_set( seeing ) then begin 
    seeing = 1.0
endif else begin 
    seeing = float( seeing ) 
endelse
min_nrad_ht    = 8.0 
min_nrad_lt    = 10.0
;; constant 
bg_thick= string( 24 )
lt_thick= string( 24 )
ht_thick= string( 24 )
filter_name    = 'default.conv'
;filter_name    = 'tophat_1.5_3x3.conv'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
image = strcompress( image, /remove_all ) 
if NOT file_test( image ) then begin 
    message, 'Can not find the original image: ' + image 
endif else begin 
    print, '#################################################################'
    print, 'IMAGE : ' + image 
    img = mrdfits( image, 0, header ) 
    temp = strsplit( image, '.', /extract )
    name_string = temp[0]
    temp = size( img, /dimension ) 
    naxis1 = temp[0] 
    naxis2 = temp[1]
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 0. Make an index array indicates the regions that are avoid from masking  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
blank = long( img - img )
if keyword_set( region ) then begin 
    region = strcompress( region, /remove_all ) 
    if NOT file_test( region ) then begin 
        message, 'Can not find the region file !! '
    endif else begin 
        input_image = blank 
        reg_file = region
        region_combine, input_image, reg_file, output_image, index = -99 
        void_mask = output_image
        write_jpeg, 'void.jpg', void_mask, quality = 100
    endelse
endif else begin 
    void_mask = blank 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 1. Run SExtractor, get the background and subtract it out 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
template = 'shuang.sex' 
if NOT file_test( template ) then begin 
    message, 'Can not find shuang.sex file ! '
endif else begin 
    num_lines = file_lines( template )
    para = strarr( num_lines ) 
    openr, 10, template
    readf, 10, para 
    close, 10 
endelse
;; Change the catalog name 
index = where( strpos( para, 'CATALOG_NAME' ) ne -1 ) 
cat_name = 'temp.fits'
para[ index ] = 'CATALOG_NAME ' + cat_name 
;; Change the background size 
index = where( strpos( para, 'BACK_SIZE' ) ne -1 ) 
para[ index ] = 'BACK_SIZE   ' + bg_size 
;; Change the background thick 
index = where( strpos( para, 'BACKPHOTO_THICK' ) ne -1 ) 
para[ index ] = 'BACKPHOTO_THICK ' + bg_thick
;; Change the background image 
back_file = name_string + '_background.fit' 
back_sex  = name_string + '_background.sex' 
index = where( strpos( para, 'CHECKIMAGE_TYPE' ) ne -1 ) 
para[ index ] = 'CHECKIMAGE_TYPE BACKGROUND'
index = where( strpos( para, 'CHECKIMAGE_NAME' ) ne -1 ) 
para[ index ] = 'CHECKIMAGE_NAME ' + back_file 
;; Change the magnitude zeropoint
index = where( strpos( para, 'MAG_ZEROPOINT' ) ne -1 ) 
para[ index ] = 'MAG_ZEROPOINT ' + magzpt
;; Change the pixel scale
index = where( strpos( para, 'PIXEL_SCALE' ) ne -1 ) 
para[ index ] = 'PIXEL_SCALE ' + pixel
;; Save the new Sex file 
openw, 20, back_sex, width = 600 
for i = 0, ( num_lines - 1 ), 1 do begin 
    printf, 20, para[i] 
endfor 
close, 20 
;; Run SExtractor to get the background image 
spawn, 'sex  ' + image + ' -c ' + back_sex 
;; Subtract the background 
if NOT file_test( back_file ) then begin 
    print, 'Can not find the background image ! '
    print, 'No subtraction is applied !!' 
    newimage = image 
endif else begin 
    img_bkg = mrdfits( back_file, 0 )
    img_new = img - img_bkg 
    new_file = name_string + '_backfree.fit' 
    mwrfits, img_new, new_file, header, /create  
    newimage = new_file
endelse 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 2. Run SExtractor, get the objects detection under low threshold
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
template = 'shuang.sex' 
if NOT file_test( template ) then begin 
    message, 'Can not find shuang.sex file ! '
endif else begin 
    num_lines = file_lines( template )
    para = strarr( num_lines ) 
    openr, 10, template
    readf, 10, para 
    close, 10 
endelse
;; Change the catalog name 
index = where( strpos( para, 'CATALOG_NAME' ) ne -1 ) 
lt_cat = name_string + '_lt.fits'
para[ index ] = 'CATALOG_NAME ' + lt_cat
;; Change the detection filter  
index = where( strpos( para, 'FILTER_NAME' ) ne -1 ) 
para[ index ] = 'FILTER_NAME ' + filter_name 
;; Change the catalog type 
index = where( strpos( para, 'CATALOG_TYPE' ) ne -1 ) 
para[ index ] = 'CATALOG_TYPE FITS_LDAC'
;; Change the detection threshold 
index = where( strpos( para, 'DETECT_THRESH' ) ne -1 ) 
para[ index ] = 'DETECT_THRESH ' + thre_l
;; Change the background size 
index = where( strpos( para, 'BACK_SIZE' ) ne -1 ) 
para[ index ] = 'BACK_SIZE ' + lt_size
;; Change the background thick 
index = where( strpos( para, 'BACKPHOTO_THICK' ) ne -1 ) 
para[ index ] = 'BACKPHOTO_THICK ' + lt_thick
;; Change the background image 
ltseg_file = name_string + '_lt_seg.fit' 
ltapr_file = name_string + '_lt_apr.fit' 
lt_sex  = name_string + '_lt.sex' 
index = where( strpos( para, 'CHECKIMAGE_TYPE' ) ne -1 ) 
para[ index ] = 'CHECKIMAGE_TYPE SEGMENTATION,APERTURES'
index = where( strpos( para, 'CHECKIMAGE_NAME' ) ne -1 ) 
para[ index ] = 'CHECKIMAGE_NAME ' + ltseg_file + '  ' + ltapr_file 
;; Change the magnitude zeropoint
index = where( strpos( para, 'MAG_ZEROPOINT' ) ne -1 ) 
para[ index ] = 'MAG_ZEROPOINT ' + magzpt
;; Change the pixel scale
index = where( strpos( para, 'PIXEL_SCALE' ) ne -1 ) 
para[ index ] = 'PIXEL_SCALE ' + pixel
;; Change the seeing fwhm
index = where( strpos( para, 'SEEING_FWHM' ) ne -1 ) 
para[ index ] = 'SEEING_FWHM ' + string( seeing )
;; Save the new Sex file 
openw, 20, lt_sex, width = 600 
for i = 0, ( num_lines - 1 ), 1 do begin 
    printf, 20, para[i] 
endfor 
close, 20 
;; Run SExtractor using the low detection threshold 
;spawn, 'sex  ' + newimage + ' -c ' + lt_sex 
spawn, 'sex  ' + image + ' -c ' + lt_sex 
;; Read in the segmentation image 
if NOT file_test( ltseg_file ) then begin 
    message, 'Can not find the low threshold output segmentation image ! '
endif else begin 
    img_lt = mrdfits( ltseg_file, 0 ) 
endelse 
;; Read in the detection catalog  
if NOT file_test( lt_cat ) then begin 
    message, 'Can not find the low threshold detection catalog ! '
endif else begin 
    cat_lt = mrdfits( lt_cat, 2 ) 
endelse 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 3. Run SExtractor, get the objects detection under high threshold
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
template = 'shuang.sex' 
if NOT file_test( template ) then begin 
    message, 'Can not find shuang.sex file ! '
endif else begin 
    num_lines = file_lines( template )
    para = strarr( num_lines ) 
    openr, 10, template
    readf, 10, para 
    close, 10 
endelse
;; Change the catalog name 
index = where( strpos( para, 'CATALOG_NAME' ) ne -1 ) 
ht_cat = name_string + '_ht.fits'
para[ index ] = 'CATALOG_NAME ' + ht_cat
;; Change the detection filter  
index = where( strpos( para, 'FILTER_NAME' ) ne -1 ) 
para[ index ] = 'FILTER_NAME ' + filter_name 
;; Change the catalog type 
index = where( strpos( para, 'CATALOG_TYPE' ) ne -1 ) 
para[ index ] = 'CATALOG_TYPE FITS_LDAC'
;; Change the detection threshold 
index = where( strpos( para, 'DETECT_THRESH' ) ne -1 ) 
para[ index ] = 'DETECT_THRESH ' + thre_h
;; Change the background size 
index = where( strpos( para, 'BACK_SIZE' ) ne -1 ) 
para[ index ] = 'BACK_SIZE  ' + ht_size 
;; Change the background thick 
index = where( strpos( para, 'BACKPHOTO_THICK' ) ne -1 ) 
para[ index ] = 'BACKPHOTO_THICK ' + ht_thick
;; Change the background image 
htseg_file = name_string + '_ht_seg.fit' 
htapr_file = name_string + '_ht_apr.fit' 
ht_sex  = name_string + '_ht.sex' 
index = where( strpos( para, 'CHECKIMAGE_TYPE' ) ne -1 ) 
para[ index ] = 'CHECKIMAGE_TYPE SEGMENTATION,APERTURES'
index = where( strpos( para, 'CHECKIMAGE_NAME' ) ne -1 ) 
para[ index ] = 'CHECKIMAGE_NAME ' + htseg_file + '  ' + htapr_file 
;; Change the magnitude zeropoint
index = where( strpos( para, 'MAG_ZEROPOINT' ) ne -1 ) 
para[ index ] = 'MAG_ZEROPOINT ' + magzpt
;; Change the pixel scale
index = where( strpos( para, 'PIXEL_SCALE' ) ne -1 ) 
para[ index ] = 'PIXEL_SCALE ' + pixel
;; Change the seeing fwhm
index = where( strpos( para, 'SEEING_FWHM' ) ne -1 ) 
para[ index ] = 'SEEING_FWHM ' + string( seeing )
;; Save the new Sex file 
openw, 20, ht_sex, width = 600 
for i = 0, ( num_lines - 1 ), 1 do begin 
    printf, 20, para[i] 
endfor 
close, 20 
;; Run SExtractor using the high detection threshold 
spawn, 'sex  ' + newimage + ' -c ' + ht_sex 
;; Read in the segmentation image 
if NOT file_test( htseg_file ) then begin 
    message, 'Can not find the high threshold output segmentation image ! '
endif else begin 
    img_ht = mrdfits( htseg_file, 0 ) 
endelse 
;; Read in the detection catalog  
if NOT file_test( ht_cat ) then begin 
    message, 'Can not find the high threshold detection catalog ! '
endif else begin 
    cat_ht = mrdfits( ht_cat, 2 ) 
endelse 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 4. Remove the segmentation belong to the void mask region 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; First remove all the segmentation in the void region 
img_ht_new = img_ht
img_lt_new = img_lt 
if keyword_set( region ) then begin 

uniq_ht = img_ht_new[ uniq( img_ht_new, sort( img_ht_new ) ) ]
uniq_lt = img_lt_new[ uniq( img_lt_new, sort( img_lt_new ) ) ]
n_uniq_ht = n_elements( uniq_ht )
n_uniq_lt = n_elements( uniq_lt )
print, 'Number of uniq segmentations on lt : ', n_uniq_lt
print, 'Number of uniq segmentations on ht : ', n_uniq_ht
index_void = where( void_mask ne 0 )
member_void_lt = img_lt_new[ index_void ]
member_void_ht = img_ht_new[ index_void ]
uniq_lt_void = member_void_lt[ uniq( member_void_lt, sort( member_void_lt ) ) ]
uniq_ht_void = member_void_ht[ uniq( member_void_ht, sort( member_void_ht ) ) ]
n_uniq_lt_void = n_elements( uniq_lt_void )
n_uniq_ht_void = n_elements( uniq_ht_void )
print, 'Number of uniq segmentations on lt in void: ', n_uniq_lt_void
print, 'Number of uniq segmentations on ht in void: ', n_uniq_ht_void
print, n_elements( index_void )
print, '### Remove segmentations from low threshold image'
for i = 0, ( n_uniq_ht_void - 1 ), 1 do begin 
    if ( uniq_ht_void[i] ne 0 ) then begin 
        aa = uniq_ht_void[i] 
        img_ht_new[ where( img_ht_new eq aa ) ] = 0 
        print, aa
    endif
endfor
print, '### Remove segmentations from high threshold image'
for i = 0, ( n_uniq_lt_void - 1 ), 1 do begin 
    if ( uniq_lt_void[i] ne 0 ) then begin 
        aa = uniq_lt_void[i] 
        img_lt_new[ where( img_lt_new eq aa ) ] = 0 
        print, aa
    endif
endfor
;; Find the stellar-like objects from the catalog 
sindex_ht = where( cat_ht.CLASS_STAR ge star_class ) 
sindex_lt = where( cat_lt.CLASS_STAR ge star_class ) 
star_ht = cat_ht[ sindex_ht ] 
star_lt = cat_lt[ sindex_lt ] 
n_star_ht = n_elements( star_ht.CLASS_STAR )
n_star_lt = n_elements( star_lt.CLASS_STAR )
;; Add back the star on the high threshold segmentation image
for i = 0, ( n_star_ht - 1 ), 1 do begin 
    img_ht_new[ where( img_ht eq ( sindex_ht[i] + 1 ) ) ] = ( sindex_ht[i] + 1 )
endfor
;; Add back the star on the low threshold segmentation image
for i = 0, ( n_star_lt - 1 ), 1 do begin 
    img_lt_new[ where( img_lt eq ( sindex_lt[i] + 1 ) ) ] = ( sindex_lt[i] + 1 )
endfor

endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 5. Convolve the segmentation image with certain gaussian to expand the mask
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
img_lt_conv = img_lt_new 
img_ht_conv = img_ht_new
img_ht_conv[ where( img_ht_conv gt 0 ) ] = 10
img_lt_conv[ where( img_lt_conv gt 0 ) ] = 10
img_ht_conv = filter_image( img_ht_new, FWHM_GAUSSIAN=lt_gauss, /All_Pixels )
img_lt_conv = filter_image( img_lt_new, FWHM_GAUSSIAN=ht_gauss, /All_Pixels )
cutoff_ht = where( img_ht_conv lt ht_gauss_cut )
cutoff_lt = where( img_lt_conv lt lt_gauss_cut )
img_ht_conv[ cutoff_ht ] = 0
img_lt_conv[ cutoff_lt ] = 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 6. Expand the mask for large non-point sources 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print, 'Number of non-star sources on low threshold mask: '
print, n_elements( where( cat_lt.CLASS_STAR le non_star_class ) )
print, 'Number of non-star sources on high threshold mask: '
print, n_elements( where( cat_ht.CLASS_STAR le non_star_class ) )

img_lt_np = img_lt_new 
img_ht_np = img_ht_new
uniq_ht = img_ht_np[ uniq( img_ht_np, sort( img_ht_np ) ) ]
uniq_lt = img_lt_np[ uniq( img_lt_np, sort( img_lt_np ) ) ]
n_uniq_ht = n_elements( uniq_ht )
n_uniq_lt = n_elements( uniq_lt )

print, 'Expand Non-Stellar Object on Low-Threshold Mask:'
for i = 0, ( n_uniq_ht - 1 ), 1 do begin 
    if ( uniq_ht[i] ne 0 ) then begin 
        aa = uniq_ht[i] 
        bb = aa - 1 
        xpeak = cat_ht[ bb ].XPEAK_IMAGE
        ypeak = cat_ht[ bb ].YPEAK_IMAGE
        class = cat_ht[ bb ].CLASS_STAR 
        ellip = cat_ht[ bb ].ELLIPTICITY
        theta = cat_ht[ bb ].THETA_IMAGE
        rad90 = cat_ht[ bb ].FLUX_RADIUS 
        ab    = ( 1.00 - ellip )
        radnew = expand_factor * rad90 

        if ( ( class lt non_star_class ) and ( rad90 gt min_nrad_ht ) ) $
            then begin 

            print, 'SEG: ' + string( aa ) + '  ' + string( rad90 ) + ' ---> ' $
                + string( radnew )
            dist_ellipse, obj_ellip, [naxis1, naxis2], xpeak, ypeak, ab, $
                theta, /double 
            img_ht_np[ where( obj_ellip le radnew ) ] = 10

        endif
    endif
endfor
print, 'Expand Non-Stellar Object on Low-Threshold Mask:'
for i = 0, ( n_uniq_lt - 1 ), 1 do begin 
    if ( uniq_lt[i] ne 0 ) then begin 
        aa = uniq_lt[i] 
        bb = aa - 1 
        xpeak = cat_lt[ bb ].XPEAK_IMAGE
        ypeak = cat_lt[ bb ].YPEAK_IMAGE
        class = cat_lt[ bb ].CLASS_STAR 
        ellip = cat_lt[ bb ].ELLIPTICITY
        theta = cat_lt[ bb ].THETA_IMAGE
        rad90 = cat_lt[ bb ].FLUX_RADIUS 
        ab    = ( 1.00 - ellip )
        radnew = expand_factor * rad90 

        if ( ( class lt non_star_class ) and ( rad90 gt min_nrad_lt ) ) $
            then begin 

            print, 'SEG: ' + string( aa ) + '  ' + string( rad90 ) + ' ---> ' $
                + string( radnew )
            dist_ellipse, obj_ellip, [naxis1, naxis2], xpeak, ypeak, ab, $
                theta, /double 
            img_lt_np[ where( obj_ellip le radnew ) ] = 10

        endif
    endif
endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 7. Combine the high and low threshold segmentation images and save the mask
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
new_mask = img_lt_conv + img_ht_conv + img_lt_np + img_ht_np
new_mask[ where( new_mask gt 0 ) ] = 1
;new_mask = long( new_mask )
;; Read in the original image for building the mask-overlap image 
img_ori = mrdfits( image, 0, header )
min_ori = min( img_ori )
;neg_val = ( 2 * min_ori )
neg_val = !Values.F_Nan
sub_mask = img_ori 
sub_mask[ where( new_mask gt 0 ) ] = neg_val
;; Save the jpeg version of mask and overlap-mask image for preview
jpeg_over = name_string + '_over.jpg'
write_jpeg, jpeg_over, sub_mask, quality=100
;; Save the mask and overlap_mask image into fits file
fits_mask = name_string + '_mask.fits'
fits_over = name_string + '_over.fits'
mwrfits, new_mask, fits_mask, header, /create 
mwrfits, sub_mask, fits_over, header, /create 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
