pro sex_to_ds9, image, bg_size = bg_size, thre = thre, magzpt = magzpt, $
    pixel = pixel, star_class = star_class, rad_thre = rad_thre, $
    color = color, width = width, satu_level = satu_level, seeing = seeing, $
    large_galaxy = large_galaxy, star_only = star_only, $
    all = all, nonstar_only = nonstar_only, brightest = brightest, $
    largest = largest, rad_flux = rad_flux, exclude_satu = exclude_satu

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; run sextractor on image, transfer the segments to DS9 regions  
;;; written by Song Huang 05.31.2012 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
on_error, 2
compile_opt idl2

if N_params() lt 1 then begin 
    print,  'Syntax - sex_to_ds9, Image, [bg_size], [thre], [magzpt], '
    print,  '        [pixel], [star_class], [rad_thre], /large_galaxy  '
    print,  '        /star_only, /all, /galaxy_only, /brightest, /largest'
    print,  '        /rad_flux     '
    return
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Parameter List
if NOT keyword_set( bg_size ) then begin 
    bg_size = string( 1024 )
endif else begin 
    bg_size = string( bg_size )
endelse
if NOT keyword_set( thre ) then begin 
    thre = string( 2.0 )
endif else begin 
    thre = string( thre ) 
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
if NOT keyword_set( star_class ) then begin 
    star_class = 0.85
endif else begin 
    star_class = float( star_class )
endelse
if NOT keyword_set( rad_thre ) then begin 
    rad_thre = 50
endif else begin 
    rad_thre = float( star_class )
endelse
if NOT keyword_set( satu_level ) then begin 
    satu_level = 60000
endif else begin 
    satu_level = long( satu_level ) 
endelse
if NOT keyword_set( width ) then begin 
    width = '2'
endif else begin 
    width = strcompress( string( long( width ) ), /remove_all )
endelse
if NOT keyword_set( color ) then begin 
    color = 'green'
endif else begin 
    color = string( color ) 
endelse
if NOT keyword_set( seeing ) then begin 
    seeing = 1.0
endif else begin 
    seeing = float( seeing ) 
endelse
;; constant 
bg_thick= string( 24 )
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
;; find the saturated pixels 
satu_index = where( img gt satu_level ) 
if ( satu_index[0] ne -1 ) then begin 
    n_satu = n_elements( satu_index ) 
    print, 'Number of Saturated Pixels : ', n_satu
endif else begin 
    n_satu = 0
    print, 'No Saturated Pixels' 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 1. Run SExtractor, get the objects detection
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
cat_name = name_string + '_detec.fits'
para[ index ] = 'CATALOG_NAME ' + cat_name
;; Change the detection filter  
index = where( strpos( para, 'FILTER_NAME' ) ne -1 ) 
para[ index ] = 'FILTER_NAME ' + filter_name 
;; Change the catalog type 
index = where( strpos( para, 'CATALOG_TYPE' ) ne -1 ) 
para[ index ] = 'CATALOG_TYPE FITS_LDAC'
;; Change the detection threshold 
index = where( strpos( para, 'DETECT_THRESH' ) ne -1 ) 
para[ index ] = 'DETECT_THRESH ' + thre
;; Change the background size 
index = where( strpos( para, 'BACK_SIZE' ) ne -1 ) 
para[ index ] = 'BACK_SIZE ' + bg_size
;; Change the background thick 
index = where( strpos( para, 'BACKPHOTO_THICK' ) ne -1 ) 
para[ index ] = 'BACKPHOTO_THICK ' + bg_thick
;; Change the background image 
sex_name  = name_string + '_detec.sex' 
seg_name  = name_string + '_detec_seg.fits'
apr_name  = name_string + '_detec_apr.fits'
index = where( strpos( para, 'CHECKIMAGE_TYPE' ) ne -1 ) 
para[ index ] = 'CHECKIMAGE_TYPE  SEGMENTATION APERTURES '
index = where( strpos( para, 'CHECKIMAGE_NAME' ) ne -1 ) 
para[ index ] = 'CHECKIMAGE_NAME  '  + seg_name + '  ' + apr_name
;; Change the magnitude zeropoint
index = where( strpos( para, 'MAG_ZEROPOINT' ) ne -1 ) 
para[ index ] = 'MAG_ZEROPOINT ' + magzpt
;; Change the pixel scale
index = where( strpos( para, 'PIXEL_SCALE' ) ne -1 ) 
para[ index ] = 'PIXEL_SCALE ' + pixel
;; Change the saturation level
index = where( strpos( para, 'SATUR_LEVEL' ) ne -1 ) 
para[ index ] = 'SATUR_LEVEL ' + string( satu_level )
;; Change the seeing fwhm
index = where( strpos( para, 'SEEING_FWHM' ) ne -1 ) 
para[ index ] = 'SEEING_FWHM ' + string( seeing )
;; Save the new Sex file 
openw, 20, sex_name, width = 600 
for i = 0, ( num_lines - 1 ), 1 do begin 
    printf, 20, para[i] 
endfor 
close, 20 
;; Run SExtractor using the low detection threshold 
;spawn, 'sex  ' + newimage + ' -c ' + lt_sex 
spawn, 'sex  ' + image + ' -c ' + sex_name 
;; Read in the segmentation image 
if NOT file_test( seg_name ) then begin 
    message, 'Can not find the segmentation image ! '
endif else begin 
    seg_img = mrdfits( seg_name, 0 ) 
endelse 
;; Read in the detection catalog  
if NOT file_test( cat_name ) then begin 
    message, 'Can not find the detection catalog ! '
endif else begin 
    cat = mrdfits( cat_name, 2 ) 
endelse 
;; find the saturated segements 
if ( n_satu ne 0 ) then begin 
    satu_seg = seg_img[ uniq( seg_img[ satu_index ], $
        sort( seg_img[ satu_index ] ) ) ]
    n_satu_seg = n_elements( satu_seg ) 
    print, 'Number of saturated segments : ', n_satu_seg
endif else begin 
    n_satu_seg = 0 
    print, 'No Saturated segment !' 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
n_detects = n_elements( cat.NUMBER ) 
regions = strarr( n_detects )
status  = make_array( n_detects, /INTEGER, VALUE=1 ) 
for i = 0, ( n_detects - 1 ), 1 do begin 
    xcen = strcompress( string( cat[i].XPEAK_IMAGE ), /remove_all )
    ycen = strcompress( string( cat[i].YPEAK_IMAGE ), /remove_all ) 
    pa   = strcompress( string( cat[i].THETA_IMAGE ), /remove_all ) 
    radf = cat[i].FLUX_RADIUS 
    ba   = ( 1.0 - cat[i].ELLIPTICITY )
    if keyword_set( rad_flux ) then begin 
        aa = strcompress( string( radf ), /remove_all ) 
        bb = strcompress( string( radf * ba ), /remove_all )
        if ( ( aa le 0.3 ) or ( aa ge ( naxis1 * 1.1 ) ) ) then begin 
            status[i] = 0 
        endif
    endif else begin 
        aa   = strcompress( string( cat[i].A_IMAGE ), /remove_all )
        bb   = strcompress( string( cat[i].B_IMAGE ), /remove_all )
        if ( ( aa le 0.3 ) or ( aa ge ( naxis1 * 1.1 ) ) ) then begin 
            status[i] = 0 
        endif
    endelse
    reg_string = 'ellipse(' + xcen + ',' + ycen + ',' + aa + ',' + bb + $
        ',' + pa + ')' 
    regions[i] = reg_string
    useful = where( status eq 1 )
endfor
;; exclude the saturated objects ? 
if ( keyword_set( exclude_satu ) and ( n_satu_seg ne 0 ) ) then begin 
    status[ satu_seg - 1 ] = 0 
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
header1 = '# Region file format: DS9 version 4.1'
header2 = '# Filename: ' + image 
header3 = 'global dashlist=8 3 font="helvetica 10 ' + $
    'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 ' + $
    'delete=1 include=1 source=1 '
header4 = 'image'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set( all ) then begin
    reg_file  = name_string + '_all.reg' 
    openw, 20, reg_file, width=200 
    printf, 20, header1 
    printf, 20, header2 
    printf, 20, header3 
    printf, 20, header4 
    for j = 0, ( n_detects - 1 ), 1 do  begin 
        if ( status[j] eq 1 ) then begin 
            printf, 20, regions[j] + '  #color=' + color + ' width=' + width
        endif
    endfor 
    close, 20
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set( large_galaxy ) then begin
    header3 = 'global color=red width=3 dashlist=8 3 font="helvetica 10 ' + $
        'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 ' + $
        'delete=1 include=1 source=1 '
    reg_file  = name_string + '_large_galaxy.reg' 
    openw, 20, reg_file, width=200 
    printf, 20, header1 
    printf, 20, header2 
    printf, 20, header3 
    printf, 20, header4 
    for j = 0, ( n_detects - 1 ), 1 do  begin 
        if ( ( status[j] eq 1 ) and ( cat[j].CLASS_STAR le star_class ) $
            and ( cat[j].FLUX_RADIUS ge rad_thre ) ) then begin 
            printf, 20, regions[j]
        endif
    endfor 
    close, 20
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set( star_only ) then begin
    header3 = 'global color=green width=5 dashlist=8 3 font="helvetica 10 ' + $
        'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 ' + $
        'delete=1 include=1 source=1 '
    reg_file  = name_string + '_star.reg' 
    openw, 20, reg_file, width=200 
    printf, 20, header1 
    printf, 20, header2 
    printf, 20, header3 
    printf, 20, header4 
    for j = 0, ( n_detects - 1 ), 1 do  begin 
        if ( ( status[j] eq 1 ) and ( cat[j].CLASS_STAR ge star_class ) ) $ 
            then begin 
            printf, 20, regions[j]
        endif
    endfor 
    close, 20
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set( nonstar_only ) then begin
    header3 = 'global color=yellow width=2 dashlist=8 3 font="helvetica 10 ' + $
        'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 ' + $
        'delete=1 include=1 source=1  '
    reg_file  = name_string + '_nonstar.reg' 
    openw, 20, reg_file, width=200 
    printf, 20, header1 
    printf, 20, header2 
    printf, 20, header3 
    printf, 20, header4 
    for j = 0, ( n_detects - 1 ), 1 do  begin 
        if ( ( status[j] eq 1 ) and ( cat[j].CLASS_STAR le star_class ) ) $ 
            then begin 
            printf, 20, regions[j]
        endif
    endfor 
    close, 20
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set( largest ) then begin
    header3 = 'global color=cyan width=4 dashlist=8 3 font="helvetica 10 ' + $
        'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 ' + $
        'delete=1 include=1 '
    reg_file  = name_string + '_largest.reg' 
    openw, 20, reg_file, width=200 
    printf, 20, header1 
    printf, 20, header2 
    printf, 20, header3 
    printf, 20, header4 
    max_rad = max( cat[ where( status eq 1 ) ].FLUX_RADIUS ) 
    max_index = where( cat.FLUX_RADIUS eq max_rad )
    n_max = n_elements( max_index )
    for j = 0, ( n_max - 1 ), 1 do  begin 
        printf, 20, regions[j]
    endfor 
    close, 20
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set( brightest ) then begin
    header3 = 'global color=magenta width=5 dashlist=8 3 font="helvetica 10 ' + $
        'normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 ' + $
        'delete=1 include=1 source=1 '
    reg_file  = name_string + '_brightest.reg' 
    openw, 20, reg_file, width=200 
    printf, 20, header1 
    printf, 20, header2 
    printf, 20, header3 
    printf, 20, header4 
    max_flux = max( cat[ where( status eq 1 ) ].FLUX_BEST ) 
    max_index = where( cat.FLUX_BEST eq max_flux )
    n_max = n_elements( max_index )
    for j = 0, ( n_max - 1 ), 1 do  begin 
        printf, 20, regions[j]
    endfor 
    close, 20
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
