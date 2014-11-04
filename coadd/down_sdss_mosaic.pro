pro down_sdss_mosaic, ra, dec, fov, filter, download = download, jpeg = jpeg  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Download image from SDSS archive to build large mosaic  
;;; written by Song Huang 03.28.2012 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
on_error, 2
compile_opt idl2

if N_params() lt 4 then begin 
    print,  'Syntax - Down_SDSS_Mosaic, RA, DEC, FOV, Filter, /Download, /JPEG'
    print,  'RA, DEC : RA in degree; both need to be in decimal format'
    print,  'FOV     : in square degree, need to be smaller than 1.0'
    print,  'Filter  : among u, g, r, i, z'
    print,  '/Download: download images instead just building the image list'
    print,  '/JPEG    : download a cutoff jpeg image of this field'
    return
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
wget = '/usr/bin/wget'
bunzip2 = '/bin/bunzip2'
ra  = float( ra ) 
dec = float( dec ) 
fov = float( fov ) 
filter = string( filter ) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if ( ( ra lt 0.0 ) or ( ra gt 360.0 ) or ( dec gt 90.0 ) or ( dec lt -90.0 ) ) $
    then begin 
    message, 'Check the format of RA or DEC !!' 
endif else begin 
    ra_string = strcompress( string( ra ), /remove_all ) 
    dec_string = strcompress( string( dec ), /remove_all )
endelse
if ( ( fov lt 0.0 ) or ( fov gt 1.0 ) ) then begin 
    message, 'The FOV should between 0.0 and 1.0 ' 
endif else begin 
    fov_string = strcompress( string( fov, format='(F4.2)' ),/ remove_all )
endelse 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
band_array = [ 'u', 'g', 'r', 'i', 'z' ]
if ( where( band_array eq filter ) eq -1 ) then begin 
    message, 'The filter should be among u, g, r, i, z'
endif else begin 
    filter = strcompress( filter, /remove_all )
endelse 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; download the mosaic tool result page from SDSS DR8 SAS 
webpage = '"http://data.sdss3.org/mosaics.html?ra=' + ra_string + '&dec=' + $ 
    dec_string + '&size=' + fov_string + '&mosaic_bands=ugriz&submit=Submit"'
output  = 'temp.html'
print, webpage
;;spawn, wget + '   ' + webpage + ' -O  ' + output
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; read in the output html file 
if NOT file_test( output ) then begin 
    message, 'Somehting wrong with the wget !!'
endif 
lines = file_lines( output ) 
source = strarr( lines ) 
openr, 10, output 
readf, 10, source 
close, 10
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
index = where( strpos( source, "this.className='highlight'" ) ne -1 ) 
nimage = n_elements( index ) 
if ( nimage eq -1 ) then begin 
    print, 'Sorry! No useful image !'
    print, 'Perhaps this coordinate is not covered by SDSS '
endif else begin 
    print, 'There are ' + string( nimage ) + ' images in this region'
    run       = strarr( nimage ) 
    camcol    = strarr( nimage ) 
    field     = strarr( nimage ) 
    ra_center = strarr( nimage )
    dec_center= strarr( nimage )
    add_fits  = strarr( nimage )
    name_fits = strarr( nimage ) 

    for i = 0, ( nimage - 1 ), 1 do begin 
        line0 = index[i]
        ;; line for the RUN 
        line1 = ( line0 + 1 )
        run_line = source[line1] 
        temp = strsplit( run_line, '<>', /extract ) 
        run[i] = strcompress( temp[2] )
        ;; line for the CAMCOL 
        line2 = ( line0 + 2 )
        cam_line = source[line2]
        temp = strsplit( cam_line, '<>', /extract ) 
        camcol[i] = strcompress( temp[2] ) 
        ;; line for the FIELD
        line3 = ( line0 + 3 )
        fld_line = source[line3]
        temp = strsplit( fld_line, '<>', /extract ) 
        field[i] = strcompress( temp[2] ) 
        ;; line for the CENTRAL RA 
        line4 = ( line0 + 4 )
        cra_line = source[line4]
        temp = strsplit( cra_line, '<>', /extract ) 
        ra_center[i] = strcompress( temp[2] ) 
        ;; line for the CENTRAL DEC
        line5 = ( line0 + 5 )
        cde_line = source[line5]
        temp = strsplit( cde_line, '<>', /extract ) 
        dec_center[i] = strcompress( temp[2] ) 

        ;; define the address for the fits image 
        add_fits[i] = '"http://data.sdss3.org/returnIms/fits?run=' + run[i] + $
            '&camcol=' + camcol[i] + '&field=' + field[i] + '&filter=' + $
            filter + '"'
        name_fits[i] = 'frame-' + filter + '-' + run[i] + '-' + camcol[i] + $
            '-' + field[i] + '.fits.bz2'
        ;; download the image if the keyword download is used 
        if keyword_set( download ) then begin 
            spawn, wget + '  ' + add_fits[i] + '  -O  ' + name_fits[i]
            if file_test( name_fits[i] ) then begin 
                spawn, bunzip2 + '  ' + name_fits[i] 
            endif
        endif

    endfor

    ;; save a summary file for this request 
    log = ra_string + '_' + dec_string + '_' + fov_string + '_' + filter + $
        '.log'
    openw, 20, log 
    printf, 20, '##  RA CENTER     = ' + ra_string
    printf, 20, '##  DEC CENTER    = ' + dec_string
    printf, 20, '##  Field of View = ' + fov_string 
    printf, 20, '##  Filter        = ' + filter
    printf, 20, '##  RUN       CAMCOL      FIELD      RA      DEC '   
    for j = 0, ( nimage - 1 ), 1 do  begin 
        printf, 20, run[j] + '    ' + camcol[j] + '    ' + field[j] + '   ' + $ 
            ra_center[j] + '    ' + dec_center[j] + ' '
    endfor
    close, 20 

    ;; save a download script 
    sh = ra_string + '_' + dec_string + '_' + fov_string + '_' + filter + $
        '.sh'
    openw, 30, sh, width = 1000 
    printf, 30, '#!/bin/sh'
    for k = 0, ( nimage - 1 ), 1 do begin 
        printf, 30, 'wget ' + add_fits[k] + ' -O ' + name_fits[k]
    endfor
    close, 30 

    ;; if keyword "jpeg" is set, download a cutoff jpeg image 
    if keyword_set( jpeg ) then begin 
        size_jpg = long( ( fov * 60.0 * 60.0 ) / 0.396 )
        size_str = strcompress( string( size_jpg ), /remove_all )
        add_jpeg = '"http://skyservice.pha.jhu.edu/DR8/ImgCutout/' + $
            'getjpeg.aspx?ra=' + ra_string + '&dec=' + dec_string + $ 
            '&scale=0.396&width=' + size_str + '&height=' + size_str  $
            + '&opt=GL&query="'
        name_jpg = ra_string + '_' + dec_string + '_' + fov_string + '_' + $
            filter + '.jpg'
        spawn, 'wget ' + add_jpeg + ' -O ' + name_jpg
    endif
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
