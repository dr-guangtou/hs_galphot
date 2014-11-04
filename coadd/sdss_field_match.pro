pro sdss_field_match, ra, dec, fov, jpeg=jpeg 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; find the sdss images in certain field  
;;; written by Song Huang 06.19.2012 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
on_error, 2
compile_opt idl2

if N_params() lt 4 then begin 
    print,  'Syntax - SDSS_Field_Match, RA, DEC, FOV /JPEG'
    print,  'RA, DEC : RA in degree; both need to be in decimal format'
    print,  'FOV     : in square degree, need to be smaller than 1.0'
    print,  '/JPEG    : download a cutoff jpeg image of this field'
    return
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; constant
cushion = 0.02
;; sdss field summary catalog 
cat = 'sdss_field_sum.fit'
if NOT file_test( cat ) then begin 
    message, 'Can not find the SDSS field summary catalog !!' 
endif else begin 
    sum = mrdfits( cat, 1, head_cat ) 
endelse
nfields = n_elements( sum.run )
status = intarr( nfields )
;; input 
ra  = float( ra ) 
dec = float( dec ) 
fov = float( fov ) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; exam and organize the input
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
;; define the search radius
;; fov is in degree 
fov_arcsec = ( fov * 60.0 * 60.0 )
rad_arcsec = ( fov_arcsec / 2.0 ) + ( cushion * 60.0 * 60.0 ) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for i = 0L, ( nfields - 1L ), 1L do begin 
    ra_min = sum[i].raMin 
    ra_max = sum[i].raMax 
    dec_min = sum[i].decMin 
    dec_max = sum[i].decMax 
    
    





end
