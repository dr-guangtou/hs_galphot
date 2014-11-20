pro cgs_color_map, galaxy, ncomp, type, standard_filter, $
    rebin_factor = rebin_factor, add_noise=add_noise, psf=psf

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Make the Color-Magnitude Map and Color-Color Map for certain model 
;;  Need images and certain files for B, V, and R band models
;;  Be careful !! Especially the location of the files !! 
;;  Written by Song Huang, 05/28/2012 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
on_error, 2
compile_opt idl2

if N_params() lt 4 then begin 
    print,  'Syntax - cgs_color_analysis, galaxy, band1, band2, ncomp, type, '
    print,  '       , /diag_plot, standard_filter = standard_filter, '
    print,  '       , comp_color = comp_color, galaxy_color = galaxy_color, '
    print,  '       , /no_input_binary, rebin_factor = rebin_factor  '
    message, 'Check the Syntax Again !!!'
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; constants
pix = 0.259
pix_area = ( 0.259 * 0.259 )
gain = 3.0 
step_size = 1.0 
expand_factor = 2.0
color_loc = '/media/Elements/cings/test/E_color/cprof/'
filter_array = [ 'B', 'V', 'R', 'I' ]
photoerr_array = [ 0.06, 0.04, 0.03, 0.04 ]
mtype_array = [ 'as0', 'bs0', 'cs0', 'af0', 'bf0', 'cf0' ]
symbols = [ 46, 14, 16, 15, 17, 18, 19, 20 ]
;; color-color array from stellar population library 

;; harris 
age1_b = [ 6.49513, 6.77489, 7.18694, 7.37398 ] 
age1_v = [ 5.80353, 6.02197, 6.35145, 6.49440 ]
age1_r = [ 5.24418, 5.43401, 5.73215, 5.85450 ] 

age2_b = [ 6.98069, 7.26779, 7.64012, 7.83423 ] 
age2_v = [ 6.25104, 6.47280, 6.78398, 6.92183 ] 
age2_r = [ 5.67124, 5.86125, 6.14875, 6.26292 ] 

age3_b = [ 7.20203, 7.50746, 7.86371, 8.06691 ] 
age3_v = [ 6.45241, 6.69583, 6.98340, 7.13071 ] 
age3_r = [ 5.86125, 6.07558, 6.33552, 6.45712 ]

age4_b = [ 7.40821, 7.70470, 8.08279, 8.29362 ]
age4_v = [ 6.63815, 6.86587, 7.17313, 7.32674 ] 
age4_r = [ 6.03563, 6.22936, 6.50737, 6.63561 ]

met1_b = [ 6.47627, 6.96038, 7.18093, 7.38620 ] 
met1_v = [ 5.78132, 6.22794, 6.42876, 6.61399 ] 
met1_r = [ 5.22942, 5.65714, 5.84594, 6.02002 ] 

;; Johnson
;age4_B = [ 7.38700, 7.67781, 8.04986, 8.25627 ]
;age4_V = [ 6.54817, 6.77035, 7.07115, 7.21825 ] 
;age4_R = [ 6.06824, 6.26604, 6.54732, 6.67607 ] 
;
;age3_B = [ 7.18198, 7.48255, 7.83289, 8.03166 ] 
;age3_V = [ 6.36441, 6.60304, 6.88459, 7.02552 ] 
;age3_R = [ 5.89297, 6.11119, 6.37409, 6.49676 ] 
;
;age2_B = [ 6.96163, 7.24339, 7.61106, 7.80077 ] 
;age2_V = [ 6.16517, 6.38160, 6.68753, 6.81923 ] 
;age2_R = [ 5.70253, 5.89621, 6.18627, 6.30150 ] 
;
;age1_B = [ 6.47760, 6.75352, 7.15938, 7.34280 ] 
;age1_V = [ 5.72134, 5.93485, 6.25735, 6.39553 ] 
;age1_R = [ 5.27483, 5.46742, 5.76840, 5.89200 ]

age1_bv = ( age1_b - age1_v ) 
age1_vr = ( age1_v - age1_r ) 
age1_br = ( age1_b - age1_r ) 
age2_bv = ( age2_b - age2_v ) 
age2_vr = ( age2_v - age2_r ) 
age2_br = ( age2_b - age2_r ) 
age3_bv = ( age3_b - age3_v ) 
age3_vr = ( age3_v - age3_r ) 
age3_br = ( age3_b - age3_r ) 
age4_bv = ( age4_b - age4_v ) 
age4_vr = ( age4_v - age4_r ) 
age4_br = ( age4_b - age4_r ) 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
location_setting, setting
;;; define positions of files
    galfit = setting.galfit
  cat_head = setting.header
  datafile = setting.header
   fileloc = setting.workplace
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; prepare the input  
;; galaxy name
galaxy = strcompress( galaxy, /remove_all ) 
;; model type 
type = strcompress( type, /remove_all )
if ( where( mtype_array eq type ) eq -1 ) then begin 
    message, 'Wrong model type ! [ as0, af0, bs0, bf0, cs0, cf0 ] '
endif 
;; number of component 
ncomp = long( ncomp )
ncomp_string = strcompress( string( ncomp ), /remove_all )
;; rebin_factor 
if keyword_set( rebin_factor ) then begin 
    rebin_factor = float( rebin_factor )
endif else begin 
    rebin_factor = 4
endelse
;; get header information for three filters 
;; B-band 
get_head_info, galaxy, 'B', header1, exist
if ( exist eq 0 ) then begin 
    message, 'Can not find the header information for B-band'
endif
;; V-band 
get_head_info, galaxy, 'V', header2, exist
if ( exist eq 0 ) then begin 
    message, 'Can not find the header information for V-band'
endif
;; R-band 
get_head_info, galaxy, 'R', header3, exist
if ( exist eq 0 ) then begin 
    message, 'Can not find the header information for R-band'
endif
;; read information for B-band
std_name = header1.std_name
   r80_b = float( header1.r80 )
zpt_gsc1 = float( header1.zpt_gsc )
zpt_lan1 = float( header1.zpt_lan )
 cen_x_b = float( header1.cen_x )
 cen_y_b = float( header1.cen_y )
   ell_b = float( header1.ell_e )
    pa_b = float( header1.ell_pa )
old_expt_b = float( header1.old_expt )
if ( zpt_lan1 gt 15.0 ) then begin 
    magzpt_b = zpt_lan1
endif else begin 
    magzpt_b = zpt_gsc1
endelse
;; read information for V-band
   r80_v = float( header2.r80 )
zpt_gsc2 = float( header2.zpt_gsc )
zpt_lan2 = float( header2.zpt_lan )
 cen_x_v = float( header2.cen_x )
 cen_y_v = float( header2.cen_y )
   ell_v = float( header2.ell_e )
    pa_v = float( header2.ell_pa )
old_expt_v = float( header2.old_expt )
if ( zpt_lan2 gt 15.0 ) then begin 
    magzpt_v = zpt_lan2
endif else begin 
    magzpt_v = zpt_gsc2
endelse
;; read information for R-band
   r80_r = float( header3.r80 )
zpt_gsc3 = float( header3.zpt_gsc )
zpt_lan3 = float( header3.zpt_lan )
 cen_x_r = float( header3.cen_x )
 cen_y_r = float( header3.cen_y )
   ell_r = float( header3.ell_e )
    pa_r = float( header3.ell_pa )
old_expt_r = float( header3.old_expt )
if ( zpt_lan3 gt 15.0 ) then begin 
    magzpt_r = zpt_lan3
endif else begin 
    magzpt_r = zpt_gsc3
endelse
;; get the galaxy information 
get_ell_summary, galaxy, galaxy_info, exist
if ( exist eq 0 ) then begin 
    message, 'Something wrong with the galaxy information'
    return 
endif else begin 
    extinc_b = float( galaxy_info.a_b )
    extinc_v = float( galaxy_info.a_v )
    extinc_r = float( galaxy_info.a_r )
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; find the necessary file 
;; first define the model string and the location 
if ( standard_filter eq 'V' ) then begin 
    model_v = galaxy + '_V_model_' + ncomp_string + '0_fs0'  
    model_b = galaxy + '_B_color_' + ncomp_string + '0_' + type
    model_r = galaxy + '_R_color_' + ncomp_string + '0_' + type
    loc_b = color_loc + 'B_' + type + '/' 
    loc_v = color_loc + 'Vband/'
    loc_r = color_loc + 'R_' + type + '/' 
    mod_loc_b = loc_b + model_b
    mod_loc_v = loc_v + model_v
    mod_loc_r = loc_r + model_r
endif else begin 
    model_v = galaxy + '_V_color_' + ncomp_string + '0_' + type 
    loc_v = color_loc + 'V_' + standard_filter + '_' + type + '/' 
    if ( standard_filter eq 'B' ) then begin 
        model_b = galaxy + '_B_model_' + ncomp_string + '0_fs0' 
        model_r = galaxy + '_R_color_' + ncomp_string + '0_' + type 
        loc_b = color_loc + 'Bband/'
        loc_r = color_loc + 'R_' + type + '/'
    endif else begin 
        model_r = galaxy + '_R_model_' + ncomp_string + '0_fs0' 
        model_b = galaxy + '_B_color_' + ncomp_string + '0_' + type 
        loc_r = color_loc + 'Rband/'
        loc_b = color_loc + 'B_' + type + '/'
    endelse
    mod_loc_b = loc_b + model_b
    mod_loc_v = loc_v + model_v
    mod_loc_r = loc_r + model_r
endelse
;; find all the necessary files 
;; read file 
b_read = mod_loc_b + '.read'
v_read = mod_loc_v + '.read'
r_read = mod_loc_r + '.read'
;; note 
b_note  = mod_loc_b + '.note'
v_note  = mod_loc_v + '.note'
r_note  = mod_loc_r + '.note'
;; original image without sky background 
b_ori = mod_loc_b + '_nosky_ori.fit'
v_ori = mod_loc_v + '_nosky_ori.fit'
r_ori = mod_loc_r + '_nosky_ori.fit' 
;; find the model image 
if NOT keyword_set( psf ) then begin 
    ;; model image without psf convolution 
    b_mod = mod_loc_b + '_nosky_mod.fit'
    v_mod = mod_loc_v + '_nosky_mod.fit'
    r_mod = mod_loc_r + '_nosky_mod.fit'
endif else begin 
    ;; model image with psf convolution 
    b_mod = mod_loc_b + '_nosky_mof.fit'
    v_mod = mod_loc_v + '_nosky_mof.fit'
    r_mod = mod_loc_r + '_nosky_mof.fit'
endelse
;; old sbp profile 
b_sbp_prof = fileloc + galaxy + '/B/' + galaxy + '_B_sbp.dat'
v_sbp_prof = fileloc + galaxy + '/V/' + galaxy + '_V_sbp.dat'
r_sbp_prof = fileloc + galaxy + '/R/' + galaxy + '_R_sbp.dat'
;; find the model prof 
b_mod_prof = mod_loc_b + '_nosky_mod.prof'
v_mod_prof = mod_loc_v + '_nosky_mod.prof'
r_mod_prof = mod_loc_r + '_nosky_mod.prof'
;; exam whether the files are available 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; note files; if exist, read in 
if NOT file_test( b_note ) then begin 
    message, 'Can not find the note file for B-band at: ' + b_note 
endif else begin 
    read_galaxy_note, b_note, model_summary=b_mod_sum, comp_summary=b_comp_sum
    temp = b_comp_sum.r
    temp_sum = b_comp_sum
    s_radius = bsort( temp ) 
    temp_sum = b_comp_sum[ s_radius ] 
    b_comp_sum = temp_sum
endelse
if NOT file_test( v_note ) then begin 
    message, 'Can not find the note file for V-band at: ' + v_note 
endif else begin 
    read_galaxy_note, v_note, model_summary=v_mod_sum, comp_summary=v_comp_sum
    temp = v_comp_sum.r
    temp_sum = v_comp_sum
    s_radius = bsort( temp ) 
    temp_sum = v_comp_sum[ s_radius ] 
    v_comp_sum = temp_sum
endelse
if NOT file_test( r_note ) then begin 
    message, 'Can not find the note file for R-band at: ' + r_note 
endif else begin 
    read_galaxy_note, r_note, model_summary=r_mod_sum, comp_summary=r_comp_sum
    temp = r_comp_sum.r
    temp_sum = r_comp_sum
    s_radius = bsort( temp ) 
    temp_sum = r_comp_sum[ s_radius ] 
    r_comp_sum = temp_sum
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; old sbp 
if NOT file_test( b_sbp_prof ) then begin 
    message, 'Can not find the old sbp profile for B-band at : ' + b_sbp_prof 
endif else begin 
    read_sbp, b_sbp_prof, b_sbp, b_sbp_line 
    b_mu_min_lim = min( b_sbp.smag ) - b_sbp[ where( b_sbp.smag eq $
        min( b_sbp.smag ) ) ].err 
    b_mu_max_lim = max( b_sbp.smag ) + ( b_sbp[ where( b_sbp.smag eq $
        max( b_sbp.smag ) ) ].err / 2.0 )
endelse
if NOT file_test( v_sbp_prof ) then begin 
    message, 'Can not find the old sbp profile for V-band at : ' + v_sbp_prof 
endif else begin 
    read_sbp, v_sbp_prof, v_sbp, v_sbp_line 
    sma_min_lim = 0.0 
    sma_max_lim = max( v_sbp.sma )
    v_mu_min_lim = min( v_sbp.smag ) - v_sbp[ where( v_sbp.smag eq $
        min( v_sbp.smag ) ) ].err 
    v_mu_max_lim = max( v_sbp.smag ) + ( v_sbp[ where( v_sbp.smag eq $
        max( v_sbp.smag ) ) ].err / 2.0 )
endelse
if NOT file_test( r_sbp_prof ) then begin 
    message, 'Can not find the old sbp profile for R-band at : ' + r_sbp_prof 
endif else begin 
    read_sbp, r_sbp_prof, r_sbp, r_sbp_line 
    r_mu_min_lim = min( r_sbp.smag ) - r_sbp[ where( r_sbp.smag eq $
        min( r_sbp.smag ) ) ].err 
    r_mu_max_lim = max( r_sbp.smag ) + ( r_sbp[ where( r_sbp.smag eq $
        max( r_sbp.smag ) ) ].err / 2.0 )
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; image ori
if NOT file_test( b_ori ) then begin 
    message, 'Can not find the B-band original image without sky at: ' + b_ori 
endif else begin 
    b_img_ori = mrdfits( b_ori, 0, b_head_ori )
    b_size = size( b_img_ori, /dimension ) 
endelse
if NOT file_test( v_ori ) then begin 
    message, 'Can not find the V-band original image without sky at: ' + v_ori 
endif else begin 
    v_img_ori = mrdfits( v_ori, 0, v_head_ori )
    v_size = size( v_img_ori, /dimension ) 
endelse
if NOT file_test( r_ori ) then begin 
    message, 'Can not find the R-band original image without sky at: ' + r_ori 
endif else begin 
    r_img_ori = mrdfits( r_ori, 0, r_head_ori )
    r_size = size( r_img_ori, /dimension ) 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; image mod
if NOT file_test( b_mod ) then begin 
    message, 'Can not find the B-band model image without PSF at: ' + b_mod 
endif else begin 
    b_img_mod = mrdfits( b_mod, 0, b_head_mod )
    b_size2 = size( b_img_mod, /dimension ) 
endelse
if NOT file_test( v_mod ) then begin 
    message, 'Can not find the V-band model image without PSF at: ' + v_mod 
endif else begin 
    v_img_mod = mrdfits( v_mod, 0, v_head_mod )
    v_size2 = size( v_img_mod, /dimension ) 
endelse
if NOT file_test( r_mod ) then begin 
    message, 'Can not find the R-band model image without PSF at: ' + r_mod 
endif else begin 
    r_img_mod = mrdfits( r_mod, 0, r_head_mod )
    r_size2 = size( r_img_mod, /dimension ) 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; exam the size of these images 
if ( ( b_size[0] ne v_size[0] ) or ( b_size[1] ne v_size[1] ) or $ 
     ( b_size[0] ne r_size[0] ) or ( b_size[1] ne r_size[1] ) or $
     ( b_size[0] ne b_size2[0] ) or ( b_size[1] ne b_size2[1] ) or $ 
     ( b_size[0] ne v_size2[0] ) or ( b_size[1] ne v_size2[1] ) or $
     ( b_size[0] ne r_size2[0] ) or ( b_size[1] ne r_size2[1] ) ) then begin 
     message, 'The original images and model images should share the same size !!' 
endif 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; find the original images for all filters, read in and get the sky information 
b_cor = fileloc + galaxy + '/B/' + galaxy + '_B_cor.fit'
v_cor = fileloc + galaxy + '/V/' + galaxy + '_V_cor.fit'
r_cor = fileloc + galaxy + '/R/' + galaxy + '_R_cor.fit'
;; read the corrected frame and find out the sky an its error 
;; b-band 
if NOT file_test( b_cor ) then begin 
    message, 'Can not find the original image for B-band'
endif else begin 
    b_img_cor = mrdfits( b_cor, 0, b_cor_head ) 
    b_sky_new = float( sxpar( b_cor_head, 'SKY_NV' ) )
    b_sky_ner = float( sxpar( b_cor_head, 'SKY_NE' ) ) 
    print, '#########################################################'
    print, ' GALAXY : ' + galaxy 
    print, ' SKY_NE ( B-band ) : ' + string( b_sky_ner )
    print, '#######################################################'
    if ( b_sky_new le 1.0 ) then begin 
        print, '#######################################################'
        print, '!!!! Be careful with the sky for B-band'
        print, '#######################################################'
    endif
endelse
;; v-band 
if NOT file_test( v_cor ) then begin 
    message, 'Can not find the original image for B-band'
endif else begin 
    v_img_cor = mrdfits( v_cor, 0, v_cor_head ) 
    v_sky_new = float( sxpar( v_cor_head, 'SKY_NV' ) )
    v_sky_ner = float( sxpar( v_cor_head, 'SKY_NE' ) ) 
    print, ' SKY_NE ( V-band ) : ' + string( v_sky_ner )
    print, '#######################################################'
    if ( v_sky_new le 1.0 ) then begin 
        print, '#######################################################'
        print, '!!!! Be careful with the sky for V-band' 
        print, '#######################################################'
    endif
endelse
;; r-band
if NOT file_test( r_cor ) then begin 
    message, 'Can not find the original image for B-band'
endif else begin 
    r_img_cor = mrdfits( r_cor, 0, r_cor_head ) 
    r_sky_new = float( sxpar( r_cor_head, 'SKY_NV' ) )
    r_sky_ner = float( sxpar( r_cor_head, 'SKY_NE' ) ) 
    print, ' SKY_NE ( R-band ) : ' + string( r_sky_ner )
    print, '#######################################################'
    if ( r_sky_new le 1.0 ) then begin 
        print, '#######################################################'
        print, '!!!! Be careful with the sky for R-band' 
        print, '#######################################################'
    endif
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; model profile 
if NOT file_test( b_mod_prof ) then begin 
    message, 'Can not find the model profile for B-band at : ' + b_mod_prof 
endif else begin 
    read_profile, b_mod_prof, b_mod_smag, b_mod_line 
    b_mod_sma = b_mod_smag.sma 
    b_mod_mu  = -2.5 * alog10( b_mod_smag.intens / ( pix_area * old_expt_b ) $ 
        ) + magzpt_b - extinc_b
    b_mod_sig = sqrt( ( b_mod_smag.intens_err )^2.0 + ( b_sky_ner )^2.0 )
    b_mod_err = 2.5 * alog10( 1.0 + ( b_mod_sig / b_mod_smag.intens ) )
endelse
if NOT file_test( v_mod_prof ) then begin 
    message, 'Can not find the model profile for B-band at : ' + v_mod_prof 
endif else begin 
    read_profile, v_mod_prof, v_mod_smag, v_mod_line 
    v_mod_sma = v_mod_smag.sma 
    v_mod_mu  = -2.5 * alog10( v_mod_smag.intens / ( pix_area * old_expt_b ) $
        ) + magzpt_b - extinc_b
    v_mod_sig = sqrt( ( v_mod_smag.intens_err )^2.0 + ( v_sky_ner )^2.0 )
    v_mod_err = 2.5 * alog10( 1.0 + ( v_mod_sig / v_mod_smag.intens ) )
endelse
if NOT file_test( r_mod_prof ) then begin 
    message, 'Can not find the model profile for B-band at : ' + r_mod_prof 
endif else begin 
    read_profile, r_mod_prof, r_mod_smag, r_mod_line 
    r_mod_sma = r_mod_smag.sma 
    r_mod_mu  = -2.5 * alog10( r_mod_smag.intens / ( pix_area * old_expt_b ) $
        ) + magzpt_b - extinc_b
    r_mod_sig = sqrt( ( r_mod_smag.intens_err )^2.0 + ( r_sky_ner )^2.0 )
    r_mod_err = 2.5 * alog10( 1.0 + ( r_mod_sig / r_mod_smag.intens ) )
endelse
;; mod color profile 
bv_mod_cprof = ( b_mod_mu - v_mod_mu ) 
bv_mod_eprof = sqrt( b_mod_err^2.0 + v_mod_err^2.0 )
vr_mod_cprof = ( v_mod_mu - r_mod_mu ) 
vr_mod_eprof = sqrt( v_mod_err^2.0 + r_mod_err^2.0 )
br_mod_cprof = ( b_mod_mu - r_mod_mu ) 
br_mod_eprof = sqrt( b_mod_err^2.0 + r_mod_err^2.0 )
;; min and max within a reasonable radius range 
bv_prof_min = min( bv_mod_cprof[ where( v_mod_sma le sma_max_lim ) ] ) 
bv_prof_max = max( bv_mod_cprof[ where( ( v_mod_sma * 0.259 )^0.25 le 3.0 ) ] ) 
vr_prof_min = min( vr_mod_cprof[ where( v_mod_sma le sma_max_lim ) ] ) 
vr_prof_max = max( vr_mod_cprof[ where( ( v_mod_sma * 0.259 )^0.25 le 3.0 ) ] ) 
br_prof_min = min( br_mod_cprof[ where( v_mod_sma le sma_max_lim ) ] ) 
br_prof_max = max( br_mod_cprof[ where( ( v_mod_sma * 0.259 )^0.25 le 3.0 ) ] ) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; find the mask image for three bans 
b_msk = fileloc + galaxy + '/B/' + galaxy + '_B_mm.fits'
v_msk = fileloc + galaxy + '/V/' + galaxy + '_V_mm.fits'
r_msk = fileloc + galaxy + '/R/' + galaxy + '_R_mm.fits'
;; test if the files are there and readin 
;; b-band 
if NOT file_test( b_cor ) then begin 
    message, 'Can not the B-band mask image at : ' + b_cor 
endif else begin 
    b_img_msk = mrdfits( b_msk, 0 )
endelse
;; v-band 
if NOT file_test( v_cor ) then begin 
    message, 'Can not the V-band mask image at : ' + v_cor 
endif else begin 
    v_img_msk = mrdfits( v_msk, 0 )
endelse
;; r-band 
if NOT file_test( r_cor ) then begin 
    message, 'Can not the R-band mask image at : ' + r_cor 
endif else begin 
    r_img_msk = mrdfits( r_msk, 0 )
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; exam the number of components from input with the one found in note file 
;; if consistent, get information for each component 
if ( ncomp ne n_elements( b_comp_sum.r ) or $
     ncomp ne n_elements( v_comp_sum.r ) or $ 
     ncomp ne n_elements( r_comp_sum.r ) ) then begin 
    message, 'The number of component is inconsistent between the input and ' +$
        'the note file, Check again !! '
endif else begin 
    ;; re of components in arcsec
    b_comp_rad = b_comp_sum.r 
    v_comp_rad = v_comp_sum.r 
    r_comp_rad = r_comp_sum.r
    ;; mue of component 
    b_comp_mue = ( b_comp_sum.mue - extinc_b )
    v_comp_mue = ( v_comp_sum.mue - extinc_v )
    r_comp_mue = ( r_comp_sum.mue - extinc_r )
    ;; luminosity fraction 
    b_comp_frac = b_comp_sum.frac 
    v_comp_frac = v_comp_sum.frac 
    r_comp_frac = r_comp_sum.frac 
    ;; bv color 
    bv_comp_color = ( ( b_comp_sum.mag - extinc_b ) - ( v_comp_sum.mag - $ 
        extinc_v ) )
    ;; vr color 
    vr_comp_color = ( ( v_comp_sum.mag - extinc_v ) - ( r_comp_sum.mag - $ 
        extinc_r ) )
    ;; br color 
    br_comp_color = ( ( b_comp_sum.mag - extinc_b ) - ( r_comp_sum.mag - $ 
        extinc_r ) )
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the size array for the symbols of different components 
;; Using V-band as standard_filter
smallest = 1.5
largest  = 3.0
sep1 = ( largest - smallest ) 
sep2 = ( alog10( max( v_comp_frac * 100 ) ) - $ 
    alog10( min( v_comp_frac )* 100 ) )
comp_symsize = ( sep1 / sep2 ) * ( alog10( v_comp_frac * 100 ) - $
    alog10( min( v_comp_frac ) * 100 ) ) + smallest
;; calculate the total model magnitude and average color for three bands 
b_mod_flux = total( b_img_mod ) 
v_mod_flux = total( v_img_mod ) 
r_mod_flux = total( r_img_mod ) 
b_mod_tmag = -2.5 * alog10( b_mod_flux / old_expt_b ) + magzpt_b - extinc_b
v_mod_tmag = -2.5 * alog10( v_mod_flux / old_expt_v ) + magzpt_v - extinc_v
r_mod_tmag = -2.5 * alog10( r_mod_flux / old_expt_r ) + magzpt_r - extinc_r
bv_mod_color_avg = ( b_mod_tmag - v_mod_tmag )
vr_mod_color_avg = ( v_mod_tmag - r_mod_tmag )
br_mod_color_avg = ( b_mod_tmag - r_mod_tmag )
print, '############################################################'
print, 'Average Model B-V : ', bv_mod_color_avg
print, 'Average Model V-R : ', vr_mod_color_avg
print, 'Average Model B-R : ', br_mod_color_avg
print, '############################################################'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; make the color map 
position_1 = [ 0.100, 0.15, 0.495, 0.85 ]
position_2 = [ 0.505, 0.15, 0.900, 0.85 ]

rb = rebin_factor
new_area = ( pix_area * rb^2.0 ) 
newx = long( b_size[0] / rb )
newy = long( b_size[1] / rb )
msk_thre = ( rb * rb ) 

;;; add Poisson noise to the model image 
b_img_mod_ns = poidev( ( b_img_mod * gain ), SEED=1.0 ) + b_img_mod
v_img_mod_ns = poidev( ( v_img_mod * gain ), SEED=1.0 ) + v_img_mod
r_img_mod_ns = poidev( ( r_img_mod * gain ), SEED=1.0 ) + r_img_mod
;; rebin the data and make the surfacer brightness map 
if ( rb gt 1 ) then begin 
    ;; b-band 
    b_img_ori_rb = frebin( b_img_ori, newx, newy, /total ) 
    b_img_mod_rb = frebin( b_img_mod, newx, newy, /total ) 
    b_img_mod_rb_ns = frebin( b_img_mod_ns, newx, newy, /total )
    b_img_msk_rb = frebin( b_img_msk, newx, newy, /total ) 
    b_img_msk_rb[ where( b_img_msk_rb lt msk_thre ) ] = 0 
    b_img_msk_rb[ where( b_img_msk_rb ge msk_thre ) ] = 1
    ;; v-band 
    v_img_ori_rb = frebin( v_img_ori, newx, newy, /total ) 
    v_img_mod_rb = frebin( v_img_mod, newx, newy, /total ) 
    v_img_mod_rb_ns = frebin( v_img_mod_ns, newx, newy, /total )
    v_img_msk_rb = frebin( v_img_msk, newx, newy, /total ) 
    v_img_msk_rb[ where( v_img_msk_rb lt msk_thre ) ] = 0 
    v_img_msk_rb[ where( v_img_msk_rb ge msk_thre ) ] = 1
    ;; r-band 
    r_img_ori_rb = frebin( r_img_ori, newx, newy, /total ) 
    r_img_mod_rb = frebin( r_img_mod, newx, newy, /total ) 
    r_img_mod_rb_ns = frebin( r_img_mod_ns, newx, newy, /total )
    r_img_msk_rb = frebin( r_img_msk, newx, newy, /total ) 
    r_img_msk_rb[ where( r_img_msk_rb lt msk_thre ) ] = 0 
    r_img_msk_rb[ where( r_img_msk_rb ge msk_thre ) ] = 1
endif else begin 
    ;; b-band 
    b_img_ori_rb = b_img_ori 
    b_img_mod_rb = b_img_mod
    b_img_mod_rb_ns = b_img_mod_ns
    ;; v-band 
    v_img_ori_rb = v_img_ori
    v_img_mod_rb = v_img_mod
    v_img_mod_rb_ns = v_img_mod_ns
    ;; r-band 
    r_img_ori_rb = r_img_ori
    r_img_mod_rb = r_img_mod
    r_img_mod_rb_ns = r_img_mod_ns
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; make the surface brightness map 
;; b-band 
b_smag_ori = -2.5 * alog10( b_img_ori_rb / ( old_expt_b * new_area ) ) + $
    magzpt_b - extinc_b 
b_smag_mod = -2.5 * alog10( b_img_mod_rb / ( old_expt_b * new_area ) ) + $
    magzpt_b - extinc_b 
b_smag_mod_ns = -2.5 * alog10( b_img_mod_rb_ns / $
        ( old_expt_b * new_area ) ) + magzpt_b - extinc_b 
;; v-band 
v_smag_ori = -2.5 * alog10( v_img_ori_rb / ( old_expt_v * new_area ) ) + $
    magzpt_v - extinc_v 
v_smag_mod = -2.5 * alog10( v_img_mod_rb / ( old_expt_v * new_area ) ) + $
    magzpt_v - extinc_v 
v_smag_mod_ns = -2.5 * alog10( v_img_mod_rb_ns / $
    ( old_expt_v * new_area ) ) + magzpt_v - extinc_v 
;; b-band 
r_smag_ori = -2.5 * alog10( r_img_ori_rb / ( old_expt_r * new_area ) ) + $
    magzpt_r - extinc_r 
r_smag_mod = -2.5 * alog10( r_img_mod_rb / ( old_expt_r * new_area ) ) + $
    magzpt_r - extinc_r 
r_smag_mod_ns = -2.5 * alog10( r_img_mod_rb_ns / $
    ( old_expt_r * new_area ) ) + magzpt_r - extinc_r 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; find the pixels without finite value (NaN) 
b_nan_mask = finite( b_smag_ori, /NaN )
v_nan_mask = finite( v_smag_ori, /NaN )
r_nan_mask = finite( r_smag_ori, /NaN )
;; update the mask image for NaN pixels 
b_img_msk_rb[ where( b_nan_mask eq 1 ) ] = 1 
v_img_msk_rb[ where( v_nan_mask eq 1 ) ] = 1 
r_img_msk_rb[ where( r_nan_mask eq 1 ) ] = 1 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; build the color map for both original image and model 
;; b-v color
bv_ori_color = ( b_smag_ori - v_smag_ori )
bv_mod_color = ( b_smag_mod - v_smag_mod )
bv_mod_color_ns = ( b_smag_mod_ns - v_smag_mod_ns )
bv_ori_color_finite = bv_ori_color 
bv_ori_color_finite[ where( b_nan_mask eq 1 ) ] = -99
bv_ori_color_finite[ where( v_nan_mask eq 1 ) ] = -99
bv_ori_color_finite[ where( r_nan_mask eq 1 ) ] = -99
bv_ori_color_finite[ where( v_img_msk_rb eq 1 ) ] = -99
;; v-r color
vr_ori_color = ( v_smag_ori - r_smag_ori )
vr_mod_color = ( v_smag_mod - r_smag_mod )
vr_mod_color_ns = ( v_smag_mod_ns - r_smag_mod_ns )
vr_ori_color_finite = vr_ori_color 
vr_ori_color_finite[ where( b_nan_mask eq 1 ) ] = -99
vr_ori_color_finite[ where( v_nan_mask eq 1 ) ] = -99
vr_ori_color_finite[ where( r_nan_mask eq 1 ) ] = -99
vr_ori_color_finite[ where( v_img_msk_rb eq 1 ) ] = -99
;; b-r color
br_ori_color = ( b_smag_ori - r_smag_ori )
br_mod_color = ( b_smag_mod - r_smag_mod )
br_mod_color_ns = ( b_smag_mod_ns - r_smag_mod_ns )
br_ori_color_finite = br_ori_color 
br_ori_color_finite[ where( b_nan_mask eq 1 ) ] = -99
br_ori_color_finite[ where( v_nan_mask eq 1 ) ] = -99
br_ori_color_finite[ where( r_nan_mask eq 1 ) ] = -99
br_ori_color_finite[ where( v_img_msk_rb eq 1 ) ] = -99
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define a region for pixels that belong to the galaxy 
ab = ( 1.00 / ( 1.00 - ell_v ) )
ncx = ( cen_x_v / rb )
ncy = ( cen_y_v / rb )
rlim = ( r80_v * expand_factor / rb ) < min( [newx, newy] ) 
r80n = ( r80_v / rb ) > max( v_comp_rad / rb ) 
;; select a region to show the color map 
x0 = ( ncx - rlim ) > 1 
x1 = ( ncx + rlim ) < ( newx - 1 ) 
y0 = ( ncy - rlim ) > 1 
y1 = ( ncy + rlim ) < ( newy - 1 )
x0 = long( x0 )
x1 = long( x1 ) 
y0 = long( y0 ) 
y1 = long( y1 )
print, '###############################################################'
print, 'The Color Map Regon : ', x0, x1, y0, y1
print, '###############################################################'
bv_ori_color_show = bv_ori_color_finite[x0:x1,y0:y1]
vr_ori_color_show = vr_ori_color_finite[x0:x1,y0:y1]
br_ori_color_show = br_ori_color_finite[x0:x1,y0:y1]
bv_mod_color_show = bv_mod_color_ns[x0:x1,y0:y1]
vr_mod_color_show = vr_mod_color_ns[x0:x1,y0:y1]
br_mod_color_show = br_mod_color_ns[x0:x1,y0:y1]
;; make the distance mask  
dist_ellipse, ellip_mask, [newx,newy], ncx, ncy, ab, pa_v, /double
ell_msk = ellip_mask
;; select useful pixels 
index_galaxy = where( ( ellip_mask le rlim ) and ( b_img_msk_rb eq 0 ) and $ 
    ( v_img_msk_rb eq 0 ) and ( r_img_msk_rb eq 0 ) )
index_galaxy_r80 = where( ( ellip_mask le r80n ) and ( b_img_msk_rb eq 0 ) and $ 
    ( v_img_msk_rb eq 0 ) and ( r_img_msk_rb eq 0 ) )
index_inside = where( ellip_mask le rlim ) 
index_maskout = where( ( ellip_mask le rlim ) and ( b_img_msk_rb ne 0 ) and $ 
    ( v_img_msk_rb ne 0 ) and ( r_img_msk_rb ne 0 ) )
mask_frac = ( n_elements( index_maskout ) / n_elements( index_inside ) )
print, '################################################################'
print, 'There are : ', n_elements( index_galaxy ), ' pixels for galaxy '
print, 'The fraction of masked out pixels are : ', ( mask_frac * 100.00 )
print, '################################################################'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; construct arrays for useful information 
;; dist 
ell_dist = ell_msk[ index_galaxy ]
;; b-band 
b_mu_ori = b_smag_ori[ index_galaxy ]
b_mu_mod = b_smag_mod[ index_galaxy ]
;; v-band 
v_mu_ori = v_smag_ori[ index_galaxy ]
v_mu_mod = v_smag_mod[ index_galaxy ]
;; r-band 
r_mu_ori = r_smag_ori[ index_galaxy ]
r_mu_mod = r_smag_mod[ index_galaxy ]
;; color 
if NOT keyword_set( add_noise ) then begin 
    bv_ori_color_use = bv_ori_color[ index_galaxy ]
    vr_ori_color_use = vr_ori_color[ index_galaxy ]
    br_ori_color_use = br_ori_color[ index_galaxy ]
    bv_mod_color_use = bv_mod_color[ index_galaxy ]
    vr_mod_color_use = vr_mod_color[ index_galaxy ]
    br_mod_color_use = br_mod_color[ index_galaxy ]
    ;; for color-color map
    bv_ori_color_r80 = bv_ori_color[ index_galaxy_r80 ]
    bv_mod_color_r80 = bv_mod_color[ index_galaxy_r80 ]
    vr_ori_color_r80 = vr_ori_color[ index_galaxy_r80 ]
    vr_mod_color_r80 = vr_mod_color[ index_galaxy_r80 ]
    br_ori_color_r80 = br_ori_color[ index_galaxy_r80 ]
    br_mod_color_r80 = br_mod_color[ index_galaxy_r80 ]
endif else begin 
    bv_ori_color_use = bv_ori_color[ index_galaxy ]
    vr_ori_color_use = vr_ori_color[ index_galaxy ]
    br_ori_color_use = br_ori_color[ index_galaxy ]
    bv_mod_color_use = bv_mod_color_ns[ index_galaxy ]
    vr_mod_color_use = vr_mod_color_ns[ index_galaxy ]
    br_mod_color_use = br_mod_color_ns[ index_galaxy ]
    ;; for color-color map
    bv_ori_color_r80 = bv_ori_color[ index_galaxy_r80 ]
    bv_mod_color_r80 = bv_mod_color_ns[ index_galaxy_r80 ]
    vr_ori_color_r80 = vr_ori_color[ index_galaxy_r80 ]
    vr_mod_color_r80 = vr_mod_color_ns[ index_galaxy_r80 ]
    br_ori_color_r80 = br_ori_color[ index_galaxy_r80 ]
    br_mod_color_r80 = br_mod_color_ns[ index_galaxy_r80 ]
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the useful range for three colors 
;; bv 
min_bv = ( bv_prof_min - 0.05 ) < ( min( bv_comp_color ) - 0.30 ) 
max_bv = ( bv_prof_max + 0.05 ) > ( max( bv_comp_color ) + 0.25 )
min_bv_r80 = ( min( bv_mod_color_r80 ) - 0.1 ) < ( min( bv_comp_color ) - 0.1 )
min_bv_r80 = min_bv_r80 > ( min( bv_comp_color ) - 0.4 )
max_bv_r80 = ( max( bv_mod_color_r80 ) + 0.1 ) > ( max( bv_comp_color ) + 0.1 )
max_bv_r80 = max_bv_r80 < ( max( bv_comp_color ) + 0.4 )
off1 = ( bv_mod_color_avg - min_bv ) 
off2 = ( max_bv - bv_mod_color_avg )
offr = ( off2 / off1 ) 
offd = ( off2 - off1 ) 
if ( ( off1 gt 0.0 ) and ( off2 gt 0.0 ) and ( offr ge 1.2 ) ) then begin 
    if ( offd ge off1 ) then begin 
        min_bv = ( min_bv - ( offd / 2.0 ) ) 
        max_bv = ( max_bv - ( offd / 2.5 ) ) 
    endif else begin 
        min_bv = ( min_bv - offd ) 
        max_bv = ( max_bv - 0.02 ) 
    endelse
endif
print, '###############################################################'
print, 'The B-V Color Range: ' + string( min_bv ) + ' --> ' + string( max_bv )
;; vr 
min_vr = ( vr_prof_min - 0.05 ) < ( min( vr_comp_color ) - 0.30 ) 
max_vr = ( vr_prof_max + 0.05 ) > ( max( vr_comp_color ) + 0.25 )
min_vr_r80 = ( min( vr_mod_color_r80 ) - 0.1 ) < ( min( vr_comp_color ) - 0.1 )
min_vr_r80 = min_vr_r80 > ( min( vr_comp_color ) - 0.4 )
max_vr_r80 = ( max( vr_mod_color_r80 ) + 0.1 ) > ( max( vr_comp_color ) + 0.1 )
max_vr_r80 = max_vr_r80 < ( max( vr_comp_color ) + 0.4 )
off1 = ( vr_mod_color_avg - min_vr ) 
off2 = ( max_vr - vr_mod_color_avg )
offr = ( off2 / off1 ) 
offd = ( off2 - off1 ) 
if ( ( off1 gt 0.0 ) and ( off2 gt 0.0 ) and ( offr ge 1.2 ) ) then begin 
    if ( offd ge off1 ) then begin 
        min_vr = ( min_vr - ( offd / 2.0 ) ) 
        max_vr = ( max_vr - ( offd / 2.5 ) ) 
    endif else begin 
        min_vr = ( min_vr - offd ) 
        max_vr = ( max_vr - 0.02 ) 
    endelse
endif
print, 'The V-R Color Range: ' + string( min_vr ) + ' --> ' + string( max_vr )
;; br 
min_br = ( br_prof_min - 0.1 ) < ( min( br_comp_color ) - 0.45 ) 
max_br = ( br_prof_max + 0.1 ) > ( max( br_comp_color ) + 0.30 )
min_br_r80 = ( min( br_mod_color_r80 ) - 0.1 ) < ( min( br_comp_color ) - 0.2 )
min_br_r80 = min_br_r80 > ( min( br_comp_color ) - 0.4 )
max_br_r80 = ( max( br_mod_color_r80 ) + 0.1 ) > ( max( br_comp_color ) + 0.2 )
max_br_r80 = max_br_r80 < ( max( br_comp_color ) + 0.4 )
off1 = ( br_mod_color_avg - min_br ) 
off2 = ( max_br - br_mod_color_avg )
offr = ( off2 / off1 ) 
offd = ( off2 - off1 ) 
if ( ( off1 gt 0.0 ) and ( off2 gt 0.0 ) and ( offr ge 1.2 ) ) then begin 
    if ( offd ge off1 ) then begin 
        min_br = ( min_br - ( offd / 2.0 ) ) 
        max_br = ( max_br - ( offd / 2.5 ) ) 
    endif else begin 
        min_br = ( min_br - offd ) 
        max_br = ( max_br - 0.02 ) 
    endelse
endif
print, 'The B-R Color Range: ' + string( min_br ) + ' --> ' + string( max_br )
print, '###############################################################'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the useful range for surface brightness 
min_mu = ( v_mu_min_lim - 0.5 ) > ( min( v_comp_mue ) - 3.0 ) 
max_mu = ( v_mu_max_lim + 0.5 ) < ( max( v_comp_mue ) + 3.5 ) 
max_mu = max_mu > max( [ 24.1, ( max( v_comp_mue ) + 0.5 ) ] )
min_mu_bin = ( min_mu + 0.5 ) < ( min( v_comp_mue ) - 0.5 ) 
max_mu_bin = ( max_mu - 1.0 ) > ( max( v_comp_mue ) + 1.0 ) 
max_mu_bin = max_mu_bin < ( max( v_comp_mue ) + 3.0 ) 
min_mu_b = ( min( b_mu_mod ) - 0.5 ) > ( min( b_comp_mue ) - 0.6 ) 
min_mu_b = ( min_mu_b + 0.5 ) < ( min( b_comp_mue ) - 0.5 )
min_mu_r = ( min( r_mu_mod ) - 0.5 ) > ( min( r_comp_mue ) - 0.6 ) 
min_mu_r = ( min_mu_r + 0.5 ) < ( min( r_comp_mue ) - 0.5 )
;; define bins along the surface brightness axis 
step0 = min_mu_bin 
step1 = 1.0 
step2 = max_mu_bin 
nbins = long( ( step2 - step0 ) / step1 ) 
print, '###############################################################' 
print, 'Min Mu V-band = ', min_mu
print, 'Max Mu V-band = ', max_mu 
print, 'Min Mu Bin = ', min_mu_bin
print, 'Max Mu Bin = ', max_mu_bin 
print, 'Bin Numbers = ', nbins
print, '###############################################################' 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; calculate the running mean of different colors along the V-ban mu 
bin_v_mu = fltarr( nbins ) 
bin_bv_ori = fltarr( nbins ) 
bin_bv_ori_s = fltarr( nbins )
bin_bv_mod = fltarr( nbins ) 
bin_bv_mod_s = fltarr( nbins )
bin_vr_ori = fltarr( nbins ) 
bin_vr_ori_s = fltarr( nbins )
bin_vr_mod = fltarr( nbins ) 
bin_vr_mod_s = fltarr( nbins )
bin_br_ori = fltarr( nbins ) 
bin_br_ori_s = fltarr( nbins )
bin_br_mod = fltarr( nbins ) 
bin_br_mod_s = fltarr( nbins )
print, '###############################################################'
for i = 0, ( nbins - 2 ), 1 do begin  
    bin_v_mu[i] = ( step0 + step1 / 2.0 ) 
    edge1 = step0 
    edge2 = ( step0 + step1 ) 
    print, ' BIN Number : ' + string( i + 1 ) 
    print, ' BIN Edge1  : ' + string( edge1 )
    print, ' BIN Edge2  : ' + string( edge2 )
    ;; b-v color 
    index_ori_use = where( ( v_mu_ori ge edge1 ) and $
        ( v_mu_ori lt edge2 ) )
    index_mod_use = where( ( v_mu_mod ge edge1 ) and $ 
        ( v_mu_mod lt edge2 ) )
    if ( ( index_ori_use[0] eq -1 ) or ( index_mod_use[0] eq -1 ) ) then begin 
        bin_bv_ori[i] = -99 
        bin_bv_ori_s[i] = -99 
        bin_bv_mod[i] = -99 
        bin_bv_mod_s[i] = -99 
    endif else begin 
        ori_use = bv_ori_color_use[ index_ori_use ]
        mod_use = bv_mod_color_use[ index_mod_use ]
        nn = n_elements( ori_use ) 
        ;; for ori 
        resistant_mean, ori_use, 2, aa, bb, cc 
        dd = bb * sqrt( ( nn - cc ) - 1 )
        bin_bv_ori[i] = aa 
        bin_bv_ori_s[i] = ( dd / 2.0 ) 
        ;; for model 
        resistant_mean, mod_use, 2, aa, bb, cc 
        dd = bb * sqrt( ( nn - cc ) - 1 )
        bin_bv_mod[i] = aa 
        bin_bv_mod_s[i] = ( dd / 2.0 ) 
    endelse 
    ;; v-r color 
    index_ori_use = where( ( v_mu_ori ge edge1 ) and $
        ( v_mu_ori lt edge2 ) )
    index_mod_use = where( ( v_mu_mod ge edge1 ) and $ 
        ( v_mu_mod lt edge2 ) )
    if ( ( index_ori_use[0] eq -1 ) or ( index_mod_use[0] eq -1 ) ) then begin 
        bin_vr_ori[i] = -99 
        bin_vr_ori_s[i] = -99 
        bin_vr_mod[i] = -99 
        bin_vr_mod_s[i] = -99 
    endif else begin 
        ori_use = vr_ori_color_use[ index_ori_use ]
        mod_use = vr_mod_color_use[ index_mod_use ]
        nn = n_elements( ori_use ) 
        ;; for ori 
        resistant_mean, ori_use, 2, aa, bb, cc 
        dd = bb * sqrt( ( nn - cc ) - 1 )
        bin_vr_ori[i] = aa 
        bin_vr_ori_s[i] = ( dd / 2.0 ) 
        ;; for model 
        resistant_mean, mod_use, 2, aa, bb, cc 
        dd = bb * sqrt( ( nn - cc ) - 1 )
        bin_vr_mod[i] = aa 
        bin_vr_mod_s[i] = ( dd / 2.0 ) 
    endelse 
    ;; b-r color 
    index_ori_use = where( ( v_mu_ori ge edge1 ) and $
        ( v_mu_ori lt edge2 ) )
    index_mod_use = where( ( v_mu_mod ge edge1 ) and $ 
        ( v_mu_mod lt edge2 ) )
    if ( ( index_ori_use[0] eq -1 ) or ( index_mod_use[0] eq -1 ) ) then begin 
        bin_br_ori[i] = -99 
        bin_br_ori_s[i] = -99 
        bin_br_mod[i] = -99 
        bin_br_mod_s[i] = -99 
    endif else begin 
        ori_use = br_ori_color_use[ index_ori_use ]
        mod_use = br_mod_color_use[ index_mod_use ]
        nn = n_elements( ori_use ) 
        ;; for ori 
        resistant_mean, ori_use, 2, aa, bb, cc 
        dd = bb * sqrt( ( nn - cc ) - 1 )
        bin_br_ori[i] = aa 
        bin_br_ori_s[i] = ( dd / 2.0 ) 
        ;; for model 
        resistant_mean, mod_use, 2, aa, bb, cc 
        dd = bb * sqrt( ( nn - cc ) - 1 )
        bin_br_mod[i] = aa 
        bin_br_mod_s[i] = ( dd / 2.0 ) 
    endelse 
    step0 = edge2 
endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; only keep the useful bin 
;; b-v color 
bin_use = where( ( bin_bv_ori ne -99 ) and ( bin_vr_ori ne -99 ) and $ 
    ( bin_br_ori ne -99 ) ) 
if ( bin_use[0] ne -1 ) then begin 
    bin_v_mu = bin_v_mu[ bin_use ]
    bin_bv_ori = bin_bv_ori[ bin_use ]
    bin_bv_ori_s = bin_bv_ori_s[ bin_use ]
    bin_bv_mod = bin_bv_mod[ bin_use ]
    bin_bv_mod_s = bin_bv_mod_s[ bin_use ]
    bin_vr_ori = bin_vr_ori[ bin_use ]
    bin_vr_ori_s = bin_vr_ori_s[ bin_use ]
    bin_vr_mod = bin_vr_mod[ bin_use ]
    bin_vr_mod_s = bin_vr_mod_s[ bin_use ]
    bin_br_ori = bin_br_ori[ bin_use ]
    bin_br_ori_s = bin_br_ori_s[ bin_use ]
    bin_br_mod = bin_br_mod[ bin_use ]
    bin_br_mod_s = bin_br_mod_s[ bin_use ]
endif else begin 
    print, '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    print, 'Somthing wrong with the Mu-Color Bins !!!! ' 
    print, '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
endelse
print, '###################################################################'
print, 'BIN_V_MU', bin_v_mu
print, 'BIN_BV_ORI', bin_bv_ori
print, 'BIN_BV_MOD', bin_bv_mod
print, 'BIN_VR_ORI', bin_vr_ori
print, 'BIN_VR_MOD', bin_vr_mod
print, 'BIN_BR_ORI', bin_br_ori
print, 'BIN_BR_MOD', bin_br_mod
print, '###################################################################'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; calculate the photometric error 
bin_v_mu_err = fltarr( nbins )
bin_bv_err = fltarr( nbins )
bin_vr_err = fltarr( nbins )
bin_br_err = fltarr( nbins )
for i = 0, ( nbins - 2 ), 1 do begin 
    print, '##############################################################'
    print, 'Bin Number : ', string( i + 1 )
    ;; mu 
    mu_b = bin_v_mu[i]  
    mu_v = ( bin_v_mu[i] + bin_bv_ori[i] )
    mu_r = ( bin_v_mu[i] - bin_vr_ori[i] )  
    print, ' Mu for B,V,R-band : ', string( mu_b ), string( mu_v ), $
        string( mu_r ) 
    ;; flux 
    f_b = ( 10.0^( ( magzpt_b - mu_b ) / 2.50 ) * old_expt_b * pix_area ) 
    f_v = ( 10.0^( ( magzpt_v - mu_v ) / 2.50 ) * old_expt_v * pix_area ) 
    f_r = ( 10.0^( ( magzpt_r - mu_r ) / 2.50 ) * old_expt_r * pix_area ) 
    ;; flux err 
    err_b = abs( 2.50 * alog10( 1.0 + ( b_sky_ner / f_b ) ) )
    err_v = abs( 2.50 * alog10( 1.0 + ( v_sky_ner / f_v ) ) )
    err_r = abs( 2.50 * alog10( 1.0 + ( r_sky_ner / f_r ) ) )
    print, ' Sky Error for B,V,R-band : ', err_b, err_v, err_r
    err_b = sqrt( err_b^2 + ( photoerr_array[0] / 2.0 )^2.0 ) 
    err_v = sqrt( err_v^2 + ( photoerr_array[1] / 2.0 )^2.0 ) 
    err_r = sqrt( err_r^2 + ( photoerr_array[2] / 2.0 )^2.0 ) 
    print, ' Combined Error for B,V,R-band : ', err_b, err_v, err_r
    bin_v_mu_err[i] = err_v 
    ;; color err 
    bin_bv_err[i] = sqrt( err_b^2.0 + err_v^2.0 ) 
    bin_vr_err[i] = sqrt( err_r^2.0 + err_v^2.0 ) 
    bin_br_err[i] = sqrt( err_b^2.0 + err_r^2.0 ) 
    print, 'B-V Error : ', bin_bv_err[i]
    print, 'V-R Error : ', bin_vr_err[i]
    print, 'B-R Error : ', bin_br_err[i]
endfor 
print, '###################################################################'
print, 'BIN_V_MU_ERR : ', bin_v_mu_err
print, 'BIN_BV_ERR : ', bin_bv_err 
print, 'BIN_VR_ERR : ', bin_vr_err 
print, 'BIN_BR_ERR : ', bin_br_err 
print, '###################################################################'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; make room for large error 
;; b-v 
max_bv_err = max( bin_bv_err )
if ( max_bv_err > 0.25 ) then begin 
    max_bv = ( max_bv + ( max_bv_err / 2.0 ) + 0.02 ) 
endif
;; v-r 
max_vr_err = max( bin_vr_err )
if ( max_vr_err > 0.25 ) then begin 
    max_vr = ( max_vr + ( max_vr_err / 2.0 ) + 0.02 ) 
endif
;; b-r 
max_br_err = max( bin_br_err )
if ( max_br_err > 0.25 ) then begin 
    max_br = ( max_br + ( max_br_err / 2.0 ) + 0.02 ) 
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; tick size 
;; b-v
range = ( max_bv - min_bv ) 
if ( range le 0.5 ) then begin 
    bv_tick = 0.1 
endif else begin 
    if ( range lt 0.9 ) then begin 
        bv_tick = 0.2 
    endif else begin 
        if ( range lt 1.2 ) then begin 
            bv_tick = 0.3 
        endif else begin 
            bv_tick = 0.4 
        endelse 
    endelse
endelse
;; v-r
range = ( max_vr - min_vr ) 
if ( range le 0.5 ) then begin 
    vr_tick = 0.1 
endif else begin 
    if ( range lt 0.9 ) then begin 
        vr_tick = 0.2 
    endif else begin 
        if ( range lt 1.2 ) then begin 
            vr_tick = 0.3 
        endif else begin 
            vr_tick = 0.4 
        endelse 
    endelse
endelse
;; b-r
range = ( max_br - min_br ) 
if ( range le 0.5 ) then begin 
    br_tick = 0.1 
endif else begin 
    if ( range lt 0.9 ) then begin 
        br_tick = 0.2 
    endif else begin 
        if ( range lt 1.2 ) then begin 
            br_tick = 0.3 
        endif else begin 
            br_tick = 0.4 
        endelse 
    endelse
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
bin_m = 0.1 
bin_c = 0.005 
bin_c2 = 0.0025
;; b-v color 
bv_ori_cmd = hist_2d( bv_ori_color_use, v_mu_ori, bin1=bin_c, bin2=bin_m, $
    min1=min_bv, min2=min_mu, max1=max_bv, max2=max_mu )
bv_mod_cmd = hist_2d( bv_mod_color_use, v_mu_mod, bin1=bin_c, bin2=bin_m, $
    min1=min_bv, min2=min_mu, max1=max_bv, max2=max_mu )
;; v-r color 
vr_ori_cmd = hist_2d( vr_ori_color_use, v_mu_ori, bin1=bin_c, bin2=bin_m, $
    min1=min_vr, min2=min_mu, max1=max_vr, max2=max_mu )
vr_mod_cmd = hist_2d( vr_mod_color_use, v_mu_mod, bin1=bin_c, bin2=bin_m, $
    min1=min_vr, min2=min_mu, max1=max_vr, max2=max_mu )
;; b-r color 
br_ori_cmd = hist_2d( br_ori_color_use, v_mu_ori, bin1=bin_c, bin2=bin_m, $
    min1=min_br, min2=min_mu, max1=max_br, max2=max_mu )
br_mod_cmd = hist_2d( br_mod_color_use, v_mu_mod, bin1=bin_c, bin2=bin_m, $
    min1=min_br, min2=min_mu, max1=max_br, max2=max_mu )
;; bv-vr ccd 
bv_vr_ori = hist_2d( bv_ori_color_use, vr_ori_color_use, bin1=bin_c, $
    bin2=bin_c2, min1=min_bv_r80, min2=min_vr_r80, max1=max_bv_r80, $
    max2=max_vr_r80 )
bv_vr_mod = hist_2d( bv_mod_color_use, vr_mod_color_use, bin1=bin_c, $
    bin2=bin_c2, min1=min_bv_r80, min2=min_vr_r80, max1=max_bv_r80, $
    max2=max_vr_r80 )
;; bv-br ccd 
bv_br_ori = hist_2d( bv_ori_color_use, br_ori_color_use, bin1=bin_c, $
    bin2=bin_c, min1=min_bv_r80, min2=min_br_r80, max1=max_bv_r80, $
    max2=max_br_r80 )
bv_br_mod = hist_2d( bv_mod_color_use, br_mod_color_use, bin1=bin_c, $
    bin2=bin_c, min1=min_bv_r80, min2=min_br_r80, max1=max_bv_r80, $
    max2=max_br_r80 )
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; start plotting !! 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; color array for components 
if ( ncomp eq 3 ) then begin 
    cname_arr = [ 'RED', 'GREEN', 'BLUE' ]
endif else begin 
    if ( ncomp eq 4 ) then begin 
        cname_arr = [ 'RED', 'GREEN', 'ORANGE', 'BLUE' ] 
    endif else begin 
        cname_arr = [ 'RED', 'GREEN', 'ORANGE', 'BLUE', 'YELLOW', 'BROWN', $
            'PURPLE', 'SALMON', 'GRAY' ]
    endelse
endelse 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; bv cmd map 
plot_bv_cmd = galaxy + '_' + standard_filter + '_' + ncomp_string + '_' + $
    type + '_B-V_cmd.ps'  
xrange = [ min_bv, max_bv ]
yrange = [ min_mu, max_mu ]
xtitle = 'B-V Color'
ytitle = textoidl( '\mu_{V} (mag/arcsec^2)' ) 
zero = make_array( nbins, VALUE=0.0 )
mydevice = !D.Name 
set_plot, 'PS' 
device, filename=plot_bv_cmd, font_size=7, /encapsul, $
    /color, set_font='HELVETICA BOLD', /tt_font, xsize=24, ysize=15
;; 1 original image
cgPlot, bv_ori_color_use, v_mu_ori, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, ytitle=ytitle, xtickinterval=bv_tick, $
    ytickinterval=2.0, charsize=2.2, charthick=4.0, xthick=4.3, ythick=4.3, $ 
    /noerase, /nodata, position=position_1
;; scale the cmd 
sky, bv_ori_cmd, bg_cmd, sig_cmd, /quiet 
min_cmd = bg_cmd + 0.5 * sig_cmd 
max_cmd = max( bv_ori_cmd ) - sig_cmd 
f_cmd = asinh( ( max_cmd - min_cmd ) / sig_cmd ) 
s_cmd = asinh( ( bv_ori_cmd - min_cmd ) / sig_cmd ) / f_cmd 
bs_cmd = 254 - bytscl( s_cmd ) 
imgunder, bs_cmd 
;; draw the line for average color 
cgPlot, [bv_mod_color_avg,bv_mod_color_avg], !Y.Crange, psym=0, linestyle=5, $ 
    thick=3.6, /overplot, color=cgColor( 'CYAN', !D.Table_Size )
;; plot the running mean and their scatters 
;; for ori 
cgPlot, bin_bv_ori, bin_v_mu, psym=0, linestyle=2, thick=4.0, /overplot, $ 
    color=cgColor( 'ORG5', !D.Table_Size ) 
cgPlot, bin_bv_ori, bin_v_mu, psym=4, symsize=1.15, thick=4.0, /overplot, $ 
    color=cgColor( 'ORG5', !D.Table_Size ) 
oploterror, bin_bv_ori, bin_v_mu, bin_v_mu_err, zero, psym=3, errthick=3.2, $ 
    errcolor=cgColor( 'ORG5', !D.Table_Size )
;; for mod 
cgPlot, bin_bv_mod, bin_v_mu, psym=0, linestyle=2, thick=4.0, /overplot, $ 
    color=cgColor( 'YELLOW', !D.Table_Size ) 
cgPlot, bin_bv_mod, bin_v_mu, psym=4, symsize=1.15, thick=4.0, /overplot, $ 
    color=cgColor( 'YELLOW', !D.Table_Size ) 
;; components
for i = 0, ( ncomp - 1 ), 1 do begin 
    mu = v_comp_mue[i] 
    c = bv_comp_color[i]
    cgPlots, c, mu, psym=symbols[i], symsize=comp_symsize[i], $ 
        color=cgColor( cname_arr[i], !D.Table_Size ), /data
endfor
;; put legend and galaxy name
cgText, 0.12, 0.18, 'DATA', charsize=2.4, charthick=4.0, $
    color=cgColor( 'BLACK', !D.Table_Size ), /normal

cgPlots, [ 0.115, 0.17 ], [ 0.280, 0.280 ], psym=0, linestyle=2, thick=3.6, $
    color=cgColor( 'ORG5', !D.Table_Size ), /normal
cgText, 0.178, 0.275, 'data', charsize=1.9, charthick=3.5, $
    alignment=0, /normal 
cgPlots, [ 0.115, 0.17 ], [ 0.250, 0.250 ], psym=0, linestyle=3, thick=3.6, $
    color=cgColor( 'YELLOW', !D.Table_Size ), /normal
cgText, 0.178, 0.245, 'model', charsize=1.9, charthick=3.5, $
    alignment=0, /normal, color=cgColor( 'BLACK', !D.Table_Size ) 
if ( strpos( galaxy, 'ESO' ) ne -1 ) then begin 
    charsize = 2.0
endif else begin 
    charsize = 2.6
endelse
cgText, 0.115, 0.340, galaxy, charsize=charsize, charthick=4.2, /normal, $ 
    color=cgColor( 'BLACK', !D.Table_Size )

cgPlot, bv_ori_color_use, v_mu_ori, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, ytitle=ytitle, xtickinterval=bv_tick, $
    ytickinterval=2.0, charsize=2.2, charthick=4.0, xthick=4.3, ythick=4.3, $ 
    /noerase, /nodata, position=position_1

;; 2 model image
cgPlot, bv_mod_color_use, v_mu_mod, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, xtickinterval=bv_tick, $
    ytickinterval=2.0, charsize=2.2, charthick=4.0, xthick=4.3, ythick=4.3, $ 
    /noerase, /nodata, position=position_2, ytickformat='(A1)'
;; scale the cmd 
sky, bv_mod_cmd, bg_cmd, sig_cmd, /quiet 
min_cmd = bg_cmd + 0.5 * sig_cmd 
max_cmd = max( bv_mod_cmd ) - sig_cmd 
f_cmd = asinh( ( max_cmd - min_cmd ) / sig_cmd ) 
s_cmd = asinh( ( bv_mod_cmd - min_cmd ) / sig_cmd ) / f_cmd 
bs_cmd = 254 - bytscl( s_cmd ) 
imgunder, bs_cmd 
;; draw the line for average color 
cgPlot, [bv_mod_color_avg, bv_mod_color_avg], !Y.Crange, psym=0, linestyle=5, $ 
    thick=3.6, /overplot, color=cgColor( 'CYAN', !D.Table_Size )

;; components
for i = 0, ( ncomp - 1 ), 1 do begin 
    mu = v_comp_mue[i] 
    c = bv_comp_color[i]
    cgPlots, c, mu, psym=symbols[i], symsize=comp_symsize[i], $ 
        color=cgColor( cname_arr[i], !D.Table_Size ), /data
endfor
;; legend and information
cgText, 0.520, 0.18, 'MODEL', charsize=2.4, charthick=4.0, $
    color=cgColor( 'BLACK', !D.Table_Size ), /normal
cgText, 0.525, 0.24, type, charsize=2.1, charthick=3.6, /normal, $ 
    color=cgColor( 'BLACK', !D.Table_Size )
cgPlots, [ 0.515, 0.568 ], [ 0.290, 0.290 ], psym=0, linestyle=5, $
    thick=3.6, color=cgColor( 'CYAN', !D.Table_Size ), /normal
temp = '<B-V>'
cgText, 0.576, 0.285, temp, charsize=1.8, charthick=3.5, $
    alignment=0, /normal, color=cgColor( 'BLACK', !D.Table_Size ) 
;; extinction arrow
x0 = xrange[0] + ( ( xrange[1] - xrange[0] ) / 13 ) 
y0 = yrange[1] - ( ( yrange[1] - yrange[0] ) / 2.2 )
diff_extinc = ( ( extinc_b / extinc_v ) - 1.0 ) * 0.25
x1 = x0 + diff_extinc
y1 = y0 + 0.2 
x2 = ( ( x0 + x1 ) / 2.0 ) * 1.01
y2 = ( y0 - 0.4  )
cgArrow, x0, y0, x1, y1,  /data, hthick=5.0, thick=4.0, $
    color=cgColor( 'RED', !D.Table_Size ), hsize=5.5
ext_string = textoidl( 'A_V=0.2' )
cgText, x2, y2, ext_string,  charsize=1.7, charthick=3.4, $
    color=cgColor( 'BLACK', !D.Table_Size ), /data, alignment=0.5
;; error ellipse 
offset = ( ( xrange[1] - xrange[0] ) / 24 ) 
if ( offset ge 0.06 ) then begin 
    offset = 0.005
endif 
max_cerr = max( bin_bv_err / 2.0 ) + offset 
xc = xrange[1] - max_cerr 
for i = 0, ( nbins - 2 ), 1 do begin 
    ra = ( bin_bv_err[i] / 2.0 )
    rb = ( bin_v_mu_err[i] / 2.0 )
    xc = xc 
    yc = bin_v_mu[i] 
    pos_ang = 0.0 
    tvellipse, ra, rb, xc, yc, pos_ang, /data, $ 
        color=cgColor( 'BLACK', !D.Table_Size ), thick=2.8 
    tvellipse, ra, rb, xc, yc, pos_ang, /data, $ 
        color=cgColor( 'GRAY', !D.Table_Size ), /fill 
endfor

cgPlot, bv_mod_color_use, v_mu_mod, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, xtickinterval=bv_tick, $
    ytickinterval=2.0, charsize=2.2, charthick=4.0, xthick=4.3, ythick=4.3, $ 
    /noerase, /nodata, position=position_2, ytickformat='(A1)'
;; close ps 
device, /close
set_plot, mydevice
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; vr cmd map 
plot_vr_cmd = galaxy + '_' + standard_filter + '_' + ncomp_string + '_' + $
    type + '_V-R_cmd.ps'  
xrange = [ min_vr, max_vr ]
yrange = [ min_mu, max_mu ]
xtitle = 'V-R Color'
ytitle = textoidl( '\mu_{V} (mag/arcsec^2)' ) 
zero = make_array( nbins, VALUE=0.0 )
mydevice = !D.Name 
set_plot, 'PS' 
device, filename=plot_vr_cmd, font_size=7, /encapsul, $
    /color, set_font='HELVETICA BOLD', /tt_font, xsize=24, ysize=15
;; 1 original image
cgPlot, vr_ori_color_use, v_mu_ori, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, ytitle=ytitle, xtickinterval=vr_tick, $
    ytickinterval=2.0, charsize=2.2, charthick=4.0, xthick=4.3, ythick=4.3, $ 
    /noerase, /nodata, position=position_1
;; scale the cmd 
sky, vr_ori_cmd, bg_cmd, sig_cmd, /quiet 
min_cmd = bg_cmd + 0.5 * sig_cmd 
max_cmd = max( vr_ori_cmd ) - sig_cmd 
f_cmd = asinh( ( max_cmd - min_cmd ) / sig_cmd ) 
s_cmd = asinh( ( vr_ori_cmd - min_cmd ) / sig_cmd ) / f_cmd 
bs_cmd = 254 - bytscl( s_cmd ) 
imgunder, bs_cmd 
;; draw the line for average color 
cgPlot, [vr_mod_color_avg,vr_mod_color_avg], !Y.Crange, psym=0, linestyle=5, $ 
    thick=3.6, /overplot, color=cgColor( 'CYAN', !D.Table_Size )
;; plot the running mean and their scatters 
;; for ori 
cgPlot, bin_vr_ori, bin_v_mu, psym=0, linestyle=2, thick=4.0, /overplot, $ 
    color=cgColor( 'ORG5', !D.Table_Size ) 
cgPlot, bin_vr_ori, bin_v_mu, psym=4, symsize=1.15, thick=4.0, /overplot, $ 
    color=cgColor( 'ORG5', !D.Table_Size ) 
oploterror, bin_vr_ori, bin_v_mu, bin_v_mu_err, zero, psym=3, errthick=3.2, $ 
    errcolor=cgColor( 'ORG5', !D.Table_Size )
;; for mod 
cgPlot, bin_vr_mod, bin_v_mu, psym=0, linestyle=2, thick=4.0, /overplot, $ 
    color=cgColor( 'YELLOW', !D.Table_Size ) 
cgPlot, bin_vr_mod, bin_v_mu, psym=4, symsize=1.15, thick=4.0, /overplot, $ 
    color=cgColor( 'YELLOW', !D.Table_Size ) 
;; components
for i = 0, ( ncomp - 1 ), 1 do begin 
    mu = v_comp_mue[i] 
    c = vr_comp_color[i]
    cgPlots, c, mu, psym=symbols[i], symsize=comp_symsize[i], $ 
        color=cgColor( cname_arr[i], !D.Table_Size ), /data
endfor
;; put legend and galaxy name
cgText, 0.12, 0.18, 'DATA', charsize=2.4, charthick=4.0, $
    color=cgColor( 'BLACK', !D.Table_Size ), /normal

cgPlots, [ 0.115, 0.17 ], [ 0.280, 0.280 ], psym=0, linestyle=2, thick=3.6, $
    color=cgColor( 'ORG5', !D.Table_Size ), /normal
cgText, 0.178, 0.275, 'data', charsize=1.9, charthick=3.5, $
    alignment=0, /normal 
cgPlots, [ 0.115, 0.17 ], [ 0.250, 0.250 ], psym=0, linestyle=3, thick=3.6, $
    color=cgColor( 'YELLOW', !D.Table_Size ), /normal
cgText, 0.178, 0.245, 'model', charsize=1.9, charthick=3.5, $
    alignment=0, /normal, color=cgColor( 'BLACK', !D.Table_Size ) 
if ( strpos( galaxy, 'ESO' ) ne -1 ) then begin 
    charsize = 2.0
endif else begin 
    charsize = 2.6
endelse
cgText, 0.115, 0.340, galaxy, charsize=charsize, charthick=4.2, /normal, $ 
    color=cgColor( 'BLACK', !D.Table_Size )

cgPlot, vr_ori_color_use, v_mu_ori, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, ytitle=ytitle, xtickinterval=vr_tick, $
    ytickinterval=2.0, charsize=2.2, charthick=4.0, xthick=4.3, ythick=4.3, $ 
    /noerase, /nodata, position=position_1

;; 2 model image
cgPlot, vr_mod_color_use, v_mu_mod, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, xtickinterval=vr_tick, $
    ytickinterval=2.0, charsize=2.2, charthick=4.0, xthick=4.3, ythick=4.3, $ 
    /noerase, /nodata, position=position_2, ytickformat='(A1)'
;; scale the cmd 
sky, vr_mod_cmd, bg_cmd, sig_cmd, /quiet 
min_cmd = bg_cmd + 0.5 * sig_cmd 
max_cmd = max( vr_mod_cmd ) - sig_cmd 
f_cmd = asinh( ( max_cmd - min_cmd ) / sig_cmd ) 
s_cmd = asinh( ( vr_mod_cmd - min_cmd ) / sig_cmd ) / f_cmd 
bs_cmd = 254 - bytscl( s_cmd ) 
imgunder, bs_cmd 
;; draw the line for average color 
cgPlot, [vr_mod_color_avg, vr_mod_color_avg], !Y.Crange, psym=0, linestyle=5, $ 
    thick=3.6, /overplot, color=cgColor( 'CYAN', !D.Table_Size )

;; components
for i = 0, ( ncomp - 1 ), 1 do begin 
    mu = v_comp_mue[i] 
    c = vr_comp_color[i]
    cgPlots, c, mu, psym=symbols[i], symsize=comp_symsize[i], $ 
        color=cgColor( cname_arr[i], !D.Table_Size ), /data
endfor
;; legend and information
cgText, 0.520, 0.18, 'MODEL', charsize=2.4, charthick=4.0, $
    color=cgColor( 'BLACK', !D.Table_Size ), /normal
cgText, 0.525, 0.24, type, charsize=2.1, charthick=3.6, /normal, $ 
    color=cgColor( 'BLACK', !D.Table_Size )
cgPlots, [ 0.515, 0.568 ], [ 0.290, 0.290 ], psym=0, linestyle=5, $
    thick=3.6, color=cgColor( 'CYAN', !D.Table_Size ), /normal
temp = '<V-R>'
cgText, 0.576, 0.285, temp, charsize=1.8, charthick=3.5, $
    alignment=0, /normal, color=cgColor( 'BLACK', !D.Table_Size ) 
;; extinction arrow
x0 = xrange[0] + ( ( xrange[1] - xrange[0] ) / 13 ) 
y0 = yrange[1] - ( ( yrange[1] - yrange[0] ) / 2.2 )
diff_extinc = ( ( extinc_v / extinc_r ) - 1.0 ) * 0.25
x1 = x0 + diff_extinc
y1 = y0 + 0.2 
x2 = ( ( x0 + x1 ) / 2.0 ) * 1.01
y2 = ( y0 - 0.4  )
cgArrow, x0, y0, x1, y1,  /data, hthick=5.0, thick=4.0, $
    color=cgColor( 'RED', !D.Table_Size ), hsize=5.5
ext_string = textoidl( 'A_V=0.2' )
cgText, x2, y2, ext_string,  charsize=1.7, charthick=3.4, $
    color=cgColor( 'BLACK', !D.Table_Size ), /data, alignment=0.5
;; error ellipse 
offset = ( ( xrange[1] - xrange[0] ) / 24 ) 
if ( offset ge 0.06 ) then begin 
    offset = 0.005
endif 
max_cerr = max( bin_vr_err / 2.0 ) + offset 
xc = xrange[1] - max_cerr 
for i = 0, ( nbins - 2 ), 1 do begin 
    ra = ( bin_vr_err[i] / 2.0 )
    rb = ( bin_v_mu_err[i] / 2.0 )
    xc = xc 
    yc = bin_v_mu[i] 
    pos_ang = 0.0 
    tvellipse, ra, rb, xc, yc, pos_ang, /data, $ 
        color=cgColor( 'BLACK', !D.Table_Size ), thick=2.8 
    tvellipse, ra, rb, xc, yc, pos_ang, /data, $ 
        color=cgColor( 'GRAY', !D.Table_Size ), /fill 
endfor

cgPlot, vr_mod_color_use, v_mu_mod, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, xtickinterval=vr_tick, $
    ytickinterval=2.0, charsize=2.2, charthick=4.0, xthick=4.3, ythick=4.3, $ 
    /noerase, /nodata, position=position_2, ytickformat='(A1)'
;; close ps 
device, /close
set_plot, mydevice
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; br cmd map 
plot_br_cmd = galaxy + '_' + standard_filter + '_' + ncomp_string + '_' + $
    type + '_B-R_cmd.ps'  
xrange = [ min_br, max_br ]
yrange = [ min_mu, max_mu ]
xtitle = 'B-R Color'
ytitle = textoidl( '\mu_{V} (mag/arcsec^2)' ) 
zero = make_array( nbins, VALUE=0.0 )
mydevice = !D.Name 
set_plot, 'PS' 
device, filename=plot_br_cmd, font_size=7, /encapsul, $
    /color, set_font='HELVETICA BOLD', /tt_font, xsize=24, ysize=15
;; 1 original image
cgPlot, br_ori_color_use, v_mu_ori, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, ytitle=ytitle, xtickinterval=br_tick, $
    ytickinterval=2.0, charsize=2.2, charthick=4.0, xthick=4.3, ythick=4.3, $ 
    /noerase, /nodata, position=position_1
;; scale the cmd 
sky, br_ori_cmd, bg_cmd, sig_cmd, /quiet 
min_cmd = bg_cmd + 0.5 * sig_cmd 
max_cmd = max( br_ori_cmd ) - sig_cmd 
f_cmd = asinh( ( max_cmd - min_cmd ) / sig_cmd ) 
s_cmd = asinh( ( br_ori_cmd - min_cmd ) / sig_cmd ) / f_cmd 
bs_cmd = 254 - bytscl( s_cmd ) 
imgunder, bs_cmd 
;; draw the line for average color 
cgPlot, [br_mod_color_avg,br_mod_color_avg], !Y.Crange, psym=0, linestyle=5, $ 
    thick=3.6, /overplot, color=cgColor( 'CYAN', !D.Table_Size )
;; plot the running mean and their scatters 
;; for ori 
cgPlot, bin_br_ori, bin_v_mu, psym=0, linestyle=2, thick=4.0, /overplot, $ 
    color=cgColor( 'ORG5', !D.Table_Size ) 
cgPlot, bin_br_ori, bin_v_mu, psym=4, symsize=1.15, thick=4.0, /overplot, $ 
    color=cgColor( 'ORG5', !D.Table_Size ) 
oploterror, bin_br_ori, bin_v_mu, bin_v_mu_err, zero, psym=3, errthick=3.2, $ 
    errcolor=cgColor( 'ORG5', !D.Table_Size )
;; for mod 
cgPlot, bin_br_mod, bin_v_mu, psym=0, linestyle=2, thick=4.0, /overplot, $ 
    color=cgColor( 'YELLOW', !D.Table_Size ) 
cgPlot, bin_br_mod, bin_v_mu, psym=4, symsize=1.15, thick=4.0, /overplot, $ 
    color=cgColor( 'YELLOW', !D.Table_Size ) 
;; components
for i = 0, ( ncomp - 1 ), 1 do begin 
    mu = v_comp_mue[i] 
    c = br_comp_color[i]
    cgPlots, c, mu, psym=symbols[i], symsize=comp_symsize[i], $ 
        color=cgColor( cname_arr[i], !D.Table_Size ), /data
endfor
;; put legend and galaxy name
cgText, 0.12, 0.18, 'DATA', charsize=2.4, charthick=4.0, $
    color=cgColor( 'BLACK', !D.Table_Size ), /normal

cgPlots, [ 0.115, 0.17 ], [ 0.280, 0.280 ], psym=0, linestyle=2, thick=3.6, $
    color=cgColor( 'ORG5', !D.Table_Size ), /normal
cgText, 0.178, 0.275, 'data', charsize=1.9, charthick=3.5, $
    alignment=0, /normal 
cgPlots, [ 0.115, 0.17 ], [ 0.250, 0.250 ], psym=0, linestyle=3, thick=3.6, $
    color=cgColor( 'YELLOW', !D.Table_Size ), /normal
cgText, 0.178, 0.245, 'model', charsize=1.9, charthick=3.5, $
    alignment=0, /normal, color=cgColor( 'BLACK', !D.Table_Size ) 
if ( strpos( galaxy, 'ESO' ) ne -1 ) then begin 
    charsize = 2.0
endif else begin 
    charsize = 2.6
endelse
cgText, 0.115, 0.340, galaxy, charsize=charsize, charthick=4.2, /normal, $ 
    color=cgColor( 'BLACK', !D.Table_Size )

cgPlot, br_ori_color_use, v_mu_ori, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, ytitle=ytitle, xtickinterval=br_tick, $
    ytickinterval=2.0, charsize=2.2, charthick=4.0, xthick=4.3, ythick=4.3, $ 
    /noerase, /nodata, position=position_1

;; 2 model image
cgPlot, br_mod_color_use, v_mu_mod, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, xtickinterval=br_tick, $
    ytickinterval=2.0, charsize=2.2, charthick=4.0, xthick=4.3, ythick=4.3, $ 
    /noerase, /nodata, position=position_2, ytickformat='(A1)'
;; scale the cmd 
sky, br_mod_cmd, bg_cmd, sig_cmd, /quiet 
min_cmd = bg_cmd + 0.5 * sig_cmd 
max_cmd = max( br_mod_cmd ) - sig_cmd 
f_cmd = asinh( ( max_cmd - min_cmd ) / sig_cmd ) 
s_cmd = asinh( ( br_mod_cmd - min_cmd ) / sig_cmd ) / f_cmd 
bs_cmd = 254 - bytscl( s_cmd ) 
imgunder, bs_cmd 
;; draw the line for average color 
cgPlot, [br_mod_color_avg, br_mod_color_avg], !Y.Crange, psym=0, linestyle=5, $ 
    thick=3.6, /overplot, color=cgColor( 'CYAN', !D.Table_Size )

;; components
for i = 0, ( ncomp - 1 ), 1 do begin 
    mu = v_comp_mue[i] 
    c = br_comp_color[i]
    cgPlots, c, mu, psym=symbols[i], symsize=comp_symsize[i], $ 
        color=cgColor( cname_arr[i], !D.Table_Size ), /data
endfor
;; legend and information
cgText, 0.520, 0.18, 'MODEL', charsize=2.4, charthick=4.0, $
    color=cgColor( 'BLACK', !D.Table_Size ), /normal
cgText, 0.525, 0.24, type, charsize=2.1, charthick=3.6, /normal, $ 
    color=cgColor( 'BLACK', !D.Table_Size )
cgPlots, [ 0.515, 0.568 ], [ 0.290, 0.290 ], psym=0, linestyle=5, $
    thick=3.6, color=cgColor( 'CYAN', !D.Table_Size ), /normal
temp = '<B-R>'
cgText, 0.576, 0.285, temp, charsize=1.8, charthick=3.5, $
    alignment=0, /normal, color=cgColor( 'BLACK', !D.Table_Size ) 
;; extinction arrow
x0 = xrange[0] + ( ( xrange[1] - xrange[0] ) / 13 ) 
y0 = yrange[1] - ( ( yrange[1] - yrange[0] ) / 2.2 )
diff_extinc = ( ( extinc_b / extinc_r ) - 1.0 ) * 0.25
x1 = x0 + diff_extinc
y1 = y0 + 0.2 
x2 = ( ( x0 + x1 ) / 2.0 ) * 1.01
y2 = ( y0 - 0.4  )
cgArrow, x0, y0, x1, y1,  /data, hthick=5.0, thick=4.0, $
    color=cgColor( 'RED', !D.Table_Size ), hsize=5.5
ext_string = textoidl( 'A_V=0.2' )
cgText, x2, y2, ext_string,  charsize=1.7, charthick=3.4, $
    color=cgColor( 'BLACK', !D.Table_Size ), /data, alignment=0.5
;; error ellipse 
offset = ( ( xrange[1] - xrange[0] ) / 24 ) 
if ( offset ge 0.06 ) then begin 
    offset = 0.005
endif 
max_cerr = max( bin_br_err / 2.0 ) + offset 
xc = xrange[1] - max_cerr 
for i = 0, ( nbins - 2 ), 1 do begin 
    ra = ( bin_br_err[i] / 2.0 )
    rb = ( bin_v_mu_err[i] / 2.0 )
    xc = xc 
    yc = bin_v_mu[i] 
    pos_ang = 0.0 
    tvellipse, ra, rb, xc, yc, pos_ang, /data, $ 
        color=cgColor( 'BLACK', !D.Table_Size ), thick=2.8 
    tvellipse, ra, rb, xc, yc, pos_ang, /data, $ 
        color=cgColor( 'GRAY', !D.Table_Size ), /fill 
endfor

cgPlot, br_mod_color_use, v_mu_mod, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, xtickinterval=br_tick, $
    ytickinterval=2.0, charsize=2.2, charthick=4.0, xthick=4.3, ythick=4.3, $ 
    /noerase, /nodata, position=position_2, ytickformat='(A1)'
;; close ps 
device, /close
set_plot, mydevice
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; bv vr color-color map 
plot_bv_vr = galaxy + '_' + standard_filter + '_' + ncomp_string + '_' + $
    type + '_BV-VR.ps'  
xrange = [ min_bv_r80, max_bv_r80 ]
yrange = [ min_vr_r80, max_vr_r80 ]
xtitle = 'B-V Color'
ytitle = 'V-R Color' 
zero = make_array( nbins, VALUE=0.0 )
mydevice = !D.Name 
set_plot, 'PS' 
device, filename=plot_bv_vr, font_size=7, /encapsul, $
    /color, set_font='HELVETICA BOLD', /tt_font, xsize=24, ysize=15
;; 1 original image
cgPlot, bv_ori_color_use, vr_ori_color_use, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, ytitle=ytitle, charsize=2.2, charthick=4.0, $
    xthick=4.3, ythick=4.3, /noerase, /nodata, position=position_1
;; scale the cmd 
sky, bv_vr_ori, bg_cmd, sig_cmd, /quiet 
min_cmd = bg_cmd + 0.5 * sig_cmd 
max_cmd = max( bv_vr_ori ) - sig_cmd 
f_cmd = asinh( ( max_cmd - min_cmd ) / sig_cmd ) 
s_cmd = asinh( ( bv_vr_ori - min_cmd ) / sig_cmd ) / f_cmd 
bs_cmd = 254 - bytscl( s_cmd ) 
imgunder, bs_cmd 
;; components
for i = 0, ( ncomp - 1 ), 1 do begin 
    c1 = bv_comp_color[i] 
    c2 = vr_comp_color[i]
    cgPlots, c1, c2, psym=symbols[i], symsize=comp_symsize[i], $ 
        color=cgColor( cname_arr[i], !D.Table_Size ), /data
endfor
;; galaxy name 
cgText, 0.130, 0.220, galaxy, charsize=3.4, charthick=5.8, /normal, $ 
    color=cgColor( 'BLACK', !D.Table_Size )
;; draw the axis again
cgPlot, bv_ori_color_use, vr_ori_color_use, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, ytitle=ytitle, charsize=2.2, charthick=4.0, $
    xthick=4.3, ythick=4.3, /noerase, /nodata, position=position_1

;; 2 model image
cgPlot, bv_mod_color_use, vr_mod_color_use, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, charsize=2.2, charthick=4.0, xthick=4.3, $
    ythick=4.3, /noerase, /nodata, position=position_2, ytickformat='(A1)'
;; scale the cmd 
sky, bv_vr_mod, bg_cmd, sig_cmd, /quiet 
min_cmd = bg_cmd + 0.5 * sig_cmd 
max_cmd = max( bv_vr_mod ) - sig_cmd 
f_cmd = asinh( ( max_cmd - min_cmd ) / sig_cmd ) 
s_cmd = asinh( ( bv_vr_mod - min_cmd ) / sig_cmd ) / f_cmd 
bs_cmd = 254 - bytscl( s_cmd ) 
imgunder, bs_cmd 
;; components
for i = 0, ( ncomp - 1 ), 1 do begin 
    c1 = bv_comp_color[i] 
    c2 = vr_comp_color[i]
    cgPlots, c1, c2, psym=symbols[i], symsize=comp_symsize[i], $ 
        color=cgColor( cname_arr[i], !D.Table_Size ), /data
endfor
;; overplot the stellar population grid 
cgPlot, age1_bv, age1_vr, psym=16, symsize=0.9, thick=3.0, /overplot, $
    color=cgColor( 'CYAN', !D.Table_Size )
cgPlot, age1_bv, age1_vr, psym=0, linestyle=2, thick=3.0, /overplot, $
    color=cgColor( 'CYAN', !D.Table_Size )
cgPlot, age2_bv, age2_vr, psym=16, symsize=0.9, thick=3.0, /overplot, $
    color=cgColor( 'ORANGE', !D.Table_Size )
cgPlot, age2_bv, age2_vr, psym=0, linestyle=2, thick=3.0, /overplot, $
    color=cgColor( 'ORANGE', !D.Table_Size )
cgPlot, age3_bv, age3_vr, psym=16, symsize=0.9, thick=3.0, /overplot, $
    color=cgColor( 'GOLD', !D.Table_Size )
cgPlot, age3_bv, age3_vr, psym=0, linestyle=2, thick=3.0, /overplot, $
    color=cgColor( 'GOLD', !D.Table_Size )
cgPlot, age4_bv, age3_vr, psym=16, symsize=0.9, thick=3.0, /overplot, $
    color=cgColor( 'MAGENTA', !D.Table_Size )
cgPlot, age4_bv, age3_vr, psym=0, linestyle=2, thick=3.0, /overplot, $
    color=cgColor( 'MAGENTA', !D.Table_Size )
;; extinction arrow
x0 = xrange[0] + ( ( xrange[1] - xrange[0] ) / 13 ) 
y0 = yrange[1] - ( ( yrange[1] - yrange[0] ) / 4.0 )
diff_extinc1 = ( ( extinc_b / extinc_v ) - 1.0 ) * 0.05
diff_extinc2 = ( 1.0 - ( extinc_r / extinc_b ) ) * 0.05 
x1 = x0 + diff_extinc1
y1 = y0 + diff_extinc2 
x2 = ( ( x0 + x1 ) / 2.0 ) * 1.01
y2 = ( y0 - ( yrange[1] - yrange[0] ) / 30 )
cgArrow, x0, y0, x1, y1,  /data, hthick=5.0, thick=4.0, $
    color=cgColor( 'RED', !D.Table_Size ), hsize=5.5
ext_string = textoidl( 'A_V=0.05' )
cgText, x2, y2, ext_string,  charsize=1.7, charthick=3.4, $
    color=cgColor( 'BLACK', !D.Table_Size ), /data, alignment=0.5

;; model type 
cgText, 0.530, 0.220, type, charsize=2.9, charthick=4.6, /normal, $ 
    color=cgColor( 'BLACK', !D.Table_Size )
;; draw the axis again
cgPlot, bv_mod_color_use, vr_mod_color_use, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, charsize=2.2, charthick=4.0, xthick=4.3, $
    ythick=4.3, /noerase, /nodata, position=position_2, ytickformat='(A1)'

;; close ps 
device, /close
set_plot, mydevice
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; bv br color-color map 
plot_bv_br = galaxy + '_' + standard_filter + '_' + ncomp_string + '_' + $
    type + '_BV-BR.ps'  
xrange = [ min_bv_r80, max_bv_r80 ]
yrange = [ min_br_r80, max_br_r80 ]
xtitle = 'B-V Color'
ytitle = 'B-R Color' 
zero = make_array( nbins, VALUE=0.0 )
mydevice = !D.Name 
set_plot, 'PS' 
device, filename=plot_bv_br, font_size=7, /encapsul, $
    /color, set_font='HELVETICA BOLD', /tt_font, xsize=24, ysize=15
;; 1 original image
cgPlot, bv_ori_color_use, br_ori_color_use, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, ytitle=ytitle, charsize=2.2, charthick=4.0, $
    xthick=4.3, ythick=4.3, /noerase, /nodata, position=position_1
;; scale the cmd 
sky, bv_br_ori, bg_cmd, sig_cmd, /quiet 
min_cmd = bg_cmd + 0.5 * sig_cmd 
max_cmd = max( bv_br_ori ) - sig_cmd 
f_cmd = asinh( ( max_cmd - min_cmd ) / sig_cmd ) 
s_cmd = asinh( ( bv_br_ori - min_cmd ) / sig_cmd ) / f_cmd 
bs_cmd = 254 - bytscl( s_cmd ) 
imgunder, bs_cmd 
;; components
for i = 0, ( ncomp - 1 ), 1 do begin 
    c1 = bv_comp_color[i] 
    c2 = br_comp_color[i]
    cgPlots, c1, c2, psym=symbols[i], symsize=comp_symsize[i], $ 
        color=cgColor( cname_arr[i], !D.Table_Size ), /data
endfor
;; galaxy name 
cgText, 0.130, 0.220, galaxy, charsize=3.4, charthick=5.8, /normal, $ 
    color=cgColor( 'BLACK', !D.Table_Size )
;; draw the axis again
cgPlot, bv_ori_color_use, br_ori_color_use, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, ytitle=ytitle, charsize=2.2, charthick=4.0, $
    xthick=4.3, ythick=4.3, /noerase, /nodata, position=position_1

;; 2 model image
cgPlot, bv_mod_color_use, br_mod_color_use, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, charsize=2.2, charthick=4.0, xthick=4.3, $
    ythick=4.3, /noerase, /nodata, position=position_2, ytickformat='(A1)'
;; scale the cmd 
sky, bv_br_mod, bg_cmd, sig_cmd, /quiet 
min_cmd = bg_cmd + 0.5 * sig_cmd 
max_cmd = max( bv_br_mod ) - sig_cmd 
f_cmd = asinh( ( max_cmd - min_cmd ) / sig_cmd ) 
s_cmd = asinh( ( bv_br_mod - min_cmd ) / sig_cmd ) / f_cmd 
bs_cmd = 254 - bytscl( s_cmd ) 
imgunder, bs_cmd 
;; components
for i = 0, ( ncomp - 1 ), 1 do begin 
    c1 = bv_comp_color[i] 
    c2 = br_comp_color[i]
    cgPlots, c1, c2, psym=symbols[i], symsize=comp_symsize[i], $ 
        color=cgColor( cname_arr[i], !D.Table_Size ), /data
endfor
;; overplot the stellar population grid 
cgPlot, age1_bv, age1_br, psym=16, symsize=0.9, thick=3.0, /overplot, $
    color=cgColor( 'CYAN', !D.Table_Size )
cgPlot, age1_bv, age1_br, psym=0, linestyle=2, thick=3.0, /overplot, $
    color=cgColor( 'CYAN', !D.Table_Size )
cgPlot, age2_bv, age2_br, psym=16, symsize=0.9, thick=3.0, /overplot, $
    color=cgColor( 'ORANGE', !D.Table_Size )
cgPlot, age2_bv, age2_br, psym=0, linestyle=2, thick=3.0, /overplot, $
    color=cgColor( 'ORANGE', !D.Table_Size )
cgPlot, age3_bv, age3_br, psym=16, symsize=0.9, thick=3.0, /overplot, $
    color=cgColor( 'GOLD', !D.Table_Size )
cgPlot, age3_bv, age3_br, psym=0, linestyle=2, thick=3.0, /overplot, $
    color=cgColor( 'GOLD', !D.Table_Size )
cgPlot, age4_bv, age3_br, psym=16, symsize=0.9, thick=3.0, /overplot, $
    color=cgColor( 'MAGENTA', !D.Table_Size )
cgPlot, age4_bv, age3_br, psym=0, linestyle=2, thick=3.0, /overplot, $
    color=cgColor( 'MAGENTA', !D.Table_Size )
;; extinction arrow
x0 = xrange[0] + ( ( xrange[1] - xrange[0] ) / 13 ) 
y0 = yrange[1] - ( ( yrange[1] - yrange[0] ) / 4.0 )
diff_extinc1 = ( ( extinc_b / extinc_v ) - 1.0 ) * 0.05
diff_extinc2 = ( ( extinc_b / extinc_r ) - 1.0 ) * 0.05 * ( extinc_r / extinc_v ) 
x1 = x0 + diff_extinc1
y1 = y0 + diff_extinc2 
x2 = ( ( x0 + x1 ) / 2.0 ) * 1.01
y2 = ( y0 - ( yrange[1] - yrange[0] ) / 30 )
cgArrow, x0, y0, x1, y1,  /data, hthick=5.0, thick=4.0, $
    color=cgColor( 'RED', !D.Table_Size ), hsize=5.5
ext_string = textoidl( 'A_V=0.05' )
cgText, x2, y2, ext_string,  charsize=1.7, charthick=3.4, $
    color=cgColor( 'BLACK', !D.Table_Size ), /data, alignment=0.5

;; model type 
cgText, 0.530, 0.220, type, charsize=2.9, charthick=4.6, /normal, $ 
    color=cgColor( 'BLACK', !D.Table_Size )
;; draw the axis again
cgPlot, bv_mod_color_use, br_mod_color_use, xstyle=1, ystyle=1, xrange=xrange, $ 
    yrange=yrange, xtitle=xtitle, charsize=2.2, charthick=4.0, xthick=4.3, $
    ythick=4.3, /noerase, /nodata, position=position_2, ytickformat='(A1)'
;; close ps 
device, /close
set_plot, mydevice
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; color map bv 
plot_br_map = galaxy + '_' + standard_filter + '_' + ncomp_string + '_' + $
    type + '_BR_map.ps'  
show = size( bv_ori_color_show, /dimension )
xrange = [ 0, show[0] ]
yrange = [ 0, show[1] ] 
xtick = ( show[0] + 2 )
ytick = ( show[1] + 2 )

;; make the ps 
position_1=[0.050, 0.0, 0.50, 1.0 ]
position_2=[0.502, 0.0, 0.952, 1.0 ]
;; cgLoadCT, 31, NColors=254
;TVLCT, palette, /Get
;mydevice = !D.Name 
;set_plot, 'PS' 
;device, filename=plot_br_map, font_size=7, /encapsul, $
;    /color, set_font='HELVETICA BOLD', /tt_font, xsize=24, ysize=12
;
;cgImage, br_ori_color_show, Stretch='Clip', Clip=4, $ 
;    background='white', position=position_1, /noerase, $
;    ;Missing_Value=-99, Missing_Color='White', $
;    MinValue=min_br_r80, MaxValue=max_br_r80
;
;cgImage, br_mod_color_show, Stretch='Clip', Clip=4, $ 
;    background='white', position=position_2, /noerase, $
;    MinValue=min_br_r80, MaxValue=max_br_r80
;;; galaxy name 
;cgText, 0.520, 0.230, galaxy, charsize=4.0, charthick=4.5, /normal, $ 
;    color=cgColor( 'BLACK', !D.Table_Size )
;
;;; close ps 
;device, /close
;set_plot, mydevice
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end 
