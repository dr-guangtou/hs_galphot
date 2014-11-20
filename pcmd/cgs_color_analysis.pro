pro cgs_color_analysis, galaxy, band1, band2, ncomp, type, $ 
    color_map = color_map, standard_filter = standard_filter, $
    comp_color = comp_color, galaxy_color = galaxy_color, $ 
    no_input_binary = no_input_binary, rebin_factor = rebin_factor, $ 
    maxsma = maxsma, add_noise=add_noise 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Color analysis based on the model from band1 and band2 
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
minsma = 0.0 
if keyword_set( maxsma ) then begin 
    maxsma = float( maxsma )
endif else begin 
    maxsma = 1100
endelse
step = 0.08
sma0 = 30
expand_factor = 2.0
color_loc = '/media/Elements/cings/test/E_color/cprof/'
if keyword_set( standard_filter ) then begin 
    standard_filter = strcompress( string( standard_filter ), /remove_all )
endif else begin 
    standard_filter = 'V'
endelse
filter_array = [ 'B', 'V', 'R', 'I' ]
photoerr_array = [ 0.06, 0.04, 0.03, 0.04 ]
mtype_array = [ 'as0', 'bs0', 'cs0', 'af0', 'bf0', 'cf0' ]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
location_setting, setting
;;; define positions of files
    galfit = setting.galfit
  cat_head = setting.header
  datafile = setting.header
   fileloc = setting.workplace
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; the parameter files for ellipse:  
   ell_par = setting.ell_par
imcopy_par = setting.imcopy_par
 tdump_par = setting.tdump_par
trebin_par = setting.trebin_par
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;iraf command:
    ximage = setting.ximage 
  isophote = setting.isophote
    ttools = setting.ttools
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; prepare the input 
galaxy = strcompress( galaxy, /remove_all ) 
;; filters 
band1 = strcompress( strupcase( band1 ), /remove_all ) 
band2 = strcompress( strupcase( band2 ), /remove_all ) 
index1 = where( filter_array eq band1 )
if ( index1 eq -1 ) then begin 
    message, 'Wrong Filter ! [ B, V, R, I ]' 
endif else begin 
    filter_index1 = index1 
    filter1 = band1
endelse
index2 = where( filter_array eq band2 )
if ( index2 eq -1 ) then begin 
    message, 'Wrong Filter ! [ B, V, R, I ]' 
endif else begin 
    filter_index2 = index2 
    filter2 = band2
endelse
;; model type 
type = strcompress( type, /remove_all )
if ( where( mtype_array eq type ) eq -1 ) then begin 
    message, 'Wrong model type ! [ as0, af0, bs0, bf0, cs0, cf0 ] '
endif 
;; number of component 
ncomp = strcompress( string( ncomp ), /remove_all )
;; rebin_factor 
if keyword_set( rebin_factor ) then begin 
    rebin_factor = float( rebin_factor )
endif else begin 
    rebin_factor = 4
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; find the galaxy in the header catalog
get_head_info, galaxy, band1, header1, exist
if ( exist eq 0 ) then begin 
    message, 'Can not find the header information for filter ' + band1  + '!!'
endif
get_head_info, galaxy, band2, header2, exist
if ( exist eq 0 ) then begin 
    message, 'Can not find the header information for filter ' + band2  + '!!'
endif
;; get the useful information
;; filter1
std_name = header1.std_name
   r20_1 = float( header1.r20 )
   r50_1 = float( header1.r50 )
   r80_1 = float( header1.r80 )
  mlim_1 = float( header1.mlim )
zpt_gsc1 = float( header1.zpt_gsc )
zpt_lan1 = float( header1.zpt_lan )
 cen_x_1 = float( header1.cen_x )
 cen_y_1 = float( header1.cen_y )
   ell_1 = float( header1.ell_e )
    pa_1 = float( header1.ell_pa )
old_expt1 = float( header1.old_expt )
 if ( zpt_lan1 gt 15.0 ) then begin 
     magzpt1 = zpt_lan1
 endif else begin 
     magzpt1 = zpt_gsc1
 endelse
;; filter2
std_name = header2.std_name
   r20_2 = float( header2.r20 )
   r50_2 = float( header2.r50 )
   r80_2 = float( header2.r80 )
  mlim_2 = float( header2.mlim )
zpt_gsc2 = float( header2.zpt_gsc )
zpt_lan2 = float( header2.zpt_lan )
 cen_x_2 = float( header2.cen_x )
 cen_y_2 = float( header2.cen_y )
   ell_2 = float( header2.ell_e )
    pa_2 = float( header2.ell_pa )
old_expt2 = float( header2.old_expt )
 if ( zpt_lan2 gt 15.0 ) then begin 
     magzpt2 = zpt_lan2
 endif else begin 
     magzpt2 = zpt_gsc2
 endelse
;; get the galaxy information 
get_ell_summary, galaxy, galaxy_info, exist
if ( exist eq 0 ) then begin 
    message, 'Something wrong with the galaxy information'
    return 
endif 
CASE band1 OF  
   'B': extinc_1 = float( galaxy_info.a_b )
   'V': extinc_1 = float( galaxy_info.a_v ) 
   'R': extinc_1 = float( galaxy_info.a_r )
   'I': extinc_1 = float( galaxy_info.a_i )
ENDCASE  
CASE band2 OF  
   'B': extinc_2 = float( galaxy_info.a_b )
   'V': extinc_2 = float( galaxy_info.a_v ) 
   'R': extinc_2 = float( galaxy_info.a_r )
   'I': extinc_2 = float( galaxy_info.a_i )
ENDCASE  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; find the old color pofile and sbp profile  
if ( filter_index1 gt filter_index2 ) then begin 
    old_prof = fileloc + galaxy + '/' + band1 + '/' + galaxy + '_' + $
        band1 + '_sbp.dat'
endif else begin 
    old_prof = fileloc + galaxy + '/' + band2 + '/' + galaxy + '_' + $
        band2 + '_sbp.dat'
endelse
if NOT file_test( old_prof ) then begin 
    print, '#######################################################'
    print, 'Can not find old sbp profile ' + old_prof + ' !!' 
    print, '#######################################################'
    old_prof_exist = 0
endif else begin 
    read_sbp, old_prof, old_sbp, old_sbo_line
    old_prof_exist = 1
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; decide the name of the note and input file 
if ( filter_index1 lt filter_index2 ) then begin 
    blue_filter = filter1 
    red_filter = filter2 
endif else begin 
    blue_filter = filter2 
    red_filter = filter1 
endelse

;; filter1
if ( filter1 eq standard_filter ) then begin 
    loc_1 = color_loc + filter1 + 'band/' 
    model_string_1 = galaxy + '_' + filter1 + '_model_' + ncomp + '0_fs0' 
    read_1 = loc_1 + model_string_1 + '.read' 
    note_1 = loc_1 + model_string_1 + '.note'
    old_ori_prof_1 = loc_1 + model_string_1 + '_ori.prof' 
    old_mod_prof_1 = loc_1 + model_string_1 + '_mod.prof' 
    model_string_loc_1 = loc_1 + model_string_1
endif else begin 
    if ( ( standard_filter ne 'V' ) and ( filter1 eq 'V' ) ) then begin 
        loc_1 = color_loc + 'V_' + filter2 + '_' + type + '/'
        model_string_1 = galaxy + '_' + filter1 + '_color_' + ncomp + '0_' + $
            type
        read_1 = loc_1 + model_string_1 + '.read'
        note_1 = loc_1 + model_string_1 + '.note'
        model_string_loc_1 = loc_1 + model_string_1
    endif else begin 
        loc_1 = color_loc + filter1 + '_' + type + '/' 
        model_string_1 = galaxy + '_' + filter1 + '_color_' + ncomp + '0_' + $
            type 
        read_1 = loc_1 + model_string_1 + '.read' 
        note_1 = loc_1 + model_string_1 + '.note'
        model_string_loc_1 = loc_1 + model_string_1
    endelse
endelse

;; filter2
if ( filter2 eq standard_filter ) then begin 
    message, 'Filter2 should not be the standard filter !' 
endif else begin 
    if ( ( standard_filter ne 'V' ) and ( filter2 eq 'V' ) ) then begin 
        loc_2 = color_loc + 'V_' + filter1 + '_' + type + '/'
        model_string_2 = galaxy + '_' + filter2 + '_color_' + ncomp + '0_' + $
            type
        read_2 = loc_2 + model_string_2 + '.read'
        note_2 = loc_2 + model_string_2 + '.note'
        model_string_loc_2 = loc_2 + model_string_2
    endif else begin 
        loc_2 = color_loc + filter2 + '_' + type + '/' 
        model_string_2 = galaxy + '_' + filter2 + '_color_' + ncomp + '0_' + $
            type 
        read_2 = loc_2 + model_string_2 + '.read' 
        note_2 = loc_2 + model_string_2 + '.note'
        model_string_loc_2 = loc_2 + model_string_2
    endelse
endelse
;; check these files, if no problem, read in ! 
;; readin files 
if NOT file_test( read_1 ) then begin 
    message, 'Can not find : ' + read_1 + '  !!' 
endif else begin 
    print, '#######################################################'
    print, 'Readin file for filter ' + filter1 + ' is : ' + read_1 
    print, '#######################################################'
    read_1_lines = file_lines( read_1 ) 
    read_1_para = strarr( read_1_lines ) 
    openr, 10, read_1 
    readf, 10, read_1_para 
    close, 10
endelse
if NOT file_test( read_2 ) then begin 
    message, 'Can not find : ' + read_2 + '  !!' 
endif else begin 
    print, '#######################################################'
    print, 'Readin file for filter ' + filter2 + ' is : ' + read_2 
    print, '#######################################################'
    read_2_lines = file_lines( read_2 ) 
    read_2_para = strarr( read_2_lines ) 
    openr, 10, read_2 
    readf, 10, read_2_para 
    close, 10
endelse
;; note file 
if NOT file_test( note_1 ) then begin 
    message, 'Can not find : ' + note_1 + '  !!' 
endif else begin 
    print, '#######################################################'
    print, 'Note file for filter ' + filter1 + ' is : ' + note_1 
    print, '#######################################################'
    note_1_lines = file_lines( note_1 ) 
    note_1_para = strarr( note_1_lines ) 
    openr, 10, note_1 
    readf, 10, note_1_para 
    close, 10
endelse
if NOT file_test( note_2 ) then begin 
    message, 'Can not find : ' + note_2 + '  !!' 
endif else begin 
    print, '#######################################################'
    print, 'Note file for filter ' + filter2 + ' is : ' + note_2 
    print, '#######################################################'
    note_2_lines = file_lines( note_2 ) 
    note_2_para = strarr( note_2_lines ) 
    openr, 10, note_2 
    readf, 10, note_2_para 
    close, 10
endelse
;; find the original images for both filter, read in and get the sky information 
ori_1 = fileloc + galaxy + '/' + filter1 + '/' + galaxy + '_' + filter1 + $
    '_cor.fit'
ori_2 = fileloc + galaxy + '/' + filter2 + '/' + galaxy + '_' + filter2 + $
    '_cor.fit'
if NOT file_test( ori_1 ) then begin 
    message, 'Can not find the original image for filter ' + filter1
endif else begin 
    img_ori_1 = mrdfits( ori_1, 0, head_ori_1 ) 
    sky_new_1 = float( sxpar( head_ori_1, 'SKY_NV' ) )
    sky_ner_1 = float( sxpar( head_ori_1, 'SKY_NE' ) ) 
    print, '#########################################################'
    print, ' GALAXY : ' + galaxy 
    print, ' SKY_NV : ' + string( sky_new_1 )
    print, ' SKY_NE : ' + string( sky_ner_1 )
    print, '#######################################################'
    if ( sky_new_1 le 1.0 ) then begin 
        print, '#######################################################'
        print, '!!!! Be careful with the sky for filter ' + filter1 
        print, '#######################################################'
    endif
endelse
if NOT file_test( ori_2 ) then begin 
    message, 'Can not find the original image for filter ' + filter2
endif else begin 
    img_ori_2 = mrdfits( ori_2, 0, head_ori_2 ) 
    sky_new_2 = float( sxpar( head_ori_2, 'SKY_NV' ) )
    sky_ner_2 = float( sxpar( head_ori_2, 'SKY_NE' ) ) 
    print, '#########################################################'
    print, ' GALAXY : ' + galaxy 
    print, ' SKY_NV : ' + string( sky_new_2 )
    print, ' SKY_NE : ' + string( sky_ner_2 )
    print, '#######################################################'
    if ( sky_new_2 le 1.0 ) then begin 
        print, '#######################################################'
        print, '!!!! Be careful with the sky for filter ' + filter2 
        print, '#######################################################'
    endif
endelse

;; find the psf images for both filter
psf_1 = fileloc + galaxy + '/' + filter1 + '/' + galaxy + '_' + filter1 + $
    '_ep.fits'
psf_2 = fileloc + galaxy + '/' + filter2 + '/' + galaxy + '_' + filter2 + $
    '_ep.fits'
if NOT file_test( psf_1 ) then begin 
    message, 'Can not find the PSF image for filter ' + filter1 + '!!' 
endif else begin 
    img_psf_1 = mrdfits( psf_1, 0 )
endelse
if NOT file_test( psf_2 ) then begin 
    message, 'Can not find the PSF image for filter ' + filter2 + '!!' 
endif else begin 
    img_psf_2 = mrdfits( psf_2, 0 )
endelse
;; find the mask images for both filter
msk_1 = fileloc + galaxy + '/' + filter1 + '/' + galaxy + '_' + filter1 + $
    '_mm.fits'
msk_2 = fileloc + galaxy + '/' + filter2 + '/' + galaxy + '_' + filter2 + $
    '_mm.fits'
if NOT file_test( msk_1 ) then begin 
    message, 'Can not find the mask image for filter ' + filter1 + '!!' 
endif
if NOT file_test( msk_2 ) then begin 
    message, 'Can not find the mask image for filter ' + filter2 + '!!' 
endif
;; find the ellipse binary file for standard filter or filter1  
standard_bin = fileloc + galaxy + '/' + filter1 + '/prof3.bin'
if NOT file_test( standard_bin ) then begin 
    standard_bin = fileloc + galaxy + '/' + filter2 + '/prof3.bin'
    if NOT file_test( standard_bin ) then begin 
        print, '#######################################################'
        print, 'Can not find the standard binary file !!' 
        print, '#######################################################'
        no_binary = 1 
    endif else begin 
        no_binary = 0 
    endelse
endif else begin 
    no_binary = 0 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; read note1 and note2 for color information of individual component
;; make sure both notes have the same number of components
temp = note_1_para[ where( strpos( note_1_para, '[  COMP_NUMBER ]' ) ne -1 ) ]
ncomp_1 = long( ( strsplit( temp, '[:] "/', /extract ) )[1] )
temp = note_2_para[ where( strpos( note_2_para, '[  COMP_NUMBER ]' ) ne -1 ) ]
ncomp_2 = long( ( strsplit( temp, '[:] "/', /extract ) )[1] )
if ( ncomp_1 ne ncomp_2 ) then begin 
    message, 'The two model have different numbers of components !' 
endif
;; create an array for the color information, one for color of every compoent, 
;; one for the average color 
  comp_color = fltarr( ncomp_1 )
galaxy_color = fltarr( 5 )
  comp_mue_1 = fltarr( ncomp_1 ) 
  comp_mue_2 = fltarr( ncomp_1 ) 
  comp_frac  = fltarr( ncomp_1 )
;; read note into summary
read_galaxy_note, note_1, model_summary=model_sum1, comp_summary=comp_sum1
read_galaxy_note, note_2, model_summary=model_sum2, comp_summary=comp_sum2
;; sort the component using radius
;; filter
radius1 = comp_sum1.r
temp_sum1 = comp_sum1
s_radius1 = bsort( radius1 ) 
temp_sum1 = comp_sum1[ s_radius1 ] 
comp_sum1 = temp_sum1
;; filter2
radius2 = comp_sum2.r
temp_sum2 = comp_sum2
s_radius2 = bsort( radius2 ) 
temp_sum2 = comp_sum2[ s_radius2 ] 
comp_sum2 = temp_sum2
;; get the color for every component; correct the Galactic extinction
for i = 0, ( ncomp - 1 ), 1 do begin 
    comp_mue_1[i] = ( comp_sum1[i].mue - extinc_1 )
    comp_mue_2[i] = ( comp_sum2[i].mue - extinc_2 )
    if ( ( comp_mue_1[i] le 0.1 ) or ( comp_mue_2[i] le 0.1 ) ) then begin 
        print, '#############################################################'
        print, 'The model does not have useful surface brightness information'
        print, '      The Mu_e is calculated from Mag !!!!!!!!             '
        print, '#############################################################'
        aa_1 = ( comp_sum1[i].mag )
        aa_2 = ( comp_sum2[i].mag )
        bb_1 = ( comp_sum1[i].r ) 
        bb_2 = ( comp_sum2[i].r ) 
        cc_1 = ( comp_sum1[i].para_a )
        cc_2 = ( comp_sum2[i].para_a )
        dd_1 = ( comp_sum1[i].para_b )
        dd_2 = ( comp_sum2[i].para_b )
        sersic_mu, aa_1, bb_1, cc_1, dd_1, ee_1, ff_1 
        sersic_mu, aa_2, bb_2, cc_2, dd_2, ee_2, ff_2 
        comp_mue_1[i] = ee_1 - extinc_1
        comp_mue_2[i] = ee_2 - extinc_2
    endif
    if ( filter_index1 gt filter_index2 ) then begin 
        comp_frac[i] = comp_sum1[i].frac 
    endif else begin 
        comp_frac[i] = comp_sum2[i].frac 
    endelse

    print, '#######################################################'
    print, 'Component ' + string( i + 1 )
    print, 'Magnitude ' + filter1 + ' : ' + string( comp_sum1[i].mag - extinc_1 ) 
    print, 'Magnitude ' + filter2 + ' : ' + string( comp_sum2[i].mag - extinc_2 ) 
    print, 'Mu_e      ' + filter1 + ' : ' + string( comp_sum1[i].mue - extinc_1 )
    print, 'Mu_e      ' + filter2 + ' : ' + string( comp_sum2[i].mue - extinc_2 )
    print, '#######################################################'
    if ( filter_index1 lt filter_index2 ) then begin 
        comp_color[i] = ( comp_sum1[i].mag - extinc_1 ) - $
            ( comp_sum2[i].mag - extinc_2 )
        print, '#######################################################'
        print, 'Color ' + filter1 + '-' + filter2 + ' : ' + $
            string( comp_color[i] )
        print, '#######################################################'
    endif
    if ( filter_index1 gt filter_index2 ) then begin 
        comp_color[i] = ( comp_sum2[i].mag - extinc_2 ) - $
            ( comp_sum1[i].mag - extinc_1 )
        print, '#######################################################'
        print, 'Color ' + filter2 + '-' + filter1 + ' : ' + $
            string( comp_color[i] )
        print, '#######################################################'
    endif
endfor
comp_mue1 = comp_mue_1
comp_mue2 = comp_mue_2

;; define the size array for the symbols of different components 
smallest = 1.5
largest  = 3.0
sep1 = ( largest - smallest ) 
sep2 = ( alog10( max( comp_frac * 100 ) ) - alog10( min( comp_frac )* 100 ) )
comp_symsize = ( sep1 / sep2 ) * ( alog10( comp_frac * 100 ) - $
    alog10( min( comp_frac ) * 100 ) ) + smallest
;; is there a sky component in GALFIT model ? 
if ( model_sum1.sky_galfit ne 0 ) then begin 
    sky_true_1 = 1 
endif else begin 
    sky_true_1 = 0 
endelse
if ( model_sum2.sky_galfit ne 0 ) then begin 
    sky_true_2 = 1 
endif else begin 
    sky_true_2 = 0 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; modified the readin file and create the sky free original image and the 
;; model image
;; model image filter1
;; model withtout psf 
model_nosky_1 = model_string_loc_1 + '_nosky_mod.fit'
if NOT file_test( model_nosky_1 ) then begin 
    model_gen_nosky, read_1, /nopsf, /run 
    if NOT file_test( model_nosky_1 ) then begin 
        message, 'Can not find : ' + model_nosky_1 + '!!' 
    endif else begin 
        img_mod_nosky_1 = mrdfits( model_nosky_1, 0 ) 
    endelse
endif else begin 
    img_mod_nosky_1 = mrdfits( model_nosky_1, 0 ) 
endelse

;; model with psf
model_nosky_1p = model_string_loc_1 + '_nosky_mof.fit'
if NOT file_test( model_nosky_1p ) then begin 
    model_gen_nosky, read_1, /run 
    if NOT file_test( model_nosky_1p ) then begin 
        message, 'Can not find : ' + model_nosky_1p + '!!' 
    endif else begin
        img_mod_nosky_1p = mrdfits( model_nosky_1p, 0 ) 
    endelse
endif else begin 
    img_mod_nosky_1p = mrdfits( model_nosky_1p, 0 ) 
endelse
;; model image filter2
;; model without psf 
model_nosky_2 = model_string_loc_2 + '_nosky_mod.fit'
if NOT file_test( model_nosky_2 ) then begin 
    model_gen_nosky, read_2, /nopsf, /run 
    if NOT file_test( model_nosky_2 ) then begin 
        message, 'Can not find : ' + model_nosky_2 + '!!' 
    endif else begin 
        img_mod_nosky_2 = mrdfits( model_nosky_2, 0 ) 
    endelse
endif else begin 
    img_mod_nosky_2 = mrdfits( model_nosky_2, 0 ) 
endelse
;; model with psf
model_nosky_2p = model_string_loc_2 + '_nosky_mof.fit'
if NOT file_test( model_nosky_2p ) then begin 
    model_gen_nosky, read_2, /run 
    if NOT file_test( model_nosky_2p ) then begin 
        message, 'Can not find : ' + model_nosky_2p + '!!' 
    endif else begin
        img_mod_nosky_2p = mrdfits( model_nosky_2p, 0 ) 
    endelse
endif else begin 
    img_mod_nosky_2p = mrdfits( model_nosky_2p, 0 ) 
endelse
;; calculate the average color 
flux_1 = total( img_mod_nosky_1p ) 
flux_2 = total( img_mod_nosky_2p )
mag_total_1 = -2.5 * alog10( flux_1 / old_expt1 ) + magzpt1 - extinc_1
mag_total_2 = -2.5 * alog10( flux_2 / old_expt2 ) + magzpt2 - extinc_2
if ( filter_index1 lt filter_index2 ) then begin 
    color_average = ( mag_total_1 - mag_total_2 )
endif else begin 
    color_average = ( mag_total_2 - mag_total_1 )
endelse
galaxy_color[0] = color_average
print, '#########################################################'
print, 'The Average Color : ' + string( color_average ) 
print, '#########################################################'
;; original image filter1
ori_nosky_1 = model_string_loc_1 + '_nosky_ori.fit'
if NOT file_test( ori_nosky_1 ) then begin 
    model_sky_subtract, read_1 
endif
if NOT file_test( ori_nosky_1 ) then begin 
    message, 'Can not find : ' + ori_nosky_1 + '!!' 
endif else begin 
    img_ori_nosky_1 = mrdfits( ori_nosky_1, 0, head_ori_nosky_1 )
endelse
;; original image filter2
ori_nosky_2 = model_string_loc_2 + '_nosky_ori.fit'
if NOT file_test( ori_nosky_2 ) then begin 
    model_sky_subtract, read_2 
endif
if NOT file_test( ori_nosky_2 ) then begin 
    message, 'Can not find : ' + ori_nosky_2 + '!!' 
endif else begin 
    img_ori_nosky_2 = mrdfits( ori_nosky_2, 0, head_ori_nosky_2 )
endelse
;; convole the original image with the PSF image of the other filter
;; filter1
print, '#######################################################'
img_ori_conv_1 = convol( img_ori_nosky_1, img_psf_2, /edge_zero, /normalize )
ori_conv_1 = model_string_loc_1 + '_ori_conv.fit'
mwrfits, img_ori_conv_1, ori_conv_1, head_ori_nosky_1, /create 
;; filter2
img_ori_conv_2 = convol( img_ori_nosky_2, img_psf_2, /edge_zero, /normalize )
ori_conv_2 = model_string_loc_2 + '_ori_conv.fit'
mwrfits, img_ori_conv_2, ori_conv_2, head_ori_nosky_2, /create 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; run ellipse using the old binary file from standard filter on the original 
;; and model images for both filters
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 1. ellipse run on original image of filter 1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the input files
 input_img = ori_nosky_1
  mask_img = msk_1
    pl_img = model_string_loc_1 + '_nosky_ori.pl'
output_bin = model_string_loc_1 + '_nosky_ori.bin'
output_asc = model_string_loc_1 + '_nosky_ori.prof'
output_head = 'head.dat'

if NOT file_test( output_asc ) then begin 

    ellip_file = model_string_loc_1 + '_nosky_ori.ellip'
    ;; convert the mask file into pl format
    if NOT file_test( pl_img ) then begin 
        para_line = file_lines( imcopy_par ) 
        para = strarr( para_line ) 
        openr, 10, imcopy_par 
        readf, 10, para 
        close, 10 
        para[0] = 'imcopy.input = "' + mask_img + '"'
        para[1] = 'imcopy.output = "' + pl_img + '"'
        openw, 10, 'imcopy.par'
        for j = 0, ( para_line - 1 ), 1 do begin
            printf, 10, para[j]
        endfor
        close, 10
        spawn, ximage + ' imcopy @imcopy.par'     
        if NOT file_test( pl_img ) then begin 
            message, 'Somthing Wrong with the .pl Mask File for filter ' + filter1
        endif
    endif
    ;; get the necessary information
    x0 = strcompress( string( cen_x_1 ), /remove_all )
    y0 = strcompress( string( cen_y_1 ), /remove_all )
    ellip0 = strcompress( string( ell_1 ), /remove_all )
    pa0 = strcompress( string( pa_1 ), /remove_all )
    sma0 = strcompress( string( sma0 ), /remove_all )
    minsma = strcompress( string( minsma ), /remove_all )
    maxsma = strcompress( string( maxsma ), /remove_all )
    step = strcompress( string( step ), /remove_all ) 
    mag0 = strcompress( string( magzpt1 ), /remove_all ) 
    
    ;; edit the ellipse parameter file
    num_ellpara = file_lines( ell_par ) 
    ell_para = strarr( num_ellpara )
    openr, 10, ell_par
    readf, 10, ell_para
    close, 10
    ell_para[0] = 'ellipse.input = "' + input_img + '"'
    ell_para[1] = 'ellipse.output = "' + output_bin + '"'
    if NOT keyword_set( no_input_binary ) then begin 
        ell_para[3] = 'ellipse.inellip = "' + standard_bin + '"'
    endif else begin 
        ell_para[3] = 'ellipse.inellip = ""'
    endelse
    ell_para[27] = 'controlpar.hellip = yes '
    ell_para[28] = 'controlpar.hpa = yes '
    ell_para[34] = 'geompar.x0 = ' + x0
    ell_para[35] = 'geompar.y0 = ' + y0
    ell_para[36] = 'geompar.ellip0 = ' + ellip0
    ell_para[37] = 'geompar.pa0 = ' + pa0
    ell_para[38] = 'geompar.sma0 = ' + sma0
    ell_para[39] = 'geompar.minsma = ' + minsma
    ell_para[40] = 'geompar.maxsma = ' + maxsma
    ell_para[41] = 'geompar.step = ' + step
    ell_para[48] = 'magpar.mag0 = ' + mag0
    openw, 10, ellip_file 
    for k = 0, num_ellpara-1, 1 do begin
        printf, 10, ell_para[k]
    endfor
    close, 10
    ;; remove the previous ellipse results
    if file_test( output_bin ) then begin 
        spawn, 'rm  ' + output_bin
    endif
    if file_test( output_asc ) then begin 
        spawn, 'rm  ' + output_asc
    endif
    if file_test( output_head ) then begin 
        spawn, 'rm  ' + output_head
    endif
    ;; Run Ellipse 
    spawn, isophote + ' ellipse @' +  ellip_file
    ;; if no binary file, use the outpu binary as inellip
    ;; Test if the ellipse fitting goes well
    if NOT file_test( output_bin ) then begin 
        Print, 'Somthing Wrong with the Ellipse Fitting using Old Binary File'
        ;; if the ellipse run using old binary file fails, start another run 
        ;; using the opposite choice 
        ;; edit the ellipse parameter file
        num_ellpara = file_lines( ell_par ) 
        ell_para = strarr( num_ellpara )
        openr, 10, ell_par
        readf, 10, ell_para
        close, 10
        if NOT keyword_set( no_input_binary ) then begin 
            ell_para[0] = 'ellipse.input = "' + input_img + '"'
            ell_para[1] = 'ellipse.output = "' + output_bin + '"'
            ell_para[3] = 'ellipse.inellip = ""'
            ell_para[27] = 'controlpar.hellip = yes '
            ell_para[28] = 'controlpar.hpa = yes '
            ell_para[34] = 'geompar.x0 = ' + x0
            ell_para[35] = 'geompar.y0 = ' + y0
            ell_para[36] = 'geompar.ellip0 = ' + ellip0
            ell_para[37] = 'geompar.pa0 = ' + pa0
            ell_para[38] = 'geompar.sma0 = ' + sma0
            ell_para[39] = 'geompar.minsma = ' + minsma
            ell_para[40] = 'geompar.maxsma = ' + maxsma
            ell_para[41] = 'geompar.step = ' + step
            ell_para[48] = 'magpar.mag0 = ' + mag0
        endif else begin 
            step = ( float( step ) + 0.02 )
            step = strcompress( string( step ), /remove_all )
            maxsma = ( float( maxsma ) * 0.9 ) 
            maxsma = strcompress( string( maxsma ), /remove_all )
            ell_para[0] = 'ellipse.input = "' + input_img + '"'
            ell_para[1] = 'ellipse.output = "' + output_bin + '"'
            ell_para[3] = 'ellipse.inellip = ""'
            ell_para[27] = 'controlpar.hellip = yes '
            ell_para[28] = 'controlpar.hpa = yes '
            ell_para[34] = 'geompar.x0 = ' + x0
            ell_para[35] = 'geompar.y0 = ' + y0
            ell_para[36] = 'geompar.ellip0 = ' + ellip0
            ell_para[37] = 'geompar.pa0 = ' + pa0
            ell_para[38] = 'geompar.sma0 = ' + sma0
            ell_para[39] = 'geompar.minsma = ' + minsma
            ell_para[40] = 'geompar.maxsma = ' + maxsma
            ell_para[41] = 'geompar.step = ' + step
            ell_para[48] = 'magpar.mag0 = ' + mag0
        endelse
        openw, 10, ellip_file 
        for k = 0, num_ellpara-1, 1 do begin
            printf, 10, ell_para[k]
        endfor
        close, 10
        ;; remove the previous ellipse results
        if file_test( output_bin ) then begin 
            spawn, 'rm  ' + output_bin
        endif
        if file_test( output_asc ) then begin 
            spawn, 'rm  ' + output_asc
        endif
        if file_test( output_head ) then begin 
            spawn, 'rm  ' + output_head
        endif
        ;; Run Ellipse 
        spawn, isophote + ' ellipse @' +  ellip_file
        ;; Test if the ellipse run goes well 
        if NOT file_test( output_bin ) then begin 
            message, 'Somthing wrong with the Ellipse Fitting, Check Again !!' 
        endif  
    endif 
    ;; change the binary table into ascii table
    print, 'Output Binary File: ' + output_bin
    num_tpara = file_lines( tdump_par )
    tdump_para = strarr( num_tpara )
    openr, 10, tdump_par
    readf, 10, tdump_para
    close, 10
    tdump_para[0] = 'tdump.table = "' +  output_bin + '"'
    tdump_para[2] = 'tdump.pfile = "' +  output_head + '"'
    tdump_para[3] = 'tdump.datafile = "' + output_asc + '"'
    openw, 10, 'tdump.par'
    for m = 0, num_tpara-1, 1 do begin
        printf, 10, tdump_para[m]
    endfor
    close, 10
    spawn, ttools + ' tdump @tdump.par'
    spawn, 'rm ' + output_head
    if file_test( output_asc ) then begin 
        ;spawn, 'sed -in-place -e "s/INDEF/-9999.9/g" ' + output_asc 
        spawn, 'sed -i "s/INDEF/-9999.9/g" ' + output_asc 
    endif

endif
;; readin the surface brightness profile 
old_ori_prof_1 = output_asc
print, old_ori_prof_1
read_profile, old_ori_prof_1, old_ori_sbp_1, old_ori_line_1 
 old_ori_sma_1 = old_ori_sbp_1.sma 
old_ori_rsma_1 = ( old_ori_sma_1 * pix )^0.25
old_ori_smag_1 = -2.5 * alog10( old_ori_sbp_1.intens / $
    ( pix_area * old_expt1 ) ) + magzpt1 - extinc_1
 old_ori_sig_1 = sqrt( ( old_ori_sbp_1.intens_err )^2.0 + ( sky_ner_1 )^2.0 ) 
old_ori_serr_1 = 2.5 * alog10( 1.0 + ( old_ori_sig_1 / $
    ( old_ori_sbp_1.intens ) ) )
;; change the standard_bin to the output one 
if keyword_set( no_input_binary ) then begin 
    standard_bin = output_bin
endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 2. ellipse run on original image of filter 2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the input files
 input_img = ori_nosky_2
  mask_img = msk_2
    pl_img = model_string_loc_2 + '_nosky_ori.pl'
output_bin = model_string_loc_2 + '_nosky_ori.bin'
output_asc = model_string_loc_2 + '_nosky_ori.prof'
output_head = 'head.dat'
if NOT file_test( output_asc ) then begin 

    ellip_file = model_string_loc_2 + '_nosky_ori.ellip'
    ;; convert the mask file into pl format
    if NOT file_test( pl_img ) then begin 
        para_line = file_lines( imcopy_par ) 
        para = strarr( para_line ) 
        openr, 10, imcopy_par 
        readf, 10, para 
        close, 10 
        para[0] = 'imcopy.input = "' + mask_img + '"'
        para[1] = 'imcopy.output = "' + pl_img + '"'
        openw, 10, 'imcopy.par'
        for j = 0, ( para_line - 1 ), 1 do begin
            printf, 10, para[j]
        endfor
        close, 10
        spawn, ximage + ' imcopy @imcopy.par'     
        if NOT file_test( pl_img ) then begin 
            message, 'Somthing Wrong with the .pl Mask File for filter ' + filter1
        endif
    endif
    ;; get the necessary information
    x0 = strcompress( string( cen_x_2 ), /remove_all )
    y0 = strcompress( string( cen_y_2 ), /remove_all )
    ellip0 = strcompress( string( ell_2 ), /remove_all )
    pa0 = strcompress( string( pa_2 ), /remove_all )
    sma0 = strcompress( string( 30.0 ), /remove_all )
    minsma = strcompress( string( 0.0 ), /remove_all )
    maxsma = strcompress( string( 1080.0 ), /remove_all )
    step = strcompress( string( 0.05 ), /remove_all ) 
    mag0 = strcompress( string( magzpt2 ), /remove_all ) 
    ;; edit the ellipse parameter file
    num_ellpara = file_lines( ell_par ) 
    ell_para = strarr( num_ellpara )
    openr, 10, ell_par
    readf, 10, ell_para
    close, 10
    ell_para[0] = 'ellipse.input = "' + input_img + '"'
    ell_para[1] = 'ellipse.output = "' + output_bin + '"'
    ell_para[3] = 'ellipse.inellip = "' + standard_bin + '"'
    ell_para[27] = 'controlpar.hellip = yes '
    ell_para[28] = 'controlpar.hpa = yes '
    ell_para[34] = 'geompar.x0 = ' + x0
    ell_para[35] = 'geompar.y0 = ' + y0
    ell_para[36] = 'geompar.ellip0 = ' + ellip0
    ell_para[37] = 'geompar.pa0 = ' + pa0
    ell_para[38] = 'geompar.sma0 = ' + sma0
    ell_para[39] = 'geompar.minsma = ' + minsma
    ell_para[40] = 'geompar.maxsma = ' + maxsma
    ell_para[41] = 'geompar.step = ' + step
    ell_para[48] = 'magpar.mag0 = ' + mag0
    openw, 10, ellip_file 
    for k = 0, num_ellpara-1, 1 do begin
        printf, 10, ell_para[k]
    endfor
    close, 10
    ;; remove the previous ellipse results
    if file_test( output_bin ) then begin 
        spawn, 'rm  ' + output_bin
    endif
    if file_test( output_asc ) then begin 
        spawn, 'rm  ' + output_asc
    endif
    if file_test( output_head ) then begin 
        spawn, 'rm  ' + output_head
    endif
    ;; Run Ellipse 
    spawn, isophote + ' ellipse @' +  ellip_file
    ;; Test if the ellipse fitting goes well
    if NOT file_test( output_bin ) then begin 
        message, 'Somthing Wrong with the Ellipse Fitting, Check Again !'
    endif else begin 
    ;; change the binary table into ascii table
        print, 'Output Binary File: ' + output_bin
        num_tpara = file_lines( tdump_par )
        tdump_para = strarr( num_tpara )
        openr, 10, tdump_par
        readf, 10, tdump_para
        close, 10
        tdump_para[0] = 'tdump.table = "' +  output_bin + '"'
        tdump_para[2] = 'tdump.pfile = "' +  output_head + '"'
        tdump_para[3] = 'tdump.datafile = "' + output_asc + '"'
        openw, 10, 'tdump.par'
        for m = 0, num_tpara-1, 1 do begin
            printf, 10, tdump_para[m]
        endfor
        close, 10
        spawn, ttools + ' tdump @tdump.par'
        spawn, 'rm ' + output_head
        if file_test( output_asc ) then begin 
            ;spawn, 'sed -in-place -e "s/INDEF/-9999.9/g" ' + output_asc 
            spawn, 'sed -i "s/INDEF/-9999.9/g" ' + output_asc 
        endif
    endelse

endif
;; readin the surface brightness profile 
old_ori_prof_2 = output_asc
read_profile, old_ori_prof_2, old_ori_sbp_2, old_ori_line_2 
 old_ori_sma_2 = old_ori_sbp_2.sma 
old_ori_rsma_2 = ( old_ori_sma_2 * pix )^0.25
old_ori_smag_2 = -2.5 * alog10( old_ori_sbp_2.intens / $
    ( pix_area * old_expt2 ) ) + magzpt2 - extinc_2
 old_ori_sig_2 = sqrt( ( old_ori_sbp_2.intens_err )^2.0 + ( sky_ner_2 )^2.0 ) 
old_ori_serr_2 = 2.5 * alog10( 1.0 + ( old_ori_sig_2 / $
    ( old_ori_sbp_2.intens ) ) )
;; get the color profile and color error profile for old_original images
if ( filter_index1 lt filter_index2 ) then begin 
    old_ori_cprof = old_ori_smag_1 - old_ori_smag_2 
endif else begin 
    old_ori_cprof = old_ori_smag_2 - old_ori_smag_1 
endelse
old_ori_cerr = sqrt( ( old_ori_serr_1 )^2.0 + ( old_ori_serr_2 )^2.0 )
old_ori_color_median = median( old_ori_cprof )
galaxy_color[1] = old_ori_color_median
print, '############################################################'
print, 'Median old original Color : ' + string( old_ori_color_median )
print, '############################################################'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 3. ellipse run on convolved original image of filter 1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the input files
 input_img = ori_conv_1
  mask_img = msk_1
    pl_img = model_string_loc_1 + '_ori_conv.pl'
output_bin = model_string_loc_1 + '_ori_conv.bin'
output_asc = model_string_loc_1 + '_ori_conv.prof'
output_head = 'head.dat'
if NOT file_test( output_asc ) then begin 

    ellip_file = model_string_loc_1 + '_ori_conv.ellip'
    ;; convert the mask file into pl format
    if NOT file_test( pl_img ) then begin 
        para_line = file_lines( imcopy_par ) 
        para = strarr( para_line ) 
        openr, 10, imcopy_par 
        readf, 10, para 
        close, 10 
        para[0] = 'imcopy.input = "' + mask_img + '"'
        para[1] = 'imcopy.output = "' + pl_img + '"'
        openw, 10, 'imcopy.par'
        for j = 0, ( para_line - 1 ), 1 do begin
            printf, 10, para[j]
        endfor
        close, 10
        spawn, ximage + ' imcopy @imcopy.par'     
        if NOT file_test( pl_img ) then begin 
            message, 'Somthing Wrong with the .pl Mask File for filter ' + $
                filter1
        endif
    endif
    ;; get the necessary information
    x0 = strcompress( string( cen_x_1 ), /remove_all )
    y0 = strcompress( string( cen_y_1 ), /remove_all )
    ellip0 = strcompress( string( ell_1 ), /remove_all )
    pa0 = strcompress( string( pa_1 ), /remove_all )
    sma0 = strcompress( string( 40.0 ), /remove_all )
    minsma = strcompress( string( 0.0 ), /remove_all )
    maxsma = strcompress( string( 1080.0 ), /remove_all )
    step = strcompress( string( 0.05 ), /remove_all ) 
    mag0 = strcompress( string( magzpt1 ), /remove_all ) 
    ;; edit the ellipse parameter file
    num_ellpara = file_lines( ell_par ) 
    ell_para = strarr( num_ellpara )
    openr, 10, ell_par
    readf, 10, ell_para
    close, 10
    ell_para[0] = 'ellipse.input = "' + input_img + '"'
    ell_para[1] = 'ellipse.output = "' + output_bin + '"'
    ell_para[3] = 'ellipse.inellip = "' + standard_bin + '"'
    ell_para[27] = 'controlpar.hellip = yes '
    ell_para[28] = 'controlpar.hpa = yes '
    ell_para[34] = 'geompar.x0 = ' + x0
    ell_para[35] = 'geompar.y0 = ' + y0
    ell_para[36] = 'geompar.ellip0 = ' + ellip0
    ell_para[37] = 'geompar.pa0 = ' + pa0
    ell_para[38] = 'geompar.sma0 = ' + sma0
    ell_para[39] = 'geompar.minsma = ' + minsma
    ell_para[40] = 'geompar.maxsma = ' + maxsma
    ell_para[41] = 'geompar.step = ' + step
    ell_para[48] = 'magpar.mag0 = ' + mag0
    openw, 10, ellip_file 
    for k = 0, num_ellpara-1, 1 do begin
        printf, 10, ell_para[k]
    endfor
    close, 10
    ;; remove the previous ellipse results
    if file_test( output_bin ) then begin 
        spawn, 'rm  ' + output_bin
    endif
    if file_test( output_asc ) then begin 
        spawn, 'rm  ' + output_asc
    endif
    if file_test( output_head ) then begin 
        spawn, 'rm  ' + output_head
    endif
    ;; Run Ellipse 
    spawn, isophote + ' ellipse @' +  ellip_file
    ;; Test if the ellipse fitting goes well
    if NOT file_test( output_bin ) then begin 
        message, 'Somthing wrong with the Ellipse Fitting, Check Again !!' 
    endif 
    ;; change the binary table into ascii table
    print, 'Output Binary File: ' + output_bin
    num_tpara = file_lines( tdump_par )
    tdump_para = strarr( num_tpara )
    openr, 10, tdump_par
    readf, 10, tdump_para
    close, 10
    tdump_para[0] = 'tdump.table = "' +  output_bin + '"'
    tdump_para[2] = 'tdump.pfile = "' +  output_head + '"'
    tdump_para[3] = 'tdump.datafile = "' + output_asc + '"'
    openw, 10, 'tdump.par'
    for m = 0, num_tpara-1, 1 do begin
        printf, 10, tdump_para[m]
    endfor
    close, 10
    spawn, ttools + ' tdump @tdump.par'
    spawn, 'rm ' + output_head
    if file_test( output_asc ) then begin 
        ;spawn, 'sed -in-place -e "s/INDEF/-9999.9/g" ' + output_asc 
        spawn, 'sed -i "s/INDEF/-9999.9/g" ' + output_asc 
    endif

endif
;; readin the surface brightness profile 
ori_prof_1 = output_asc
read_profile, ori_prof_1, ori_sbp_1, ori_line_1 
 ori_sma_1 = ori_sbp_1.sma 
ori_rsma_1 = ( ori_sma_1 * pix )^0.25
ori_smag_1 = -2.5 * alog10( ori_sbp_1.intens / ( pix_area * old_expt1 ) ) + $
    magzpt1 - extinc_1
 ori_sig_1 = sqrt( ( ori_sbp_1.intens_err )^2.0 + ( sky_ner_1 )^2.0 ) 
ori_serr_1 = 2.5 * alog10( 1.0 + ( ori_sig_1 / ( ori_sbp_1.intens ) ) )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 4. ellipse run on convolved original image of filter 2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the input files
 input_img = ori_conv_2
  mask_img = msk_2
    pl_img = model_string_loc_2 + '_ori_conv.pl'
output_bin = model_string_loc_2 + '_ori_conv.bin'
output_asc = model_string_loc_2 + '_ori_conv.prof'
output_head = 'head.dat'
if NOT file_test( output_asc ) then begin 

    ellip_file = model_string_loc_2 + '_ori_conv.ellip'
    ;; convert the mask file into pl format
    if NOT file_test( pl_img ) then begin 
        para_line = file_lines( imcopy_par ) 
        para = strarr( para_line ) 
        openr, 10, imcopy_par 
        readf, 10, para 
        close, 10 
        para[0] = 'imcopy.input = "' + mask_img + '"'
        para[1] = 'imcopy.output = "' + pl_img + '"'
        openw, 10, 'imcopy.par'
        for j = 0, ( para_line - 1 ), 1 do begin
            printf, 10, para[j]
        endfor
        close, 10
        spawn, ximage + ' imcopy @imcopy.par'     
        if NOT file_test( pl_img ) then begin 
            message, 'Somthing Wrong with the .pl Mask File for filter ' + filter1
        endif
    endif
    ;; get the necessary information
    x0 = strcompress( string( cen_x_2 ), /remove_all )
    y0 = strcompress( string( cen_y_2 ), /remove_all )
    ellip0 = strcompress( string( ell_2 ), /remove_all )
    pa0 = strcompress( string( pa_2 ), /remove_all )
    sma0 = strcompress( string( 30.0 ), /remove_all )
    minsma = strcompress( string( 0.0 ), /remove_all )
    maxsma = strcompress( string( 1080.0 ), /remove_all )
    step = strcompress( string( 0.05 ), /remove_all ) 
    mag0 = strcompress( string( magzpt2 ), /remove_all ) 
    ;; edit the ellipse parameter file
    num_ellpara = file_lines( ell_par ) 
    ell_para = strarr( num_ellpara )
    openr, 10, ell_par
    readf, 10, ell_para
    close, 10
    ell_para[0] = 'ellipse.input = "' + input_img + '"'
    ell_para[1] = 'ellipse.output = "' + output_bin + '"'
    ell_para[3] = 'ellipse.inellip = "' + standard_bin + '"'
    ell_para[27] = 'controlpar.hellip = yes '
    ell_para[28] = 'controlpar.hpa = yes '
    ell_para[34] = 'geompar.x0 = ' + x0
    ell_para[35] = 'geompar.y0 = ' + y0
    ell_para[36] = 'geompar.ellip0 = ' + ellip0
    ell_para[37] = 'geompar.pa0 = ' + pa0
    ell_para[38] = 'geompar.sma0 = ' + sma0
    ell_para[39] = 'geompar.minsma = ' + minsma
    ell_para[40] = 'geompar.maxsma = ' + maxsma
    ell_para[41] = 'geompar.step = ' + step
    ell_para[48] = 'magpar.mag0 = ' + mag0
    openw, 10, ellip_file 
    for k = 0, num_ellpara-1, 1 do begin
        printf, 10, ell_para[k]
    endfor
    close, 10
    ;; remove the previous ellipse results
    if file_test( output_bin ) then begin 
        spawn, 'rm  ' + output_bin
    endif
    if file_test( output_asc ) then begin 
        spawn, 'rm  ' + output_asc
    endif
    if file_test( output_head ) then begin 
        spawn, 'rm  ' + output_head
    endif
    ;; Run Ellipse 
    spawn, isophote + ' ellipse @' +  ellip_file
    ;; Test if the ellipse fitting goes well
    if NOT file_test( output_bin ) then begin 
        message, 'Somthing Wrong with the Ellipse Fitting, Check Again !'
    endif else begin 
    ;; change the binary table into ascii table
        print, 'Output Binary File: ' + output_bin
        num_tpara = file_lines( tdump_par )
        tdump_para = strarr( num_tpara )
        openr, 10, tdump_par
        readf, 10, tdump_para
        close, 10
        tdump_para[0] = 'tdump.table = "' +  output_bin + '"'
        tdump_para[2] = 'tdump.pfile = "' +  output_head + '"'
        tdump_para[3] = 'tdump.datafile = "' + output_asc + '"'
        openw, 10, 'tdump.par'
        for m = 0, num_tpara-1, 1 do begin
            printf, 10, tdump_para[m]
        endfor
        close, 10
        spawn, ttools + ' tdump @tdump.par'
        spawn, 'rm ' + output_head
        if file_test( output_asc ) then begin 
            ;spawn, 'sed -in-place -e "s/INDEF/-9999.9/g" ' + output_asc 
            spawn, 'sed -i "s/INDEF/-9999.9/g" ' + output_asc 
        endif
    endelse

endif
;; readin the surface brightness profile 
ori_prof_2 = output_asc
read_profile, ori_prof_2, ori_sbp_2, ori_line_2 
 ori_sma_2 = ori_sbp_2.sma 
ori_rsma_2 = ( ori_sma_2 * pix )^0.25
ori_smag_2 = -2.5 * alog10( ori_sbp_2.intens / ( pix_area * old_expt2 ) ) + $
    magzpt2 - extinc_2
 ori_sig_2 = sqrt( ( ori_sbp_2.intens_err )^2.0 + ( sky_ner_2 )^2.0 ) 
ori_serr_2 = 2.5 * alog10( 1.0 + ( ori_sig_2 / ( ori_sbp_2.intens ) ) )
;; get the color profile and color error profile for original images
if ( filter_index1 lt filter_index2 ) then begin 
    ori_cprof = ori_smag_1 - ori_smag_2 
endif else begin 
    ori_cprof = ori_smag_2 - ori_smag_1 
endelse
ori_cerr = sqrt( ( ori_serr_1 )^2.0 + ( ori_serr_2 )^2.0 )
ori_color_median = median( ori_cprof )
galaxy_color[2] = ori_color_median
print, '############################################################'
print, 'Median Original Color (Convolved): ' + string( ori_color_median )
print, '############################################################'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 5. ellipse run on model image of filter 1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the input files
 input_img = model_nosky_1p
output_bin = model_string_loc_1 + '_nosky_mof.bin'
output_asc = model_string_loc_1 + '_nosky_mof.prof'
ellip_file = model_string_loc_1 + '_nosky_mof.ellip'
output_head = 'head.dat'

if NOT file_test( output_asc ) then begin 

    ;; get the necessary information
    x0 = strcompress( string( cen_x_1 ), /remove_all )
    y0 = strcompress( string( cen_y_1 ), /remove_all )
    ellip0 = strcompress( string( ell_1 ), /remove_all )
    pa0 = strcompress( string( pa_1 ), /remove_all )
    sma0 = strcompress( string( 30.0 ), /remove_all )
    minsma = strcompress( string( 0.0 ), /remove_all )
    maxsma = strcompress( string( 1080.0 ), /remove_all )
    step = strcompress( string( 0.05 ), /remove_all ) 
    mag0 = strcompress( string( magzpt1 ), /remove_all ) 
    ;; edit the ellipse parameter file
    num_ellpara = file_lines( ell_par ) 
    ell_para = strarr( num_ellpara )
    openr, 10, ell_par
    readf, 10, ell_para
    close, 10
    ell_para[0] = 'ellipse.input = "' + input_img + '"'
    ell_para[1] = 'ellipse.output = "' + output_bin + '"'
    ell_para[3] = 'ellipse.inellip = "' + standard_bin + '"'
    ell_para[27] = 'controlpar.hellip = yes '
    ell_para[28] = 'controlpar.hpa = yes '
    ell_para[34] = 'geompar.x0 = ' + x0
    ell_para[35] = 'geompar.y0 = ' + y0
    ell_para[36] = 'geompar.ellip0 = ' + ellip0
    ell_para[37] = 'geompar.pa0 = ' + pa0
    ell_para[38] = 'geompar.sma0 = ' + sma0
    ell_para[39] = 'geompar.minsma = ' + minsma
    ell_para[40] = 'geompar.maxsma = ' + maxsma
    ell_para[41] = 'geompar.step = ' + step
    ell_para[48] = 'magpar.mag0 = ' + mag0
    openw, 10, ellip_file 
    for k = 0, num_ellpara-1, 1 do begin
        printf, 10, ell_para[k]
    endfor
    close, 10
    ;; remove the previous ellipse results
    if file_test( output_bin ) then begin 
        spawn, 'rm  ' + output_bin
    endif
    if file_test( output_asc ) then begin 
        spawn, 'rm  ' + output_asc
    endif
    if file_test( output_head ) then begin 
        spawn, 'rm  ' + output_head
    endif
    ;; Run Ellipse 
    spawn, isophote + ' ellipse @' +  ellip_file
    ;; Test if the ellipse fitting goes well
    if NOT file_test( output_bin ) then begin 
        message, 'Somthing Wrong with the Ellipse Fitting, Check Again !'
    endif else begin 
    ;; change the binary table into ascii table
        print, 'Output Binary File: ' + output_bin
        num_tpara = file_lines( tdump_par )
        tdump_para = strarr( num_tpara )
        openr, 10, tdump_par
        readf, 10, tdump_para
        close, 10
        tdump_para[0] = 'tdump.table = "' +  output_bin + '"'
        tdump_para[2] = 'tdump.pfile = "' +  output_head + '"'
        tdump_para[3] = 'tdump.datafile = "' + output_asc + '"'
        openw, 10, 'tdump.par'
        for m = 0, num_tpara-1, 1 do begin
            printf, 10, tdump_para[m]
        endfor
        close, 10
        spawn, ttools + ' tdump @tdump.par'
        spawn, 'rm ' + output_head
        if file_test( output_asc ) then begin 
            ;spawn, 'sed -in-place -e "s/INDEF/-9999.9/g" ' + output_asc 
            spawn, 'sed -i "s/INDEF/-9999.9/g" ' + output_asc 
        endif
    endelse

endif
;; readin the surface brightness profile 
old_mod_prof_1 = output_asc
read_profile, old_mod_prof_1, old_mod_sbp_1, old_mod_line_1
 old_mod_sma_1 = old_mod_sbp_1.sma 
old_mod_rsma_1 = ( old_mod_sma_1 * pix )^0.25
old_mod_smag_1 = -2.5 * alog10( old_mod_sbp_1.intens / $
    ( pix_area * old_expt1 ) ) + magzpt1 - extinc_1
 old_mod_sig_1 = sqrt( ( old_mod_sbp_1.intens_err )^2.0 + ( sky_ner_1 )^2.0 ) 
old_mod_serr_1 = 2.5 * alog10( 1.0 + ( old_mod_sig_1 / ( old_mod_sbp_1.intens ) ) )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 6. ellipse run on model image of filter 2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the input files
 input_img = model_nosky_2p
output_bin = model_string_loc_2 + '_nosky_mof.bin'
output_asc = model_string_loc_2 + '_nosky_mof.prof'
ellip_file = model_string_loc_2 + '_nosky_mof.ellip'
output_head = 'head.dat'

if NOT file_test( output_asc ) then begin 

    ;; get the necessary information
    x0 = strcompress( string( cen_x_2 ), /remove_all )
    y0 = strcompress( string( cen_y_2 ), /remove_all )
    ellip0 = strcompress( string( ell_2 ), /remove_all )
    pa0 = strcompress( string( pa_2 ), /remove_all )
    sma0 = strcompress( string( 30.0 ), /remove_all )
    minsma = strcompress( string( 0.0 ), /remove_all )
    maxsma = strcompress( string( 1080.0 ), /remove_all )
    step = strcompress( string( 0.05 ), /remove_all ) 
    mag0 = strcompress( string( magzpt2 ), /remove_all ) 
    ;; edit the ellipse parameter file
    num_ellpara = file_lines( ell_par ) 
    ell_para = strarr( num_ellpara )
    openr, 10, ell_par
    readf, 10, ell_para
    close, 10
    ell_para[0] = 'ellipse.input = "' + input_img + '"'
    ell_para[1] = 'ellipse.output = "' + output_bin + '"'
    ell_para[3] = 'ellipse.inellip = "' + standard_bin + '"'
    ell_para[27] = 'controlpar.hellip = yes '
    ell_para[28] = 'controlpar.hpa = yes '
    ell_para[34] = 'geompar.x0 = ' + x0
    ell_para[35] = 'geompar.y0 = ' + y0
    ell_para[36] = 'geompar.ellip0 = ' + ellip0
    ell_para[37] = 'geompar.pa0 = ' + pa0
    ell_para[38] = 'geompar.sma0 = ' + sma0
    ell_para[39] = 'geompar.minsma = ' + minsma
    ell_para[40] = 'geompar.maxsma = ' + maxsma
    ell_para[41] = 'geompar.step = ' + step
    ell_para[48] = 'magpar.mag0 = ' + mag0
    openw, 10, ellip_file 
    for k = 0, num_ellpara-1, 1 do begin
        printf, 10, ell_para[k]
    endfor
    close, 10
    ;; remove the previous ellipse results
    if file_test( output_bin ) then begin 
        spawn, 'rm  ' + output_bin
    endif
    if file_test( output_asc ) then begin 
        spawn, 'rm  ' + output_asc
    endif
    if file_test( output_head ) then begin 
        spawn, 'rm  ' + output_head
    endif
    ;; Run Ellipse 
    spawn, isophote + ' ellipse @' +  ellip_file
    ;; Test if the ellipse fitting goes well
    if NOT file_test( output_bin ) then begin 
        message, 'Somthing Wrong with the Ellipse Fitting, Check Again !'
    endif else begin 
    ;; change the binary table into ascii table
        print, 'Output Binary File: ' + output_bin
        num_tpara = file_lines( tdump_par )
        tdump_para = strarr( num_tpara )
        openr, 10, tdump_par
        readf, 10, tdump_para
        close, 10
        tdump_para[0] = 'tdump.table = "' +  output_bin + '"'
        tdump_para[2] = 'tdump.pfile = "' +  output_head + '"'
        tdump_para[3] = 'tdump.datafile = "' + output_asc + '"'
        openw, 10, 'tdump.par'
        for m = 0, num_tpara-1, 1 do begin
            printf, 10, tdump_para[m]
        endfor
        close, 10
        spawn, ttools + ' tdump @tdump.par'
        spawn, 'rm ' + output_head
        if file_test( output_asc ) then begin 
            ;spawn, 'sed -in-place -e "s/INDEF/-9999.9/g" ' + output_asc 
            spawn, 'sed -i "s/INDEF/-9999.9/g" ' + output_asc 
        endif
    endelse

endif

;; readin the surface brightness profile 
old_mod_prof_2 = output_asc
read_profile, old_mod_prof_2, old_mod_sbp_2, old_mod_line_2
 old_mod_sma_2 = old_mod_sbp_2.sma 
old_mod_rsma_2 = ( old_mod_sma_2 * pix )^0.25
old_mod_smag_2 = -2.5 * alog10( old_mod_sbp_2.intens / $
    ( pix_area * old_expt2 ) ) + magzpt2 - extinc_2
 old_mod_sig_2 = sqrt( ( old_mod_sbp_2.intens_err )^2.0 + ( sky_ner_2 )^2.0 ) 
old_mod_serr_2 = 2.5 * alog10( 1.0 + ( old_mod_sig_2 / $
    ( old_mod_sbp_2.intens ) ) )
;; get the color profile and color error profile for old_model images
if ( filter_index1 lt filter_index2 ) then begin 
    old_mod_cprof = old_mod_smag_1 - old_mod_smag_2 
endif else begin 
    old_mod_cprof = old_mod_smag_2 - old_mod_smag_1 
endelse
old_mod_cerr = sqrt( ( old_mod_serr_1 )^2.0 + ( old_mod_serr_2 )^2.0 )
old_mod_color_median = median( old_mod_cprof )
print, '############################################################'
print, 'Median old model Color : ' + string( old_mod_color_median )
print, '############################################################'
galaxy_color[3] = old_mod_color_median

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 7. ellipse run on model image without psf convolution for filter 1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the input files
 input_img = model_nosky_1
output_bin = model_string_loc_1 + '_nosky_mod.bin'
output_asc = model_string_loc_1 + '_nosky_mod.prof'
ellip_file = model_string_loc_1 + '_nosky_mod.ellip'
output_head = 'head.dat'

if NOT file_test( output_asc ) then begin 

    ;; get the necessary information
    x0 = strcompress( string( cen_x_1 ), /remove_all )
    y0 = strcompress( string( cen_y_1 ), /remove_all )
    ellip0 = strcompress( string( ell_1 ), /remove_all )
    pa0 = strcompress( string( pa_1 ), /remove_all )
    sma0 = strcompress( string( 30.0 ), /remove_all )
    minsma = strcompress( string( 0.0 ), /remove_all )
    maxsma = strcompress( string( 1080.0 ), /remove_all )
    step = strcompress( string( 0.05 ), /remove_all ) 
    mag0 = strcompress( string( magzpt1 ), /remove_all ) 
    ;; edit the ellipse parameter file
    num_ellpara = file_lines( ell_par ) 
    ell_para = strarr( num_ellpara )
    openr, 10, ell_par
    readf, 10, ell_para
    close, 10
    ell_para[0] = 'ellipse.input = "' + input_img + '"'
    ell_para[1] = 'ellipse.output = "' + output_bin + '"'
    ell_para[3] = 'ellipse.inellip = "' + standard_bin + '"'
    ell_para[27] = 'controlpar.hellip = yes '
    ell_para[28] = 'controlpar.hpa = yes '
    ell_para[34] = 'geompar.x0 = ' + x0
    ell_para[35] = 'geompar.y0 = ' + y0
    ell_para[36] = 'geompar.ellip0 = ' + ellip0
    ell_para[37] = 'geompar.pa0 = ' + pa0
    ell_para[38] = 'geompar.sma0 = ' + sma0
    ell_para[39] = 'geompar.minsma = ' + minsma
    ell_para[40] = 'geompar.maxsma = ' + maxsma
    ell_para[41] = 'geompar.step = ' + step
    ell_para[48] = 'magpar.mag0 = ' + mag0
    openw, 10, ellip_file 
    for k = 0, num_ellpara-1, 1 do begin
        printf, 10, ell_para[k]
    endfor
    close, 10
    ;; remove the previous ellipse results
    if file_test( output_bin ) then begin 
        spawn, 'rm  ' + output_bin
    endif
    if file_test( output_asc ) then begin 
        spawn, 'rm  ' + output_asc
    endif
    if file_test( output_head ) then begin 
        spawn, 'rm  ' + output_head
    endif
    ;; Run Ellipse 
    spawn, isophote + ' ellipse @' +  ellip_file
    ;; Test if the ellipse fitting goes well
    if NOT file_test( output_bin ) then begin 
        message, 'Somthing Wrong with the Ellipse Fitting, Check Again !'
    endif else begin 
    ;; change the binary table into ascii table
        print, 'Output Binary File: ' + output_bin
        num_tpara = file_lines( tdump_par )
        tdump_para = strarr( num_tpara )
        openr, 10, tdump_par
        readf, 10, tdump_para
        close, 10
        tdump_para[0] = 'tdump.table = "' +  output_bin + '"'
        tdump_para[2] = 'tdump.pfile = "' +  output_head + '"'
        tdump_para[3] = 'tdump.datafile = "' + output_asc + '"'
        openw, 10, 'tdump.par'
        for m = 0, num_tpara-1, 1 do begin
            printf, 10, tdump_para[m]
        endfor
        close, 10
        spawn, ttools + ' tdump @tdump.par'
        spawn, 'rm ' + output_head
        if file_test( output_asc ) then begin 
            ;spawn, 'sed -in-place -e "s/INDEF/-9999.9/g" ' + output_asc 
            spawn, 'sed -i "s/INDEF/-9999.9/g" ' + output_asc 
        endif
    endelse

endif
;; readin the surface brightness profile 
mod_prof_1 = output_asc
read_profile, mod_prof_1, mod_sbp_1, mod_line_1
 mod_sma_1 = mod_sbp_1.sma 
mod_rsma_1 = ( mod_sma_1 * pix )^0.25
mod_smag_1 = -2.5 * alog10( mod_sbp_1.intens / ( pix_area * old_expt1 ) ) + $
    magzpt1 - extinc_1
 mod_sig_1 = sqrt( ( mod_sbp_1.intens_err )^2.0 + ( sky_ner_1 )^2.0 ) 
mod_serr_1 = 2.5 * alog10( 1.0 + ( mod_sig_1 / ( mod_sbp_1.intens ) ) )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 8. ellipse run on model image wihtout psf convolution for filter 2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the input files
 input_img = model_nosky_2
output_bin = model_string_loc_2 + '_nosky_mod.bin'
output_asc = model_string_loc_2 + '_nosky_mod.prof'
ellip_file = model_string_loc_2 + '_nosky_mod.ellip'
output_head = 'head.dat'

if NOT file_test( output_asc ) then begin 

    ;; get the necessary information
    x0 = strcompress( string( cen_x_2 ), /remove_all )
    y0 = strcompress( string( cen_y_2 ), /remove_all )
    ellip0 = strcompress( string( ell_2 ), /remove_all )
    pa0 = strcompress( string( pa_2 ), /remove_all )
    sma0 = strcompress( string( 30.0 ), /remove_all )
    minsma = strcompress( string( 0.0 ), /remove_all )
    maxsma = strcompress( string( 1080.0 ), /remove_all )
    step = strcompress( string( 0.05 ), /remove_all ) 
    mag0 = strcompress( string( magzpt2 ), /remove_all ) 
    ;; edit the ellipse parameter file
    num_ellpara = file_lines( ell_par ) 
    ell_para = strarr( num_ellpara )
    openr, 10, ell_par
    readf, 10, ell_para
    close, 10
    ell_para[0] = 'ellipse.input = "' + input_img + '"'
    ell_para[1] = 'ellipse.output = "' + output_bin + '"'
    ell_para[3] = 'ellipse.inellip = "' + standard_bin + '"'
    ell_para[27] = 'controlpar.hellip = yes '
    ell_para[28] = 'controlpar.hpa = yes '
    ell_para[34] = 'geompar.x0 = ' + x0
    ell_para[35] = 'geompar.y0 = ' + y0
    ell_para[36] = 'geompar.ellip0 = ' + ellip0
    ell_para[37] = 'geompar.pa0 = ' + pa0
    ell_para[38] = 'geompar.sma0 = ' + sma0
    ell_para[39] = 'geompar.minsma = ' + minsma
    ell_para[40] = 'geompar.maxsma = ' + maxsma
    ell_para[41] = 'geompar.step = ' + step
    ell_para[48] = 'magpar.mag0 = ' + mag0
    openw, 10, ellip_file 
    for k = 0, num_ellpara-1, 1 do begin
        printf, 10, ell_para[k]
    endfor
    close, 10
    ;; remove the previous ellipse results
    if file_test( output_bin ) then begin 
        spawn, 'rm  ' + output_bin
    endif
    if file_test( output_asc ) then begin 
        spawn, 'rm  ' + output_asc
    endif
    if file_test( output_head ) then begin 
        spawn, 'rm  ' + output_head
    endif
    ;; Run Ellipse 
    spawn, isophote + ' ellipse @' +  ellip_file
    ;; Test if the ellipse fitting goes well
    if NOT file_test( output_bin ) then begin 
        message, 'Somthing Wrong with the Ellipse Fitting, Check Again !'
    endif else begin 
    ;; change the binary table into ascii table
        print, 'Output Binary File: ' + output_bin
        num_tpara = file_lines( tdump_par )
        tdump_para = strarr( num_tpara )
        openr, 10, tdump_par
        readf, 10, tdump_para
        close, 10
        tdump_para[0] = 'tdump.table = "' +  output_bin + '"'
        tdump_para[2] = 'tdump.pfile = "' +  output_head + '"'
        tdump_para[3] = 'tdump.datafile = "' + output_asc + '"'
        openw, 10, 'tdump.par'
        for m = 0, num_tpara-1, 1 do begin
            printf, 10, tdump_para[m]
        endfor
        close, 10
        spawn, ttools + ' tdump @tdump.par'
        spawn, 'rm ' + output_head
        if file_test( output_asc ) then begin 
            ;spawn, 'sed -in-place -e "s/INDEF/-9999.9/g" ' + output_asc 
            spawn, 'sed -i "s/INDEF/-9999.9/g" ' + output_asc 
        endif
    endelse

endif

;; readin the surface brightness profile 
mod_prof_2 = output_asc
read_profile, mod_prof_2, mod_sbp_2, mod_line_2
 mod_sma_2 = mod_sbp_2.sma 
mod_rsma_2 = ( mod_sma_2 * pix )^0.25
mod_smag_2 = -2.5 * alog10( mod_sbp_2.intens / $
    ( pix_area * old_expt2 ) ) + magzpt2 - extinc_2
 mod_sig_2 = sqrt( ( mod_sbp_2.intens_err )^2.0 + ( sky_ner_2 )^2.0 ) 
mod_serr_2 = 2.5 * alog10( 1.0 + ( mod_sig_2 / $
    ( mod_sbp_2.intens ) ) )
;; get the color profile and color error profile for model images
if ( filter_index1 lt filter_index2 ) then begin 
    mod_cprof = mod_smag_1 - mod_smag_2 
endif else begin 
    mod_cprof = mod_smag_2 - mod_smag_1 
endelse
mod_cerr = sqrt( ( mod_serr_1 )^2.0 + ( mod_serr_2 )^2.0 )
mod_color_median = median( mod_cprof )
print, '############################################################'
print, 'Median model Color (No PSF): ' + string( mod_color_median )
print, '############################################################'
galaxy_color[4] = mod_color_median

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; get the residual profile of these two model 
old_res_1 = ( old_ori_smag_1 - old_mod_smag_1 ) 
old_res_2 = ( old_ori_smag_2 - old_mod_smag_2 ) 
res_prof_1 = ( ori_smag_1 - mod_smag_1 ) 
res_prof_2 = ( ori_smag_2 - mod_smag_2 ) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Save a summary file for this color analysis 
color_sum_file = galaxy + '_' + filter1 + '_' + filter2 + '_' + ncomp + '_' + $
    type + '_color_sum.dat'
openw, 20, color_sum_file, width = 900 
printf, 20, '##  MODEL        ' + filter1 + ' : ' + model_string_1 
printf, 20, '##  MODEL        ' + filter2 + ' : ' + model_string_2
printf, 20, '##  SKY_VAL/_ERR ' + filter1 + ' : ' + string( sky_new_1 ) + $
    ' +/- ' + string( sky_ner_1 ) 
printf, 20, '##  SKY_VAL/_ERR ' + filter2 + ' : ' + string( sky_new_2 ) + $
    ' +/- ' + string( sky_ner_2 ) 
printf, 20, '##  Extinction   ' + filter1 + ' : ' + string( extinc_1 )
printf, 20, '##  Extinction   ' + filter2 + ' : ' + string( extinc_2 )
printf, 20, '##  Total Mag    ' + filter1 + ' : ' + string( mag_total_1 ) 
printf, 20, '##  Total Mag    ' + filter2 + ' : ' + string( mag_total_2 ) 
printf, 20, '##  Avg Color    ' + string( color_average )
for i = 0, ( ncomp_1 - 1 ), 1 do begin 
    index = strcompress( string( i ), /remove_all ) 
    printf, 20, '## Color Component ' + index + ' : ' + $
        string( comp_color[i] ) 
endfor
printf, 20, '## Median Color ORI (OLD) : ' + string( old_ori_color_median ) 
printf, 20, '## Median Color MOD (OLD) : ' + string( old_mod_color_median ) 
printf, 20, '## Median Color ORI (CONV)   : ' + string( ori_color_median ) 
printf, 20, '## Median Color MOD (NO PSF) : ' + string( mod_color_median ) 
 key1 = 'SMA   '
 key2 = 'RSMA  '
 key3 = 'ORI_CONV_SBP_' + filter1 
 key4 = 'ORI_CONV_SBP_ERR_' + filter1 
 key5 = 'ORI_CONV_SBP_' + filter2 
 key6 = 'ORI_CONV_SBP_ERR_' + filter2 
 key7 = 'MOD_NOPSF_SBP_' + filter1 
 key8 = 'MOD_NOPSF_SBP_ERR_' + filter1 
 key9 = 'MOD_NOPSF_SBP_' + filter2 
key10 = 'MOD_NOPSF_SBP_ERR_' + filter2 
key11 = 'COLOR_ORI' 
key12 = 'COLOR_ERR_ORI'
key13 = 'COLOR_MOD'
key14 = 'COLOR_ERR_MOD'
key15 = 'ORI_OLD_SBP_' + filter1
key16 = 'ORI_OLD_SBP_ERR_' + filter1
key17 = 'ORI_OLD_SBP_' + filter2
key18 = 'ORI_OLD_SBP_ERR_' + filter2
key19 = 'MOD_OLD_SBP_' + filter1
key20 = 'MOD_OLD_SBP_ERR_' + filter1
key21 = 'MOD_OLD_SBP_' + filter2
key22 = 'MOD_OLD_SBP_ERR_' + filter2
key23 = 'OLD_COLOR_ORI' 
key24 = 'OLD_COLOR_ERR_ORI'
key25 = 'OLD_COLOR_MOD'
key26 = 'OLD_COLOR_ERR_MOD'
key27 = 'OLD_RESIDUL_' + filter1 
key28 = 'OLD_RESIDUL_' + filter2 

printf, 20, '#' + key1 + '  ' + key2 + '  ' + key3 + '  ' + $ 
    key4 + '  ' + key5 + '  ' + key6 + '  ' + key7 + '  ' + key8 + '  ' + $ 
    key9 + '  ' + key10 + '  ' + key11 + '  ' + key12 + '  ' + key13 + '  ' + $ 
    key14 + '  ' + key15 + '  ' + key16 + '  ' + key17 + '  ' + key18 + $
    key19 + '  ' + key20 + '  ' + key21 + '  ' + key22 + '  ' + key23 + $ 
    key24 + '  ' + key25 + '  ' + key26 + '  ' + key27 + '  ' + key28
space = '   '
for i  = 0, ( ori_line_1 - 1 ), 1 do begin 
    val1 = string( ori_sma_1[i], FORMAT='(F8.3)' )
    val2 = string( ori_rsma_1[i], FORMAT='(F8.3)' )
    val3 = string( ori_smag_1[i], FORMAT='(F8.3)' )
    val4 = string( ori_serr_1[i], FORMAT='(F8.3)' )
    val5 = string( ori_smag_2[i], FORMAT='(F8.3)' )
    val6 = string( ori_serr_2[i], FORMAT='(F8.3)' )
    val7 = string( mod_smag_1[i], FORMAT='(F8.3)' )
    val8 = string( mod_serr_1[i], FORMAT='(F8.3)' )
    val9 = string( mod_smag_2[i], FORMAT='(F8.3)' )
   val10 = string( mod_serr_2[i], FORMAT='(F8.3)' )
   val11 = string( ori_cprof[i], FORMAT='(F8.3)' )
   val12 = string( ori_cerr[i], FORMAT='(F8.3)' )
   val13 = string( mod_cprof[i], FORMAT='(F8.3)' )
   val14 = string( mod_cerr[i], FORMAT='(F8.3)' )

   val15 = string( old_ori_smag_1[i], FORMAT='(F8.3)' )
   val16 = string( old_ori_serr_1[i], FORMAT='(F8.3)' )
   val17 = string( old_ori_smag_2[i], FORMAT='(F8.3)' )
   val18 = string( old_ori_serr_2[i], FORMAT='(F8.3)' )
   val19 = string( old_mod_smag_1[i], FORMAT='(F8.3)' )
   val20 = string( old_mod_serr_1[i], FORMAT='(F8.3)' )
   val21 = string( old_mod_smag_2[i], FORMAT='(F8.3)' )
   val22 = string( old_mod_serr_2[i], FORMAT='(F8.3)' )
   val23 = string( old_ori_cprof[i], FORMAT='(F8.3)' )
   val24 = string( old_ori_cerr[i], FORMAT='(F8.3)' )
   val25 = string( old_mod_cprof[i], FORMAT='(F8.3)' )
   val26 = string( old_mod_cerr[i], FORMAT='(F8.3)' )

   val27 = string( old_res_1[i], FORMAT='(F8.3)' )
   val28 = string( old_res_2[i], FORMAT='(F8.3)' )

   printf, 20, val1 + space + val2 + space + val3 + space + val4 + space + $ 
       val5 + space + val6 + space + val7 + space + val8 + space + val9 + $ 
       space + val10 + space + val11 + space + val12 + space + val13 + space + $ 
       val14 + space + val15 + space + val16 + space + val17 + space + val18 + $
       space + val19 + space + val20 + space + val21 + space + val22 + space + $ 
       val23 + space + val24 + space + val25 + space + val26 + space + val27 + $ 
       space + val28

endfor
close, 20

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; make a surface brightness profile and color profile comparison plot
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
color_sum_plot = galaxy + '_' + filter1 + '_' + filter2 + '_' + ncomp + '_' + $
    type + '_color_sum.ps'
;; define plotting range 
;; rsma range 
min_rsma = 0.50
if ( old_prof_exist eq 1 ) then begin 
    max_rsma = max( ( old_sbp.sma )^0.25 ) + 0.25
endif else begin 
    max_rsma = max( ori_rsma_1 )
endelse
xrange = [ min_rsma, max_rsma ]
;; sbp range 
smag_lim = where( ori_rsma_1 le max_rsma )
min_smag = min( [ min( ori_smag_1 ), min( ori_smag_1), min( mod_smag_1 ), $
    min( mod_smag_2 ) ] )
min_smag = min_smag > 10.0
max_smag = max( [ max( old_ori_smag_1[ smag_lim ] ), $
    max( old_ori_smag_2[ smag_lim ] ), max( old_mod_smag_1[ smag_lim ] ), $ 
    max( old_mod_smag_1[ smag_lim ] ) ] )
max_smag = max_smag < 29.5
yrange1 = [ max_smag + 0.9, min_smag - 0.1 ]
;; color range 
min_color = min( mod_cprof[ smag_lim ] )
min_color = min_color < min( comp_color )
min_color = min_color - 0.05 
min_color = min_color > 0.01
max_color = max( ori_cprof[ where( ori_rsma_1 le 2.0 ) ] ) > $
    max( mod_cprof[ where( ori_rsma_1 le 2.0 ) ] )
max_color = max_color > max( comp_color )
max_color = max_color + 0.07
yrange2 = [ min_color, max_color ]
;; residual range
yrange3 = [ -0.19999, 0.19999 ]
;; define plotting regions 
position_1 = [ 0.18, 0.18, 0.98, 0.60 ]
position_2 = [ 0.18, 0.60, 0.98, 0.80 ]
position_4 = [ 0.18, 0.80, 0.98, 0.99 ]
position_3 = [ 0.18, 0.07, 0.98, 0.18 ]
;; define the title 
 xtitle = textoidl( 'R^{1/4} (arcsec^{1/4})' )
ytitle1 = textoidl( '\mu (mag/arcsec^2)' )
if ( filter_index1 lt filter_index2 ) then begin 
    ytitle2 = filter1 + '-' + filter2 + ' Color'
endif else begin 
    ytitle2 = filter2 + '-' + filter1 + ' Color'
endelse
ytitle3 = 'Residual'
;; define the size of the image 
psxsize = 20 
psysize = 26
;; plot the surface brightness profile 
mydevice = !D.NAME
set_plot, 'PS'
device, filename=color_sum_plot, font_size=8, /encapsul, $
    /color, set_font='HELVETICA BOLD', /tt_font, xsize=psxsize,$
    ysize=psysize
cgPlot, ori_rsma_1, ori_smag_1, xstyle=1, ystyle=1, xrange=xrange, $
    yrange=yrange1, /nodata, ytitle=ytitle1, position=position_1, $
    /noerase, xtickinterval=1.0, charsize=2.2, $
    charthick=3.8, xthick=4.3, xtickformat="(A1)", $
    ythick=4.3,  ytickinterval=3.0 
if ( filter_index1 lt filter_index2 ) then begin 
    cgPlot, old_ori_rsma_1, old_ori_smag_1, psym=9, symsize=1.1, thick=3.5, $
        color=cgColor( 'CYAN', !D.Table_Size ), /overplot
    oploterror, old_ori_rsma_1, old_ori_smag_1, old_ori_serr_1, psym=3, $
        errcolor=cgColor('CYAN', !D.Table_Size), errthick=3.5
    cgPlot, old_mod_rsma_1, old_ori_smag_1, psym=16, symsize=0.9, thick=3.5, $
        color=cgColor( 'NAVY', !D.Table_Size ), /overplot
    cgPlot, ori_rsma_1, ori_smag_1, psym=0, linestyle=1, thick=4.0, $
        color=cgColor('BLU4', !D.Table_Size), /overplot
    cgPlot, mod_rsma_1, mod_smag_1, psym=0, linestyle=2, thick=3.8, $
        color=cgColor('BLU7', !D.Table_Size), /overplot

    cgPlot, old_ori_rsma_2, old_ori_smag_2, psym=4, symsize=1.3, thick=3.5, $
        color=cgColor( 'YELLOW', !D.Table_Size ), /overplot
    oploterror, old_ori_rsma_2, old_ori_smag_2, old_ori_serr_2, psym=3, $
        errcolor=cgColor('YELLOW', !D.Table_Size), errthick=3.5
    cgPlot, old_mod_rsma_2, old_ori_smag_2, psym=14, symsize=1.2, thick=3.5, $
        color=cgColor( 'RED4', !D.Table_Size ), /overplot
    cgPlot, ori_rsma_2, ori_smag_2, psym=0, linestyle=3, thick=4.0, $
        color=cgColor('ORANGE', !D.Table_Size), /overplot
    cgPlot, mod_rsma_2, mod_smag_2, psym=0, linestyle=4, thick=3.8, $
        color=cgColor('RED7', !D.Table_Size), /overplot
endif else begin 
    cgPlot, old_ori_rsma_2, old_ori_smag_2, psym=9, symsize=1.1, thick=3.5, $
        color=cgColor( 'CYAN', !D.Table_Size ), /overplot
    oploterror, old_ori_rsma_2, old_ori_smag_2, old_ori_serr_2, psym=3, $
        errcolor=cgColor('CYAN', !D.Table_Size), errthick=3.5
    cgPlot, old_mod_rsma_2, old_mod_smag_2, psym=16, symsize=0.9, thick=3.5, $
        color=cgColor( 'NAVY', !D.Table_Size ), /overplot
    cgPlot, ori_rsma_2, ori_smag_2, psym=0, linestyle=1, thick=4.0, $
        color=cgColor('BLU4', !D.Table_Size), /overplot
    cgPlot, mod_rsma_2, mod_smag_2, psym=0, linestyle=2, thick=3.8, $
        color=cgColor('BLU7', !D.Table_Size), /overplot

    cgPlot, old_ori_rsma_1, old_ori_smag_1, psym=4, symsize=1.3, thick=3.5, $
        color=cgColor( 'YELLOW', !D.Table_Size ), /overplot
    oploterror, old_ori_rsma_1, old_ori_smag_1, old_ori_serr_1, psym=3, $
        errcolor=cgColor('YELLOW', !D.Table_Size), errthick=3.5
    cgPlot, old_mod_rsma_1, old_mod_smag_1, psym=14, symsize=1.2, thick=3.5, $
        color=cgColor( 'RED4', !D.Table_Size ), /overplot
    cgPlot, ori_rsma_1, ori_smag_1, psym=0, linestyle=3, thick=4.0, $
        color=cgColor('ORANGE', !D.Table_Size), /overplot
    cgPlot, mod_rsma_1, mod_smag_1, psym=0, linestyle=4, thick=3.8, $
        color=cgColor('RED7', !D.Table_Size), /overplot
endelse
cgText, 0.76, 0.54, galaxy, charsize=4.2, alignment=0.5, $
    charthick=8.0, /normal, color=cgColor( 'BLACK', !D.Table_Size )
model_index = strcompress( 'Model: ' + ncomp + '0_' + type, /remove_all )
cgText, 0.76, 0.49, model_index, charsize=2.6, alignment=0.5, charthick=4.0, $
    /normal
if ( filter_index1 lt filter_index2 ) then begin 
    avg_color_string = '<' + filter1 + '-' + filter2 +'>=' + $
        string( color_average, format='(F5.2)' )
endif else begin 
    avg_color_string = '<' + filter2 + '-' + filter1 +'>=' + $
        strcompress( string( color_average, format='(F5.2)' ), /remove_all )
endelse
cgText, 0.76, 0.45, avg_color_string, charsize=2.6, alignment=0.5, $
    charthick=4.0, /normal, color=cgColor( 'BLACK', !D.Table_Size )
;; label
cgPlots, 0.26, 0.425, psym=9, symsize=1.1, thick=3.5, /normal, $ 
    color=cgColor( 'CYAN', !D.Table_Size ) 
cgText, 0.32, 0.42, blue_filter + ' data', charsize=1.8, charthick=3.5, $
    alignment=0, /normal 
cgPlots, 0.26, 0.395, psym=16, symsize=1.0, thick=3.5, /normal, $ 
    color=cgColor( 'NAVY', !D.Table_Size ) 
cgText, 0.32, 0.39, blue_filter + ' model', charsize=1.8, charthick=3.5, $
    alignment=0, /normal 
cgPlots, 0.26, 0.365, psym=4, symsize=1.3, thick=3.5, /normal, $ 
    color=cgColor( 'YELLOW', !D.Table_Size ) 
cgText, 0.32, 0.36, red_filter + ' data', charsize=1.8, charthick=3.5, $
    alignment=0, /normal 
cgPlots, 0.26, 0.335, psym=14, symsize=1.3, thick=3.5, /normal, $ 
    color=cgColor( 'RED4', !D.Table_Size ) 
cgText, 0.32, 0.33, red_filter + ' model', charsize=1.8, charthick=3.5, $
    alignment=0, /normal 

cgPlots, [ 0.22, 0.30 ], [ 0.305, 0.305 ], psym=0, linestyle=1, thick=2.6, $
    color=cgColor( 'BLU4', !D.Table_Size ), /normal
cgText, 0.32, 0.30, blue_filter + ' conv. data', charsize=1.9, charthick=3.5, $
    alignment=0, /normal 
cgPlots, [ 0.22, 0.30 ], [ 0.275, 0.275 ], psym=0, linestyle=2, thick=2.6, $
    color=cgColor( 'BLU7', !D.Table_Size ), /normal
cgText, 0.32, 0.27, blue_filter + ' nopsf model', charsize=1.9, charthick=3.5, $
    alignment=0, /normal 

cgPlots, [ 0.22, 0.30 ], [ 0.245, 0.245 ], psym=0, linestyle=3, thick=2.6, $
    color=cgColor( 'ORANGE', !D.Table_Size ), /normal
cgText, 0.32, 0.24, red_filter + ' conv. data', charsize=1.9, charthick=3.5, $
    alignment=0, /normal 
cgPlots, [ 0.22, 0.30 ], [ 0.215, 0.215 ], psym=0, linestyle=4, thick=2.6, $
    color=cgColor( 'RED7', !D.Table_Size ), /normal
cgText, 0.32, 0.21, red_filter + ' nopsf model', charsize=1.9, charthick=3.5, $
    alignment=0, /normal 

re_array = ( comp_sum1.r )^0.25  
for i = 0, ( ncomp_1 - 1 ), 1 do begin 
    mu1 = comp_mue1[i] 
    mu2 = comp_mue2[i] 
    r = re_array[i]
    symsize = ( comp_symsize[i] * 1.2 ) 
    if ( filter_index1 lt filter_index2 ) then begin 
        cgPlots, r, mu1, psym=28, symsize=symsize, $
            color=cgColor( 'BLU8', !D.Table_Size ), /data
        cgPlots, r, mu2, psym=28, symsize=symsize, $
            color=cgColor( 'RED5', !D.Table_Size ), /data
    endif else begin 
        cgPlots, r, mu2, psym=28, symsize=symsize, $
            color=cgColor( 'BLU8', !D.Table_Size ), /data
        cgPlots, r, mu1, psym=28, symsize=symsize, $
            color=cgColor( 'RED5', !D.Table_Size ), /data
    endelse
endfor

;; plot the color profile ( new ) 
if ( ( yrange2[1] - yrange2[0] ) ge 0.5 ) then begin 
    if ( ( yrange2[1] - yrange2[0] ) lt 0.8 ) then begin 
        ytick = 0.2 
    endif else begin 
        if ( ( yrange2[1] - yrange2[0] ) lt 1.1 ) then begin 
            ytick = 0.3
        endif else begin 
            ytick = 0.4 
        endelse
    endelse
endif else begin 
    ytick = 0.1
endelse
    
cgPlot, ori_rsma_1, ori_cprof, xstyle=1, ystyle=1, xrange=xrange, $
    yrange=yrange2, /nodata, ytitle=ytitle2, position=position_2, $
    /noerase, xtickinterval=1.0, charsize=2.2, $
    charthick=3.8, xthick=4.3, xtickformat="(A1)",  $
    ythick=4.3,  yminor=-1, ytickinterval=ytick
cgPlot, ori_rsma_1, ori_cprof, psym=16, symsize=1.1, thick=3.5, $
    color=cgColor( 'BLACK', !D.Table_Size ), /overplot
oploterror, ori_rsma_1, ori_cprof, ori_cerr, psym=3, $
    errcolor=cgColor('BLACK', !D.Table_Size), errthick=3.5
cgPlot, mod_rsma_1, mod_cprof, psym=9, symsize=1.2, thick=3.1, $
    color=cgColor( 'RED', !D.Table_Size ), /overplot
oploterror, mod_rsma_1, mod_cprof, mod_cerr, psym=3, $
    errcolor=cgColor('RED', !D.Table_Size), errthick=3.4
cgPlot, !X.Crange, [ori_color_median, ori_color_median], psym=0, linestyle=2, $
    thick=3.9, color=cgColor( 'DARK GRAY', !D.Table_Size ), /overplot
for i = 0, ( ncomp_1 - 1 ), 1 do begin 
    c = comp_color[i] 
    r = re_array[i]
    symsize = ( comp_symsize[i] * 1.2 ) 
    cgPlots, r, c, psym=28, symsize=symsize, $ 
        color=cgColor( 'NAVY', !D.Table_Size ), /data
endfor
color_middle = ( ( yrange2[1] - yrange2[0] ) / 2.0 ) + yrange2[0] 
if ( color_average lt color_middle ) then begin 
    cgText, 0.42, 0.780, 'solid : convol data', charsize=1.5, charthick=3.3, $
        alignment=0, /normal
    cgText, 0.42, 0.765, 'open : nopsf model', charsize=1.5, charthick=3.3, $
        alignment=0, /normal, color=cgColor( 'RED', !D.Table_Size )
endif else begin 
    cgText, 0.40, 0.609, 'solid : convol data', charsize=1.5, charthick=3.3, $
        alignment=0, /normal
    cgText, 0.40, 0.624, 'open : nopsf model', charsize=1.5, charthick=3.3, $
        alignment=0, /normal, color=cgColor( 'RED', !D.Table_Size )
endelse

;; plot the color profile ( old ) 
cgPlot, old_ori_rsma_1, old_ori_cprof, xstyle=1, ystyle=1, xrange=xrange, $
    yrange=yrange2, /nodata, ytitle=ytitle2, position=position_4, $
    /noerase, xtickinterval=1.0, charsize=2.2, $
    charthick=3.8, xthick=4.3, xtickformat="(A1)",  $
    ythick=4.3,  yminor=-1, ytickinterval=ytick
cgPlot, old_ori_rsma_1, old_ori_cprof, psym=14, symsize=1.2, thick=3.5, $
    color=cgColor( 'BLACK', !D.Table_Size ), /overplot
oploterror, old_ori_rsma_1, old_ori_cprof, old_ori_cerr, psym=3, $
    errcolor=cgColor('BLACK', !D.Table_Size), errthick=3.5
cgPlot, old_mod_rsma_1, old_mod_cprof, psym=4, symsize=1.3, thick=3.1, $
    color=cgColor( 'RED', !D.Table_Size ), /overplot
oploterror, old_mod_rsma_1, old_mod_cprof, old_mod_cerr, psym=3, $
    errcolor=cgColor('RED', !D.Table_Size), errthick=3.4
cgPlot, !X.Crange, [ori_color_median, ori_color_median], psym=0, linestyle=2, $
    thick=3.9, color=cgColor( 'DARK GRAY', !D.Table_Size ), /overplot
re_array = ( comp_sum1.r )^0.25  
for i = 0, ( ncomp_1 - 1 ), 1 do begin 
    c = comp_color[i] 
    r = re_array[i]
    symsize = ( comp_symsize[i] * 1.2 ) 
    cgPlots, r, c, psym=28, symsize=symsize, $ 
        color=cgColor( 'NAVY', !D.Table_Size ), /data
endfor 
color_middle = ( ( yrange2[1] - yrange2[0] ) / 2.0 ) + yrange2[0] 
if ( color_average lt color_middle ) then begin 
    cgText, 0.42, 0.970, 'solid : data', charsize=1.5, charthick=3.3, $
        alignment=0, /normal
    cgText, 0.42, 0.955, 'open : model', charsize=1.5, charthick=3.3, $
        alignment=0, /normal, color=cgColor( 'RED', !D.Table_Size )
endif else begin 
    cgText, 0.42, 0.809, 'solid : data', charsize=1.5, charthick=3.3, $
        alignment=0, /normal
    cgText, 0.42, 0.824, 'open : model', charsize=1.5, charthick=3.3, $
        alignment=0, /normal, color=cgColor( 'RED', !D.Table_Size )
endelse

;; plot the residual profile 
cgPlot, old_ori_rsma_1, old_res_1, xstyle=1, ystyle=1, xrange=xrange, $
    yrange=yrange3, /nodata, ytitle=ytitle3, position=position_3, $
    /noerase, xtickinterval=0.5, charsize=2.2, $
    charthick=3.8, xthick=4.3, $
    ythick=4.3,  yminor=-1, xtitle=xtitle
if ( filter_index1 lt filter_index2 ) then begin 
    cgPlot, old_ori_rsma_1, old_res_1, psym=0, linestyle=2, thick=4.5, $ 
        color=cgColor( 'BLUE', !D.Table_Size ), /overplot 
    cgPlot, old_ori_rsma_2, old_res_2, psym=0, linestyle=3, thick=4.5, $ 
        color=cgColor( 'RED', !D.Table_Size ), /overplot
endif else begin 
    cgPlot, old_ori_rsma_1, old_res_1, psym=0, linestyle=3, thick=4.5, $ 
        color=cgColor( 'RED', !D.Table_Size ), /overplot 
    cgPlot, old_ori_rsma_2, old_res_2, psym=0, linestyle=2, thick=4.5, $ 
        color=cgColor( 'BLUE', !D.Table_Size ), /overplot
endelse
cgPlot, !X.Crange, [0.00,0.00], psym=0, linestyle=0, thick=3.5, /overplot, $
    color=cgColor( 'DARK GRAY', !D.Table_Size )
cgPlot, ori_rsma_1, abs( ori_serr_1 ), psym=0, linestyle=4, thick=3.0, $
    /overplot, color=cgColor( 'Dark GRAY', !D.Table_Size )
cgPlot, ori_rsma_1, -1.0 * abs( ori_serr_1 ), psym=0, linestyle=4, thick=3.0, $
    /overplot, color=cgColor( 'Dark GRAY', !D.Table_Size )

device, /close
set_plot, mydevice
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set( color_map ) then begin
    ;; readin mask image 
    position_4 = [ 0.100, 0.15, 0.475, 0.85 ]
    position_5 = [ 0.485, 0.15, 0.860, 0.85 ]
    position_6 = [ 0.870, 0.15, 0.995, 0.85 ]

    rebin_factor = rebin_factor 
    new_area = ( pix * rebin_factor * pix * rebin_factor )
    img_msk_1 = mrdfits( msk_1, 0 ) 
    img_msk_2 = mrdfits( msk_2, 0 )
    temp = size( img_ori_conv_1, /dimension ) 
    naxis1 = temp[0] 
    naxis2 = temp[1]
    newx = long( naxis1 / rebin_factor )
    newy = long( naxis2 / rebin_factor )

    mue_lim1 = max( comp_mue_1 ) + 3.5 
    mue_lim2 = max( comp_mue_2 ) + 3.5 

    ;; original image color map
    if ( rebin_factor gt 1 ) then begin 
        img_ori_rebin_1 = frebin( img_ori_nosky_1, newx, newy, /total )
        img_ori_rebin_2 = frebin( img_ori_nosky_2, newx, newy, /total )
    endif else begin 
        img_ori_rebin_1 = img_ori_nosky_1
        img_ori_rebin_2 = img_ori_nosky_2
    endelse
    smag_ori_1 = -2.5 * alog10( img_ori_rebin_1 / ( old_expt1 * new_area ) ) $
        + magzpt1 - extinc_1
    smag_ori_2 = -2.5 * alog10( img_ori_rebin_2 / ( old_expt2 * new_area ) ) $
        + magzpt2 - extinc_2
    nan_mask_1 = finite( smag_ori_1, /NaN ) 
    nan_mask_2 = finite( smag_ori_2, /NaN )
    if ( filter_index1 lt filter_index2 ) then begin 
        ori_color = ( smag_ori_1 - smag_ori_2 )
    endif else begin 
        ori_color = ( smag_ori_2 - smag_ori_1 )
    endelse

    ;; add Poisson noise to the model image 
    img_mod_nosky_noise_1 = poidev( ( img_mod_nosky_1 * gain ), SEED=1.0 ) + $ 
        img_mod_nosky_1 
    img_mod_nosky_noise_2 = poidev( ( img_mod_nosky_2 * gain ), SEED=1.0 ) + $ 
        img_mod_nosky_2

    ;; model image color map
    if ( rebin_factor gt 1 ) then begin 
        img_mod_rebin_1 = frebin( img_mod_nosky_1, newx, newy, /total )
        img_mod_rebin_2 = frebin( img_mod_nosky_2, newx, newy, /total )
        img_mod_rebin_noise_1 = frebin( img_mod_nosky_noise_1, newx, newy, $
            /total )
        img_mod_rebin_noise_2 = frebin( img_mod_nosky_noise_2, newx, newy, $
            /total )
    endif else begin 
        img_mod_rebin_1 = img_mod_nosky_1 
        img_mod_rebin_2 = img_mod_nosky_2 
        img_mod_rebin_noise_1 = img_mod_nosky_noise_1 
        img_mod_rebin_noise_2 = img_mod_nosky_noise_2 
    endelse
    smag_mod_1 = -2.5 * alog10( img_mod_rebin_1 / ( old_expt1 * new_area ) ) $
        + magzpt1 - extinc_1
    smag_mod_2 = -2.5 * alog10( img_mod_rebin_2 / ( old_expt2 * new_area ) ) $
        + magzpt2 - extinc_2
    smag_mod_noise_1 = -2.5 * alog10( img_mod_rebin_noise_1 / ( old_expt1 * $
        new_area ) ) + magzpt1 - extinc_1
    smag_mod_noise_2 = -2.5 * alog10( img_mod_rebin_noise_2 / ( old_expt2 * $ 
        new_area ) ) + magzpt2 - extinc_2
    if ( filter_index1 lt filter_index2 ) then begin 
        mod_color = ( smag_mod_1 - smag_mod_2 )
        mod_color_noise = ( smag_mod_noise_1 - smag_mod_noise_2 )
    endif else begin 
        mod_color = ( smag_mod_2 - smag_mod_1 )
        mod_color_noise = ( smag_mod_noise_2 - smag_mod_noise_1 )
    endelse

    ;; whether or not use the noise added model image for color analysis 
    if keyword_set( add_noise ) then begin 
        mod_color = mod_color_noise 
    endif 

    ;; rebin the mask image 
    if ( rebin_factor gt 1 ) then begin
        img_msk_rebin_1 = frebin( img_msk_1, newx, newy, /total )
        img_msk_rebin_2 = frebin( img_msk_2, newx, newy, /total )
        msk_thre = long( rebin_factor * rebin_factor )
        img_msk_rebin_1[ where( img_msk_rebin_1 lt msk_thre ) ] = 0
        img_msk_rebin_2[ where( img_msk_rebin_2 lt msk_thre ) ] = 0
        img_msk_rebin_1[ where( img_msk_rebin_1 ge msk_thre ) ] = 1
        img_msk_rebin_2[ where( img_msk_rebin_2 ge msk_thre ) ] = 1
    endif else begin 
        img_msk_rebin_1 = img_msk_1 
        img_msk_rebin_2 = img_msk_2
    endelse

    ;; save the color image (using the noise added version !) 
    ori_color_file = galaxy + '_' + filter1 + '_' + filter2 + '_' + ncomp + $
        '_' + type + '_ori_color.fits'
    mwrfits, ori_color, ori_color_file, /create 
    mod_color_file = galaxy + '_' + filter1 + '_' + filter2 + '_' + ncomp + $
        '_' + type + '_mod_color.fits'
    mwrfits, mod_color_noise, mod_color_file, /create 

    ;; define a region that contains the useful part of galaxy 
    ab_ratio = 1.00 / ( 1.00 - ell_1 )
    new_cen_x = ( cen_x_1 / rebin_factor )
    new_cen_y = ( cen_y_1 / rebin_factor )
    dist_ellipse, ellip_mask, [ newx, newy ], new_cen_x, new_cen_y, ab_ratio, $
        pa_1, /double
    rlim = ( r80_1 * expand_factor ) < naxis1
    rlim = ( rlim / rebin_factor )

    ;; select the useful pixels
    print, '########################################################'
    print, 'R_LIMIT : ' + string( rlim / rebin_factor ) + ' pixels'
    print, 'R_LIMIT : ' + string( rlim * rebin_factor * 0.259 ) + ' arcsec'
    if ( filter_index1 gt filter_index2 ) then begin 
        ;index_galaxy = ( where( smag_ori_1 ge min( ori_smag_1 + 0.2 ) ) and $
        ;    where( smag_ori_1 le max( ori_smag_1 + 1.2 ) ) and $ 
        ;    where( ori_color le ( max( ori_color ) + 0.2 ) ) and $ 
        ;    where( ori_color ge ( min( ori_color ) - 0.2 ) ) $ 
        ;    and where( img_msk_rebin_1 eq 0 ) and where( img_msk_rebin_2 eq 0 ) )  
            ;)
        index_galaxy = where( ( ellip_mask le rlim ) and $
            ( nan_mask_1 eq 0 ) and ( nan_mask_2 eq 0 ) and $
            ( img_msk_rebin_1 eq 0 ) ) 
        index_inregion = where( ellip_mask le rlim )
        index_maskout  = where( ( ellip_mask le rlim ) and $ 
            ( nan_mask_1 eq 0 ) and ( nan_mask_2 eq 0 ) and $
            ( img_msk_rebin_1 ne 0 ) )

    endif else begin 
        ;index_galaxy = ( where( smag_ori_2 ge min( ori_smag_2 + 0.2 ) ) and $
        ;    where( smag_ori_2 le max( ori_smag_2 + 1.2 ) ) and $ 
        ;    where( ori_color le ( max( ori_color ) + 0.2 ) ) and $ 
        ;    where( ori_color ge ( min( ori_color ) - 0.2 ) ) $
        ;    and where( img_msk_rebin_1 eq 0 ) and where( img_msk_rebin_2 eq 0 ) )  
            ;) 
        index_galaxy = where( ( ellip_mask le rlim ) and $
            ( nan_mask_1 eq 0 ) and ( nan_mask_2 eq 0 ) and $
            ( img_msk_rebin_1 eq 0 ) ) 
        index_inregion = where( ellip_mask le rlim )
        index_maskout  = where( ( ellip_mask le rlim ) and $ 
            ( nan_mask_1 eq 0 ) and ( nan_mask_2 eq 0 ) and $
            ( img_msk_rebin_1 ne 0 ) )
    endelse
    num_pixels = n_elements( index_galaxy )
    num_inregion = n_elements( index_inregion )
    num_maskout = n_elements( index_maskout )
    frac_maskout = ( num_maskout / num_inregion ) * 100.0 
    ;; create a color map array 
    print, '########################################################'
    print, 'Number of Total Pixels : ', ( newx * newy )
    print, 'Number of Useful Pixels : ', num_pixels 
    print, 'Number of Pixels in useful reiong : ' , num_inregion 
    print, 'Number of Pixels have been masked out : ' , num_maskout
    print, 'Mask Fraction : ' , frac_maskout
    print, '########################################################'
    cmap_array = make_array( num_pixels, 7, /float ) 
    for k = 0L, ( num_pixels - 1 ), 1L do begin 
        index = index_galaxy[ k ] 
        cmap_array[k,0] = ellip_mask[index]
        cmap_array[k,1] = smag_ori_1[index]
        cmap_array[k,2] = smag_ori_2[index]
        cmap_array[k,3] = smag_mod_1[index]
        cmap_array[k,4] = smag_mod_2[index]
        cmap_array[k,5] = ori_color[index]
        cmap_array[k,6] = mod_color[index]
    endfor
    print, min( cmap_array[*,1] ), max( cmap_array[*,1] ) 
    print, min( cmap_array[*,2] ), max( cmap_array[*,2] )
    print, min( cmap_array[*,3] ), max( cmap_array[*,3] )
    print, min( cmap_array[*,4] ), max( cmap_array[*,4] )
    print, min( cmap_array[*,5] ), max( cmap_array[*,5] )
    print, min( cmap_array[*,6] ), max( cmap_array[*,6] )

    ;; define the min and max of x/y axis
    c1 = ( min_color - 0.05 ) < ( min( comp_color ) - 0.45 )
    c2 = ( max_color - 0.05 ) > ( max( comp_color ) + 0.25 )
    ;c1 = c1 > 0.01

    ;; xrange  
    print, '################################################'
    print, 'OLD C1/C2 : ', c1, c2 
    offset1 = ( color_average - c1 ) 
    offset2 = ( c2 - color_average )
    off_ratio = ( offset2 / offset1 )
    off_diff  = ( offset2 - offset1 )
    if ( ( offset1 gt 0.0 ) and ( offset2 ge 0.0 ) ) then begin 
        if ( off_ratio ge 1.2 ) then begin 
            if ( off_diff ge offset1 ) then begin 
                c1 = ( c1 - ( off_diff / 2.0 ) ) 
                c2 = ( c2 - ( off_diff / 2.5 ) )
            endif else begin 
                c1 = ( c1 - off_diff ) 
                c2 = ( c2 - 0.02 ) 
            endelse
        endif
    endif
    print, 'NEW C1/C2 : ', c1, c2
    print, '################################################'
    xrange = [ c1, c2 ]
    
    ;; yrange
    ytitle = textoidl( '\mu_' + red_filter + ' (mag/arcsec^2)' )
    step1 = step_size
    smag_lim2 = where( old_ori_rsma_1 le ( rlim^0.25 ) )
    if ( filter_index1 lt filter_index2 ) then begin 
        min_y = ( min( old_ori_smag_2[ smag_lim2 ] ) - 0.3 )
        max_y  = ( max( old_ori_smag_2[ smag_lim2 ] ) + 1.0 ) < mue_lim2
        max_y = max_y > max( [ 24.1, ( max( comp_mue_2 ) + 0.3 ) ] )
        max_y2 = ( max( old_ori_smag_2[ smag_lim ] ) + ( step_size / 2.0 ) )
        max_y2 = max_y2 > ( max( comp_mue_2 ) + ( step_size / 2.0 ) )
        max_y2 = max_y2 < mue_lim2
        if ( max_y2 gt max_y ) then begin 
            max_y = ( max_y + ( max_y2 - max_y ) + 0.02 )
        endif
        yrange = [ min_y, max_y ]
        xtitle = filter1 + '-' + filter2 + ' Color'
        step0 = ( min_y + 1.0 ) < ( min( comp_mue_2 ) - ( step_size / 2.0 ) )
        step2 = ( max_y2 - ( step_size / 2.0 ) )
        nbin = long( ( step2 - step0 - 0.5 * step1 ) / step1 )
    endif else begin 
        min_y = ( min( old_ori_smag_1[ smag_lim2 ] ) - 0.3 )
        max_y  = ( max( old_ori_smag_1[ smag_lim2 ] ) + 1.0 ) < mue_lim2
        max_y = max_y > max( [ 24.1, ( max( comp_mue_2 ) + 0.3 ) ] )
        max_y2 = ( max( old_ori_smag_1[ smag_lim ] ) + ( step_size / 2.0 ) ) 
        max_y2 = max_y2 > ( max( comp_mue_1 ) + ( step_size / 2.0 ) )
        max_y2 = max_y2 < mue_lim1
        if ( max_y2 gt max_y ) then begin 
            max_y = ( max_y + ( max_y2 - max_y ) + 0.02 )
        endif
        yrange = [ min_y, max_y ]
        xtitle = filter2 + '-' + filter1 + ' Color'
        step0 = ( min_y + 1.0 ) < ( min( comp_mue_1 ) - ( step_size / 2.0 ) )
        step2 = ( max_y2 - ( step_size / 2.0 ) )
        nbin = long( ( step2 - step0 - 0.5 * step1 ) / step1 )
    endelse
    print, 'STEP0    STEP1    STEP2    N_BINS'
    print, step0, step1, step2, nbin

    bin_ceny = fltarr( nbin )  
    bin_cenx_ori = fltarr( nbin ) 
    bin_cenx_mod = fltarr( nbin ) 
    bin_scat = fltarr( nbin )
    bin_scat_2 = make_array( nbin, value=0.0 )
    for m = 0, ( nbin - 1 ), 1 do begin
        print, '###########################################'
        print, 'Bin Number : ' + string( m + 1 )
        bin_ceny[m] = step0 + step1 / 2.0
        print, 'Bin Center : ' + string( bin_ceny[m] )
        bin_edge1 = step0
        print, 'Bin Edge1 : ' + string( bin_edge1 )
        bin_edge2 = step0 + step1 
        print, 'Bin Edge2 : ' + string( bin_edge2 )
        if ( filter_index1 lt filter_index2 ) then begin 
            index1 = where( ( cmap_array[*,2] ge bin_edge1 ) and $ 
                ( cmap_array[*,2] lt bin_edge2 ) )
            index2 = where( ( cmap_array[*,4] ge bin_edge1 ) and $ 
                ( cmap_array[*,4] lt bin_edge2 ) )
            if ( ( index1[0] eq -1 ) or ( index2[0] eq -1 ) ) then begin 
                bin_ceny[m] = 0.0
                bin_c_ori = [ 0.0, 0.0, 0.0, 0.0, 0.0 ]
                bin_c_mod = [ 0,0, 0.0, 0.0, 0.0, 0.0 ]
            endif else begin 
                bin_c_ori = cmap_array[ index1, 5 ]
                bin_c_mod = cmap_array[ index2, 6 ]
            endelse
        endif else begin  
            index1 = where( ( cmap_array[*,1] ge bin_edge1 ) and $ 
                ( cmap_array[*,1] lt bin_edge2 ) )
            index2 = where( ( cmap_array[*,3] ge bin_edge1 ) and $ 
                ( cmap_array[*,3] lt bin_edge2 ) )
            if ( ( index1[0] eq -1 ) or ( index2[0] eq -1 ) ) then begin 
                bin_ceny[m] = 0.0
                bin_c_ori = [ 0.0, 0.0, 0.0, 0.0, 0.0 ]
                bin_c_mod = [ 0,0, 0.0, 0.0, 0.0, 0.0 ]
            endif else begin 
                bin_c_ori = cmap_array[ index1, 5 ]
                bin_c_mod = cmap_array[ index2, 6 ]
            endelse
        endelse

        print, 'Number of usef ori data : ' + string( n_elements( bin_c_ori ) )
        print, 'Number of usef mod data : ' + string( n_elements( bin_c_mod ) )
        nn = n_elements( bin_c_ori )
        resistant_mean, bin_c_ori, 2, aa, bb, cc
        dd = bb * sqrt( ( nn - cc ) - 1 )
        ;dd = robust_sigma( bin_c_ori )
        print, 'Mean ori color ' + string( aa ) + ' +/- ' + string( bb )
        bin_cenx_ori[m] = aa
        bin_scat[m] = ( dd / 2.0 )
        resistant_mean, bin_c_mod, 2, aa, bb, cc
        bin_cenx_mod[m] = aa
        print, 'Mean mod color ' + string( aa ) + ' +/- ' + string( bb )
        step0 = bin_edge2
    endfor
    temp = where( bin_ceny ne 0.0 ) 
    bin_ceny = bin_ceny[ temp ]
    bin_cenx_ori = bin_cenx_ori[ temp ]
    bin_cenx_mod = bin_cenx_mod[ temp ]
    bin_scat = bin_scat[ temp ]
    bin_scat_2 = bin_scat_2[ temp ]
    nbin = n_elements( bin_ceny )

    ;; calculate the photometric error from sky uncertainties 
    bin_perr = fltarr( nbin )
    bin_cerr = fltarr( nbin )
    ;a1 = cmap_array[*,1]
    ;a2 = cmap_array[*,2]
    ;index = where( ( a1 ge 0.0 ) and ( a2 ge 0.0 ) )
    ;s1 = ( min( a1[ index ] ) + 0.1 )
    ;s2 = ( min( a2[ index ] ) + 0.1 )
    print, '###########################'
    for m = 0, ( nbin - 1 ), 1 do begin 
        if ( filter_index1 gt filter_index2 ) then begin 
            m1 = bin_ceny[m]
            m2 = ( m1 + bin_cenx_ori[m] )
        endif else begin 
            m2 = bin_ceny[m] 
            m1 = ( m2 + bin_cenx_ori[m] )
        endelse
        ;m1 = s1 + ( step1 / 2.0 ) 
        ;m2 = s2 + ( step1 / 2.0 ) 
        f1 = ( 10.0^( ( magzpt1 - m1 ) / 2.50 ) * old_expt1 * pix_area ) 
        f2 = ( 10.0^( ( magzpt2 - m2 ) / 2.50 ) * old_expt2 * pix_area )
        e1 = abs( 2.50 * alog10( 1.0 + ( sky_ner_1 / f1 ) ) )
        e2 = abs( 2.50 * alog10( 1.0 + ( sky_ner_2 / f2 ) ) )
        print, '#############################################'
        print, 'Mu ' + filter1 + ' ' + filter2, m1, m2 
        print, 'Sky Err ' + filter1 + ' ' + filter2, sky_ner_1, sky_ner_2
        print, 'Err 1 ' + filter1 + ' ' + filter2, e1, e2
        e1 = sqrt( e1^2.0 + ( photoerr_array[ filter_index1 ] / 2.0 )^2.0 )
        e2 = sqrt( e2^2.0 + ( photoerr_array[ filter_index2 ] / 2.0 )^2.0 )
        print, 'Err 2 ' + filter1 + ' ' + filter2, e1, e2
        e3 = sqrt( e1^2.0 + e2^2.0 )
        print, filter1 + ' ' + filter2 + ' Color Err : ', e3
        if ( filter_index1 gt filter_index2 ) then begin 
            bin_perr[m] = e1 
        endif else begin 
            bin_perr[m] = e2 
        endelse
        bin_cerr[m] = e3 
        ;s1 = s1 + step1 
        ;s2 = s2 + step1
    endfor

    ;;
    comp_cerr = fltarr( ncomp )
    for l = 0, ( ncomp_1 - 1 ), 1 do begin 
        if ( filter_index1 lt filter_index2 ) then begin
            mu0 = comp_mue_2[l] 
            mu1 = mu0 - 0.35 
            mu2 = mu0 + 0.35 
            pix_useful = where( ( cmap_array[*,2] ge mu1 ) and $ 
                ( cmap_array[*,2] lt mu2 ) )
            if ( ( pix_useful[0] ne -1 ) and ( n_elements( pix_useful ) ge 3 ) ) $
                then begin 
                c = cmap_array[ pix_useful, 5 ]
                nn = n_elements( c )
                resistant_mean, c, 5, aa, bb, cc
                dd = bb * sqrt( ( nn - cc ) - 1 )
                comp_cerr[l] = ( dd / 2.0 )
            endif else begin 
                comp_cerr[l] = 0 
            endelse 
        endif else begin 
            mu0 = comp_mue_1[l] 
            mu1 = mu0 - 0.35 
            mu2 = mu0 + 0.35 
            pix_useful = where( ( cmap_array[*,1] ge mu1 ) and $ 
                ( cmap_array[*,1] lt mu2 ) )
            if ( ( pix_useful[0] ne -1 ) and ( n_elements( pix_useful ) ge 3 ) ) $
                then begin 
                c = cmap_array[ pix_useful, 5 ]
                nn = n_elements( c )
                resistant_mean, c, 5, aa, bb, cc
                dd = bb * sqrt( ( nn - cc ) - 1 )
                comp_cerr[l] = ( dd / 2.0 )
            endif else begin 
                comp_cerr[l] = 0 
            endelse 
        endelse 
    endfor

    ;; 
    max_bcerr = max( bin_cerr ) 
    if ( max_bcerr > 0.25 ) then begin 
        xrange[1] = ( xrange[1] + ( max_bcerr / 2.0 ) + 0.02 ) 
    endif

    ;; tick size 
    if ( ( xrange[1] - xrange[0] ) ge 0.5 ) then begin 
        if ( ( xrange[1] - xrange[0] ) le 0.9 ) then begin 
            xtick = 0.2
        endif else begin 
            if ( ( xrange[1] - xrange[0] ) ge 1.2 ) then begin 
                xtick = 0.4 
            endif else begin 
                xtick = 0.3
            endelse
        endelse
    endif else begin 
        xtick = 0.1
    endelse
    print, 'XTICK : ' + string( xtick )

    ;; correct for the situation where the first or last x-tick is too close to 
    ;; the edge of the plot 
    dist_tick_1 = ( xrange[1] - ( indgen( 1000 ) * xtick ) )
    change_1 = where( dist_tick_1 le 0.0 ) 
    edge_1 = ( xtick * ( change_1[0] - 1 ) )
    dist_1 = ( xrange[1] - edge_1 )
    if ( edge_1 le ( xtick / 4.0 + 0.01 ) ) then begin 
        xrange[1] = ( edge_1 - 0.001 ) 
    endif 

;    dist_tick_2 = ( xrange[0] - ( indgen( 1000 ) * xtick ) )
;    change_2 = where( dist_tick_2 le 0.0 ) 
;    edge_2 = ( xtick * ( change_2[0] ) )
;    dist_2 = ( edge_2 - xrange[0] )
;    if ( edge_2 le ( xtick / 4.0 - 0.01 ) ) then begin 
;        xrange[0] = ( edge_2 + 0.001 ) 
;    endif 
    
    ;; make a 2d-hist 
    bin1 = 0.01 
    bin2 = 0.2
    if ( filter_index1 gt filter_index2 ) then begin 
        cmd_hist = hist_2d( cmap_array[*,5], cmap_array[*,1], bin1=bin1, $
            bin2=bin2, min1=xrange[0], min2=yrange[0], max1=xrange[1], $
            ;max2=( yrange[1] + 2 * bin2 ) )
            max2=yrange[1] )
        mod_hist = hist_2d( cmap_array[*,6], cmap_array[*,3], bin1=bin1, $
            bin2=bin2, min1=xrange[0], min2=yrange[0], max1=xrange[1], $
            ;max2=( yrange[1] + 2 * bin2 ) )
            max2=yrange[1] )
    endif else begin 
        cmd_hist = hist_2d( cmap_array[*,5], cmap_array[*,2], bin1=bin1, $
            bin2=bin2, min1=xrange[0], min2=yrange[0], max1=xrange[1], $
            ;max2=( yrange[1] + 2 * bin2 ) )
            max2=yrange[1] )
        mod_hist = hist_2d( cmap_array[*,6], cmap_array[*,4], bin1=bin1, $
            bin2=bin2, min1=xrange[0], min2=yrange[0], max1=xrange[1], $
            ;max2=( yrange[1] + 2 * bin2 ) )
            max2=yrange[1] )
    endelse
    max_ori_hist = max( cmd_hist )
    max_mod_hist = min( mod_hist )

    ;; save the cmd into two ascii files 
    cmd_ori_file = galaxy + '_' + filter1 + '_' + filter2 + '_' + ncomp + $
        '_' + type +'_ori_cmd.fits' 
    cmd_mod_file = galaxy + '_' + filter1 + '_' + filter2 + '_' + ncomp + $
        '_' + type +'_mod_cmd.fits' 
    dimen = size( cmd_hist, /dimension ) 
    mwrfits, cmd_hist, cmd_ori_file, /create 
    mwrfits, mod_hist, cmd_mod_file, /create 
    temp = mrdfits( cmd_ori_file, 0, head )
    sxaddpar, head, 'MIN_X', xrange[0]
    sxaddpar, head, 'MAX_X', xrange[1]
    sxaddpar, head, 'MIN_Y', yrange[0]
    sxaddpar, head, 'MAX_Y', yrange[1]
    sxaddpar, head, 'BIN_1', bin1
    sxaddpar, head, 'BIN_2', bin2
    mwrfits, cmd_hist, cmd_ori_file, head, /create 
    mwrfits, mod_hist, cmd_mod_file, head, /create 

    ;; plot a diagnostic fiture   
    diag_plot = galaxy + '_' + filter1 + '_' + filter2 + '_' + ncomp + $
        '_' + type + '_smag_color.ps'
    mydevice = !D.NAME
    set_plot, 'PS'
    device, filename=diag_plot, font_size=8, /encapsul, $
      /color, set_font='HELVETICA BOLD', /tt_font, xsize=24, ysize=16
    ;; original image
    cgPlot, cmap_array[*,1], cmap_array[*,5], xstyle=1, ystyle=1, $
        xrange=xrange, yrange=yrange, xtitle=xtitle, ytitle=ytitle, $
        xtickinterval=xtick, ytickinterval=2.0, charsize=2.2, charthick=3.8, $
        xthick=4.2, ythick=4.2, /noerase, /nodata, position=position_4
    sky, cmd_hist, bg_cmd, sig_cmd, /quiet
    min_cmd = bg_cmd + 0.5 * sig_cmd 
    max_cmd = max( cmd_hist ) - sig_cmd
    factor_cmd = asinh( ( max_cmd - min_cmd ) / sig_cmd )
    scale_cmd = asinh( ( cmd_hist - min_cmd ) / sig_cmd ) / factor_cmd
    bscale_cmd = 254 - bytscl( scale_cmd ) 
    imgunder, bscale_cmd

    print, '############################################################'
    print, bin_ceny
    print, bin_cenx_ori
    print, bin_cenx_mod 
    print, bin_scat
    print, bin_perr 
    print, bin_cerr
    print, '############################################################'
    print, comp_cerr 
    print, '############################################################'

    cgPlot, [color_average,color_average], !Y.Crange, psym=0, linestyle=5, $
        thick=3.5, /overplot, color=cgColor( 'CYAN', !D.Table_Size )

    cgPlot, bin_cenx_ori, bin_ceny, psym=0, linestyle=2, thick=4.0, $
        color=cgColor( 'ORG5', !D.Table_Size ), /overplot
    cgPlot, bin_cenx_ori, bin_ceny, psym=4, symsize=1.15, thick=4.0, $
        color=cgColor( 'ORG5', !D.Table_Size ), /overplot
    oploterror, bin_cenx_ori, bin_ceny, bin_scat, bin_scat_2, psym=3, $ 
        errcolor=cgColor( 'ORG5', !D.Table_Size ), errthick=3.5
    cgPlot, bin_cenx_mod, bin_ceny, psym=0, linestyle=3, thick=4.0, $
        color=cgColor( 'YELLOW', !D.Table_Size ), /overplot
    cgPlot, bin_cenx_mod, bin_ceny, psym=9, symsize=1.05, thick=4.0, $
        color=cgColor( 'YELLOW', !D.Table_Size ), /overplot
    if ( ncomp_1 eq 3 ) then begin 
        cname_arr = [ 'RED', 'GREEN', 'BLUE' ]
    endif else begin 
        if ( ncomp_1 eq 4 ) then begin 
            cname_arr = [ 'RED', 'GREEN', 'ORANGE', 'BLUE' ] 
        endif else begin 
            cname_arr = [ 'RED', 'GREEN', 'ORANGE', 'BLUE', 'YELLOW', 'BROWN', $
                'PURPLE', 'SALMON', 'GRAY' ]
        endelse
    endelse 

    symbols = [ 46, 14, 16, 15, 17, 18, 19, 20 ]
    for i = 0, ( ncomp_1 - 1 ), 1 do begin 
        if ( filter_index1 lt filter_index2 ) then begin 
            mu = comp_mue_2[i]
        endif else begin 
            mu = comp_mue_1[i] 
        endelse
        c = comp_color[i] 
        e = comp_cerr[i]
        cgPlots, c, mu, psym=symbols[i], symsize=comp_symsize[i], $ 
            color=cgColor( cname_arr[i], !D.Table_Size ), /data
        oploterror, [c], [mu], [e], [0], psym=3, errthick=3.5, $
            errcolor=cgColor( cname_arr[i], !D.Table_Size ), /data
    endfor
    cgPlot, cmap_array[*,1], cmap_array[*,5], xstyle=1, ystyle=1, $
        xrange=xrange, yrange=yrange, $
        xtickinterval=xtick, ytickinterval=2.0, charsize=2.2, charthick=3.8, $
        xthick=4.2, ythick=4.2, /noerase, /nodata, position=position_4
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

    ;; model image 
    cgPlot, cmap_array[*,1], cmap_array[*,5], xstyle=1, ystyle=1, $
        xrange=xrange, yrange=yrange, xtitle=xtitle, ytickformat='(A1)',  $
        xtickinterval=xtick, ytickinterval=2.0, charsize=2.2, charthick=3.8, $
        xthick=3.8, ythick=3.8, /noerase, /nodata, position=position_5
    mod_res = 253 - bytscl( alog10( mod_hist ) )
    sky, mod_hist, bg_mod, sig_mod, /quiet
    min_mod = bg_cmd + 0.5 * sig_cmd 
    max_mod = max( mod_hist ) - sig_cmd
    factor_mod = asinh( ( max_mod - min_mod ) / sig_cmd )
    scale_mod = asinh( ( mod_hist ) / sig_cmd ) / factor_mod
    bscale_mod = 254 - bytscl( scale_mod ) 
    imgunder, bscale_mod
    cgText, 0.510, 0.18, 'MODEL', charsize=2.4, charthick=4.0, $
        color=cgColor( 'BLACK', !D.Table_Size ), /normal
    cgText, 0.515, 0.24, type, charsize=2.1, charthick=3.6, /normal, $ 
        color=cgColor( 'BLACK', !D.Table_Size )
    cgPlots, [ 0.499, 0.553 ], [ 0.290, 0.290 ], psym=0, linestyle=5, $
        thick=3.6, color=cgColor( 'CYAN', !D.Table_Size ), /normal
    temp = '<' + blue_filter + '-' + red_filter + '>'
    cgText, 0.560, 0.285, temp, charsize=1.8, charthick=3.5, $
        alignment=0, /normal, color=cgColor( 'BLACK', !D.Table_Size ) 


    cgPlot, [color_average,color_average], !Y.Crange, psym=0, linestyle=5, $
        thick=3.5, /overplot, color=cgColor( 'CYAN', !D.Table_Size )

    for i = 0, ( ncomp - 1 ), 1 do begin 
        if ( filter_index1 lt filter_index2 ) then begin 
            mu = comp_mue2[i]
        endif else begin 
            mu = comp_mue1[i] 
        endelse
        c = comp_color[i] 
        e = comp_cerr[i]
        cgPlots, c, mu, psym=symbols[i], symsize=comp_symsize[i], $ 
            color=cgColor( cname_arr[i], !D.Table_Size ), /data, thick=3.5
        oploterror, [c], [mu], [e], [0], psym=3, errthick=3.5, $
            errcolor=cgColor( cname_arr[i], !D.Table_Size ), /data
    endfor

    offset = ( ( xrange[1] - xrange[0] ) / 24 ) 
    if ( offset ge 0.06 ) then begin 
        offset = 0.005
    endif 
    max_cerr = max( bin_cerr / 2.0 ) + offset 
    xc = xrange[1] - max_cerr 
    for i = 0, ( nbin - 1 ), 1 do begin 
        ra = ( bin_cerr[i] / 2.0 )
        rb = ( bin_perr[i] / 2.0 )
        xc = xc 
        yc = bin_ceny[i] 
        pos_ang = 0.0 
        tvellipse, ra, rb, xc, yc, pos_ang, /data, $ 
            color=cgColor( 'BLACK', !D.Table_Size ), thick=2.8 
        tvellipse, ra, rb, xc, yc, pos_ang, /data, $ 
            color=cgColor( 'GRAY', !D.Table_Size ), /fill 
    endfor

    x0 = xrange[0] + ( ( xrange[1] - xrange[0] ) / 13 ) 
    y0 = yrange[1] - ( ( yrange[1] - yrange[0] ) / 2.5 )
    if ( filter_index1 gt filter_index2 ) then begin 
        diff_extinc = ( ( extinc_2 / extinc_1 ) - 1.0 ) * 0.25 
    endif else begin   
        diff_extinc = ( ( extinc_1 / extinc_2 ) - 1.0 ) * 0.25
    endelse
    x1 = x0 + diff_extinc
    y1 = y0 + 0.2 
    x2 = ( ( x0 + x1 ) / 2.0 ) * 1.01
    y2 = ( y0 - 0.4  )
    cgArrow, x0, y0, x1, y1,  /data, hthick=5.0, thick=4.0, $
        color=cgColor( 'RED', !D.Table_Size ), hsize=5.5
    ext_string = textoidl( 'A_' + red_filter + '=0.2' )
    cgText, x2, y2, ext_string,  charsize=1.7, charthick=3.4, $
        color=cgColor( 'BLACK', !D.Table_Size ), /data, alignment=0.5

    cgPlot, cmap_array[*,1], cmap_array[*,5], xstyle=1, ystyle=1, $
        xrange=xrange, yrange=yrange, ytickformat='(A1)', $
        xtickinterval=xtick, ytickinterval=2.0, charsize=2.2, charthick=3.8, $
        xthick=3.8, ythick=3.8, /noerase, /nodata, position=position_5
    
    xrange = [ -0.199, 0.199 ]
    yrange = yrange
    ori_serr_1 = ori_serr_1[ 1:( n_elements(old_ori_rsma_1) - 1 ) ]
    ori_serr_2 = ori_serr_2[ 1:( n_elements(old_ori_rsma_1) - 1 ) ]
    old_ori_rsma_1 = old_ori_rsma_1[ 1:( n_elements(old_ori_rsma_1) - 1 ) ] 
    old_ori_rsma_2 = old_ori_rsma_2[ 1:( n_elements(old_ori_rsma_1) - 1 ) ] 
    old_ori_smag_1 = old_ori_smag_1[ 1:( n_elements(old_ori_rsma_1) - 1 ) ]
    old_ori_smag_2 = old_ori_smag_2[ 1:( n_elements(old_ori_rsma_1) - 1 ) ]
    old_res_1 = old_res_1[ 1:( n_elements(old_ori_rsma_1) - 1 ) ]
    old_res_2 = old_res_2[ 1:( n_elements(old_ori_rsma_1) - 1 ) ]

    ;; plot the residual profile 
    cgPlot, old_res_1, old_ori_smag_1, xstyle=1, ystyle=1, xrange=xrange, $
        yrange=yrange, /nodata, position=position_6, $
        /noerase, xtickv=[-0.1,0.1], charsize=2.2, $
        charthick=3.8, xthick=4.2, xminor=-1, ytickformat='(A1)', $
        ythick=4.2, xtickformat='(A1)', ytickinterval=2.0

    p0 = ( position_6[2] - position_6[0] ) / 2.0 + position_6[0] 
    cgText, p0, 0.058, 'Residual', charsize=2.2, charthick=3.8, /normal, $
        alignment = 0.5
    cgText, ( p0 - 0.035 ), 0.10, '-0.1', charsize=2.2, charthick=3.8, /normal, $
        alignment = 0.5 
    cgText, ( p0 + 0.035 ), 0.10, '0.1',  charsize=2.2, charthick=3.8, /normal, $
        alignment = 0.5 

    cgPlot, [0.00,0.00], !Y.Crange, psym=0, linestyle=0, thick=3.5, /overplot, $
        color=cgColor( 'BLK2', !D.Table_Size )

    if ( filter_index1 gt filter_index2 ) then begin 
        cgPlot, abs( ori_serr_1 ), old_ori_smag_1, psym=0, linestyle=4, $
            thick=3.0, /overplot, color=cgColor( 'Dark GRAY', !D.Table_Size )
        cgPlot, -1.0 * abs( ori_serr_1 ), old_ori_smag_1, psym=0, linestyle=4, $ 
            thick=3.0, /overplot, color=cgColor( 'Dark GRAY', !D.Table_Size )
        cgPlot, abs( ori_serr_2 ), old_ori_smag_2, psym=0, linestyle=5, $
            thick=3.0, /overplot, color=cgColor( 'GRAY', !D.Table_Size )
        cgPlot, -1.0 * abs( ori_serr_2 ), old_ori_smag_1, psym=0, linestyle=5, $ 
            thick=3.0, /overplot, color=cgColor( 'GRAY', !D.Table_Size )
        cgPlot, old_res_1, old_ori_smag_1, psym=0, linestyle=2, thick=4.5, $ 
            color=cgColor( 'RED', !D.Table_Size ), /overplot 
        cgPlot, old_res_2, old_ori_smag_1, psym=0, linestyle=3, thick=4.5, $ 
            color=cgColor( 'BLUE', !D.Table_Size ), /overplot
    endif else begin 
        cgPlot, abs( ori_serr_1 ), old_ori_smag_2, psym=0, linestyle=4, $
            thick=3.0, /overplot, color=cgColor( 'Dark GRAY', !D.Table_Size )
        cgPlot, -1.0 * abs( ori_serr_1 ), old_ori_smag_2, psym=0, linestyle=4, $ 
            thick=3.0, /overplot, color=cgColor( 'Dark GRAY', !D.Table_Size )
        cgPlot, abs( ori_serr_2 ), old_ori_smag_2, psym=0, linestyle=5, $
            thick=3.0, /overplot, color=cgColor( 'GRAY', !D.Table_Size )
        cgPlot, -1.0 * abs( ori_serr_2 ), old_ori_smag_2, psym=0, linestyle=5, $ 
            thick=3.0, /overplot, color=cgColor( 'GRAY', !D.Table_Size )
        cgPlot, old_res_1, old_ori_smag_2, psym=0, linestyle=3, thick=4.5, $ 
            color=cgColor( 'BLUE', !D.Table_Size ), /overplot 
        cgPlot, old_res_2, old_ori_smag_2, psym=0, linestyle=2, thick=4.5, $ 
            color=cgColor( 'RED', !D.Table_Size ), /overplot
    endelse

    device, /close
    set_plot, mydevice
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
