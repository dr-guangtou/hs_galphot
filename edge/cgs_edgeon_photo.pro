pro cgs_edgeon_photo, model 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;Do photometry on edge-on galaxy using the major/minor axis cut 
;;written by SHUANG, 06.15.2012 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
on_error, 2
compile_opt idl2

if N_params() lt 1 then begin 
    print,  'Syntax - cgs_edgeon_photo, model ' 
    return
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pix = 0.259 
pix_area = ( pix^2.0 )
bad_thre  = 0.35
;; for boxes 
;; major axis 
ma_height = 12 
ma_cenbox = 1 
ma_width  = 2 
ma_offset = 1 
ma_step1  = 6
ma_step2  = 3
;; minor axis 
mi_height = 2 
mi_cenbox = 1 
mi_width  = 10 
mi_offset = 1
mi_step1  = 2
mi_step2  = 6
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
location_setting, setting
;;; define positions of files
    galfit = setting.galfit
  cat_head = setting.header
  datafile = setting.header
   fileloc = setting.workplace
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; find the model 
model = strcompress( model, /remove_all ) 
if NOT file_test( model ) then begin 
    message, 'Can not find the model ! '
endif
name_string = strsplit( model, '._', /extract ) 
dim = n_elements( name_string )
if ( dim eq 6 ) then begin
    galaxy_name = name_string[0]
    band_name = name_string[1]
    model_number = name_string[3]
endif else begin
    galaxy_name = name_string[0] + '_' + name_string[1]
    band_name = name_string[2]
    model_number = name_string[4]
endelse
comp_number = strmid( model_number, 0, 1 )
ncomp = long( comp_number )
temp = strsplit( model, '.', /extract ) 
model_string = strcompress( temp[0], /remove_all ) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; find galaxy information 
get_head_info, galaxy_name, band_name, header, exist
if ( exist eq 0 ) then begin 
    message, 'Can not find the header information for filter ' + $
        band_name  + '!!'
endif
;; get the useful information
std_name = header.std_name
r50 = float( header.r50 )
r80 = float( header.r80 )
mlim = float( header.mlim )
zpt_gsc = float( header.zpt_gsc )
zpt_lan = float( header.zpt_lan )
ell = float( header.ell_e )
pa = float( header.ell_pa )
exptime = float( header.old_expt )
if ( zpt_lan gt 15.0 ) then begin 
    magzpt = zpt_lan
endif else begin 
    magzpt = zpt_gsc
endelse
;; find the original image 
cor_file = fileloc + galaxy_name + '/' + band_name + '/' + galaxy_name + '_' + $ 
    band_name + '_cor.fit' 
if NOT file_test( cor_file ) then begin 
    message, 'Can not find the corrected image : ' + cor_file 
endif else begin 
    img_cor = mrdfits( cor_file, 0, head_cor ) 
    sky_new = fxpar( head_cor, 'SKY_NV' ) 
    sky_err = fxpar( head_cor, 'SKY_NE' )
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; find the read-in file and the note file  
temp = strsplit( model, '.', /extract ) 
mod_string = strcompress( temp[0], /remove_all ) 
read_file = mod_string + '.read' 
note_file = mod_string + '.note' 
;; read in 
if NOT file_test( read_file ) then begin 
    message, 'Can not find the read-in file : ' + read_file 
endif else begin 

    mod_nosky_file = mod_string + '_nosky_mof.fit' 
    if NOT file_test( mod_nosky_file ) then begin 
        model_gen_nosky, read_file, /run 
    endif
    if NOT file_test( mod_nosky_file ) then begin 
        message, 'Can not find the nosky model image : ' + mod_nosky_file 
    endif else begin 
        img_mod = mrdfits( mod_nosky_file, 0 ) 
    endelse

    ori_nosky_file = mod_string + '_nosky_ori.fit'
    ori_rot_file = mod_string + '_ori_rot.fit'
    if NOT file_test( ori_nosky_file ) then begin 
        model_sky_subtract, read_file
    endif 
    if NOT file_test( ori_nosky_file ) then begin 
        message, 'Can not find the nosky original image : ' + ori_nosky_file 
    endif else begin 
        img_ori = mrdfits( ori_nosky_file, 0 ) 
    endelse

endelse
;; note file 
if NOT file_test( note_file ) then begin 
    message, 'Can not find the note file : ' + note_file 
endif else begin 
    read_galaxy_note, note_file, model_summary=model_sum, comp_summary=comp_sum 
    xcen = long( comp_sum[0].cenx )
    ycen = long( comp_sum[0].ceny )
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; find the mask image 
msk_file = fileloc + galaxy_name + '/' + band_name + '/' + galaxy_name + $ 
    '_' + band_name + '_mm.fits' 
if NOT file_test( msk_file ) then begin 
    message, 'Can not find the mask image : ' + msk_file 
endif else begin 
    img_msk = mrdfits( msk_file, 0 ) 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print, '##########################################################'
print, 'Model     :  ' + mod_string 
print, 'Ncomp     :  ', ncomp
print, 'XCEN/YCEN :  ', xcen, ycen 
print, 'R80       :  ', r80 
print, 'PA        :  ', pa 
;; rotate galaxy image 
cgs_galrot, img_ori, xcen, ycen, r80, pa, ori_rot, mask_name=img_msk, $
    /savefits, fits_name=ori_rot_file 
cgs_galrot, img_mod, xcen, ycen, r80, pa, mod_rot
msk_rot_file = mod_string + '_ori_rot_mask.fit'
if NOT file_test( msk_rot_file ) then begin 
    message, 'Can not find the rotated mask image : ' + msk_rot_file 
endif else begin 
    msk_rot = mrdfits( msk_rot_file, 0 ) 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
temp = size( img_ori, /dimension ) 
naxis1a = temp[0] 
naxis2a = temp[1] 
temp = size( ori_rot, /dimension ) 
naxis1 = temp[0] 
naxis2 = temp[1] 
new_xcen = long( naxis1 / 2.0 ) 
new_ycen = long( naxis2 / 2.0 - 2 ) 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; generate model image for each component, read in and rotate them  
sub_file = mod_string + '_subcomps.fits' 
if NOT file_test( sub_file ) then begin 
    spawn, galfit + ' -o3 ' + read_file
    if NOT file_test( 'subcomps.fits' ) then begin 
        message, 'Can not find the subcomps.fits !' 
    endif else begin 
        spawn, 'mv subcomps.fits ' + sub_file
    endelse
endif 
;;
print, 'There are ' + comp_number + ' components !!!'
img_comps = { img:fltarr( naxis1a, naxis2a ) }
img_comps = replicate( img_comps, ncomp ) 
for k = 0, ( ncomp - 1 ), 1 do begin 
    img_comps[k].img = mrdfits( sub_file, k+1 )
endfor
comps_rot = img_comps
for k = 0, ( ncomp - 1 ), 1 do begin 
    cgs_galrot, img_comps[k].img, xcen, ycen, r80, pa, rot_temp 
    comps_rot[k].img = rot_temp  
endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; release some memory 
img_ori = 0 
img_msk = 0 
img_mod = 0 
img_cor = 0
img_comps = 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define parameters for major axis cut 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ma_nbox = long( ( new_xcen - ( ma_width ) ) / ma_width ) 
ma_struc = { ma_x0p:lonarr( ma_nbox ), ma_y0p:lonarr( ma_nbox ), $ 
    ma_x0n:lonarr( ma_nbox ), ma_y0n:lonarr( ma_nbox ), $  
    ma_x1p:lonarr( ma_nbox ), ma_y1p:lonarr( ma_nbox ), $  
    ma_x1n:lonarr( ma_nbox ), ma_y1n:lonarr( ma_nbox ), $  
    ma_x2p:lonarr( ma_nbox ), ma_y2p:lonarr( ma_nbox ), $  
    ma_x2n:lonarr( ma_nbox ), ma_y2n:lonarr( ma_nbox ), $ 
    ma_x3p:lonarr( ma_nbox ), ma_y3p:lonarr( ma_nbox ), $  
    ma_x3n:lonarr( ma_nbox ), ma_y3n:lonarr( ma_nbox ), $
    ma_pix_p:lonarr( ma_nbox ), ma_use_p:lonarr( ma_nbox ), $ 
    ma_dis_p:fltarr( ma_nbox ), ma_fra_p:lonarr( ma_nbox ), $ 
    ma_pix_n:lonarr( ma_nbox ), ma_use_n:lonarr( ma_nbox ), $ 
    ma_dis_n:fltarr( ma_nbox ), ma_fra_n:lonarr( ma_nbox ), $
    ma_mea_ori_p:fltarr( ma_nbox ), ma_mea_ori_n:fltarr( ma_nbox ), $
    ma_tot_ori_p:fltarr( ma_nbox ), ma_tot_ori_n:fltarr( ma_nbox ), $
    ma_sca_ori_p:fltarr( ma_nbox ), ma_sca_ori_n:fltarr( ma_nbox ), $
    ma_med_ori_p:fltarr( ma_nbox ), ma_med_ori_n:fltarr( ma_nbox ), $
    ma_tmu_ori_p:fltarr( ma_nbox ), ma_tmu_ori_n:fltarr( ma_nbox ), $
    ma_amu_ori_p:fltarr( ma_nbox ), ma_amu_ori_n:fltarr( ma_nbox ), $
    ma_sta_ori_p:fltarr( ma_nbox ), ma_sta_ori_n:fltarr( ma_nbox ), $
    ma_mea_mod_p:fltarr( ma_nbox ), ma_mea_mod_n:fltarr( ma_nbox ), $
    ma_tot_mod_p:fltarr( ma_nbox ), ma_tot_mod_n:fltarr( ma_nbox ), $
    ma_sca_mod_p:fltarr( ma_nbox ), ma_sca_mod_n:fltarr( ma_nbox ), $
    ma_med_mod_p:fltarr( ma_nbox ), ma_med_mod_n:fltarr( ma_nbox ), $
    ma_tmu_mod_p:fltarr( ma_nbox ), ma_tmu_mod_n:fltarr( ma_nbox ), $
    ma_amu_mod_p:fltarr( ma_nbox ), ma_amu_mod_n:fltarr( ma_nbox ), $
    ma_res_p:fltarr( ma_nbox ), ma_res_n:fltarr( ma_nbox ), $ 
    ma_pix:lonarr( ma_nbox ), ma_use:lonarr( ma_nbox ), $
    ma_dis:lonarr( ma_nbox ), ma_fra:lonarr( ma_nbox ), $
    ma_mea_ori:fltarr( ma_nbox ), ma_mea_mod:fltarr( ma_nbox ), $ 
    ma_tot_ori:fltarr( ma_nbox ), ma_tot_mod:fltarr( ma_nbox ), $ 
    ma_med_ori:fltarr( ma_nbox ), ma_med_mod:fltarr( ma_nbox ), $ 
    ma_sca_ori:fltarr( ma_nbox ), ma_sca_mod:fltarr( ma_nbox ), $ 
    ma_tmu_ori:fltarr( ma_nbox ), ma_tmu_mod:fltarr( ma_nbox ), $ 
    ma_amu_ori:fltarr( ma_nbox ), ma_amu_mod:fltarr( ma_nbox ), $ 
    ma_sta_ori:fltarr( ma_nbox ), ma_sta_mod:fltarr( ma_nbox ), $ 
    ma_res: fltarr( ma_nbox ), ma_mod_log_flux:fltarr( ma_nbox ), $
    ma_reg_p:strarr( ma_nbox ), ma_reg_n:strarr( ma_nbox ) }
ma_comp = { ma_mea_p:fltarr( ma_nbox ), ma_tot_p:fltarr( ma_nbox ), $ 
    ma_med_p:fltarr( ma_nbox ), ma_sca_p:fltarr( ma_nbox ), $
    ma_tmu_p:fltarr( ma_nbox ), ma_amu_p:fltarr( ma_nbox ), $
    ma_mea_n:fltarr( ma_nbox ), ma_tot_n:fltarr( ma_nbox ), $ 
    ma_med_n:fltarr( ma_nbox ), ma_sca_n:fltarr( ma_nbox ), $
    ma_tmu_n:fltarr( ma_nbox ), ma_amu_n:fltarr( ma_nbox ), $
    ma_mea:fltarr( ma_nbox ), ma_tot:fltarr( ma_nbox ), $ 
    ma_med:fltarr( ma_nbox ), ma_sca:fltarr( ma_nbox ), $
    ma_tmu:fltarr( ma_nbox ), ma_amu:fltarr( ma_nbox )  }
ma_comp = replicate( ma_comp, ncomp )
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
for i = 0, ( ma_nbox - 1 ), 1 do begin  

    h = ma_height 
    c = ma_cenbox 
    w = ma_width 
    o = ma_offset
   s1 = ma_step1 
   s2 = ma_step2 

    print, '#########################################################'
    print, '#########' + string( i ) + '         ##########'
    print, '#########################################################'

    if ( i eq 0 ) then begin
        frac = 0.0 
        ;; central box 
        ma_struc.ma_x0p[i] = long( new_xcen - c / 2 )
        ma_struc.ma_y0p[i] = long( new_ycen - h / 2 ) + o
        ma_struc.ma_x1p[i] = long( new_xcen + c / 2 ) 
        ma_struc.ma_y1p[i] = long( new_ycen - h / 2 ) + o
        ma_struc.ma_x2p[i] = long( new_xcen + c / 2 ) 
        ma_struc.ma_y2p[i] = long( new_ycen + h / 2 ) + o
        ma_struc.ma_x3p[i] = long( new_xcen - c / 2 )
        ma_struc.ma_y3p[i] = long( new_ycen + h / 2 ) + o

        ma_struc.ma_x0n[i] = long( new_xcen - c / 2 )
        ma_struc.ma_y0n[i] = long( new_ycen - h / 2 ) + o
        ma_struc.ma_x1n[i] = long( new_xcen + c / 2 ) 
        ma_struc.ma_y1n[i] = long( new_ycen - h / 2 ) + o
        ma_struc.ma_x2n[i] = long( new_xcen + c / 2 )
        ma_struc.ma_y2n[i] = long( new_ycen + h / 2 ) + o
        ma_struc.ma_x3n[i] = long( new_xcen - c / 2 )
        ma_struc.ma_y3n[i] = long( new_ycen + h / 2 ) + o

        xx = string( ( ma_struc.ma_x0p[i] + ma_struc.ma_x1p[i] ) / 2.0 )
        yy = string( ( ma_struc.ma_y0p[i] + ma_struc.ma_y2p[i] ) / 2.0 )
        xw = string( abs( ma_struc.ma_x1p[i] - ma_struc.ma_x0p[i] ) )
        yw = string( abs( ma_struc.ma_y2p[i] - ma_struc.ma_y0p[i] ) )
        ma_struc.ma_reg_p[i] = 'box( ' + xx + ' , ' + yy + ' , ' + $
            xw + ' , ' + yw + ' )'

        xx = string( ( ma_struc.ma_x0n[i] + ma_struc.ma_x1n[i] ) / 2.0 )
        yy = string( ( ma_struc.ma_y0n[i] + ma_struc.ma_y2n[i] ) / 2.0 )
        xw = string( abs( ma_struc.ma_x1n[i] - ma_struc.ma_x0n[i] ) )
        yw = string( abs( ma_struc.ma_y2n[i] - ma_struc.ma_y0n[i] ) ) 
        ma_struc.ma_reg_n[i] = 'box( ' + xx + ' , ' + yy + ' , ' + $
            xw + ' , ' + yw + ' )'

        new_xcen_2 = ( float( ma_struc.ma_x0p[i] + ma_struc.ma_x1p[i] ) / 2.0 )

        ma_struc.ma_pix[i] = long( ( c + 1 ) * ( h + 1 ) ) 
        ma_struc.ma_pix_p[i] = ma_struc.ma_pix[i]
        ma_struc.ma_pix_n[i] = ma_struc.ma_pix[i] 

        ma_struc.ma_dis[i] = 0 
        ma_struc.ma_dis_p[i] = ma_struc.ma_dis[i] 
        ma_struc.ma_dis_n[i] = ma_struc.ma_dis[i] 

        if ( ma_struc.ma_x0p[i] ne ma_struc.ma_x1p[i] ) then begin 
            ori_box = ori_rot[ ma_struc.ma_x0p[i]:ma_struc.ma_x1p[i], $
                ma_struc.ma_y0p[i]:ma_struc.ma_y2p[i] ]  
            mod_box = mod_rot[ ma_struc.ma_x0p[i]:ma_struc.ma_x1p[i], $
                ma_struc.ma_y0p[i]:ma_struc.ma_y2p[i] ]  
            msk_box = msk_rot[ ma_struc.ma_x0p[i]:ma_struc.ma_x1p[i], $
                ma_struc.ma_y0p[i]:ma_struc.ma_y2p[i] ]  
        endif else begin 
            ori_box = ori_rot[ ma_struc.ma_x0p[i], $
                ma_struc.ma_y0p[i]:ma_struc.ma_y2p[i] ]  
            mod_box = mod_rot[ ma_struc.ma_x0p[i], $
                ma_struc.ma_y0p[i]:ma_struc.ma_y2p[i] ]  
            msk_box = msk_rot[ ma_struc.ma_x0p[i], $
                ma_struc.ma_y0p[i]:ma_struc.ma_y2p[i] ]  
        endelse

        index_use = where( msk_box lt 1 ) 

        if ( index_use[0] eq -1 ) then begin 
            frac = 0.0 
        endif else begin 
            num_use = n_elements( index_use ) 
            frac = float( num_use ) / float( ma_struc.ma_pix[i] ) 
            ma_struc.ma_use[i] = num_use 
            ma_struc.ma_use_p[i] = num_use
            ma_struc.ma_use_n[i] = num_use 

            ;; ori 
            sum_ori = total( ori_box[ index_use ], /double )
            ma_struc.ma_tot_ori[i] = ( sum_ori / num_use ) * $
                ma_struc.ma_pix[i]
            ma_struc.ma_mea_ori[i] = ( sum_ori / num_use ) 
            ma_struc.ma_med_ori[i] = median( ori_box[ index_use ], /double )
            useful_ori = ori_box[ index_use ]
            ma_struc.ma_sca_ori[i] = robust_sigma( useful_ori ) 
            ma_struc.ma_tmu_ori[i] = -2.50 * alog10( ma_struc.ma_mea_ori[i] / $
                ( pix_area * exptime ) ) + magzpt
            ma_struc.ma_amu_ori[i] = -2.50 * alog10( ma_struc.ma_med_ori[i] / $
                ( pix_area * exptime ) ) + magzpt
            ;; mod 
            sum_mod = total( mod_box[ index_use ], /double )
            ma_struc.ma_tot_mod[i] = ( sum_mod / num_use ) * $
                ma_struc.ma_pix[i]
            ma_struc.ma_mea_mod[i] = ( sum_mod / num_use ) 
            ma_struc.ma_med_mod[i] = median( mod_box[ index_use ], /double )
            useful = mod_box[ index_use ]
            ma_struc.ma_sca_mod[i] = robust_sigma( useful ) 
            ma_struc.ma_tmu_mod[i] = -2.50 * alog10( ma_struc.ma_mea_mod[i] / $
                ( pix_area * exptime ) ) + magzpt
            ma_struc.ma_amu_mod[i] = -2.50 * alog10( ma_struc.ma_med_mod[i] / $
                ( pix_area * exptime ) ) + magzpt
            ;ma_struc.ma_res[i] = ( ma_struc.ma_tmu_ori[i] - $
            ;    ma_struc.ma_tmu_mod[i] )
            ma_struc.ma_res[i] = ( ma_struc.ma_amu_ori[i] - $
                ma_struc.ma_amu_mod[i] )

            ma_struc.ma_mod_log_flux[i] = alog10( total( mod_box ) )

            ma_struc.ma_tot_ori_p[i] = ma_struc.ma_tot_ori[i]  
            ma_struc.ma_tot_ori_n[i] = ma_struc.ma_tot_ori[i] 
            ma_struc.ma_mea_ori_p[i] = ma_struc.ma_mea_ori[i]
            ma_struc.ma_mea_ori_n[i] = ma_struc.ma_mea_ori[i] 
            ma_struc.ma_med_ori_p[i] = ma_struc.ma_med_ori[i]  
            ma_struc.ma_med_ori_n[i] = ma_struc.ma_med_ori[i] 
            ma_struc.ma_sca_ori_p[i] = ma_struc.ma_sca_ori[i] 
            ma_struc.ma_sca_ori_n[i] = ma_struc.ma_sca_ori[i]
            ma_struc.ma_tmu_ori_p[i] = ma_struc.ma_tmu_ori[i] 
            ma_struc.ma_tmu_ori_n[i] = ma_struc.ma_tmu_ori[i] 
            ma_struc.ma_amu_ori_p[i] = ma_struc.ma_amu_ori[i]
            ma_struc.ma_amu_ori_n[i] = ma_struc.ma_amu_ori[i]

            ma_struc.ma_tot_mod_p[i] = ma_struc.ma_tot_mod[i]  
            ma_struc.ma_tot_mod_n[i] = ma_struc.ma_tot_mod[i] 
            ma_struc.ma_mea_mod_p[i] = ma_struc.ma_mea_mod[i]
            ma_struc.ma_mea_mod_n[i] = ma_struc.ma_mea_mod[i] 
            ma_struc.ma_med_mod_p[i] = ma_struc.ma_med_mod[i]  
            ma_struc.ma_med_mod_n[i] = ma_struc.ma_med_mod[i] 
            ma_struc.ma_sca_mod_p[i] = ma_struc.ma_sca_mod[i] 
            ma_struc.ma_sca_mod_n[i] = ma_struc.ma_sca_mod[i]
            ma_struc.ma_tmu_mod_p[i] = ma_struc.ma_tmu_mod[i] 
            ma_struc.ma_tmu_mod_n[i] = ma_struc.ma_tmu_mod[i] 
            ma_struc.ma_amu_mod_p[i] = ma_struc.ma_amu_mod[i]
            ma_struc.ma_amu_mod_n[i] = ma_struc.ma_amu_mod[i]

        endelse

        if ( frac gt bad_thre ) then begin 
            ma_struc.ma_sta_ori_p[i] = 1 
            ma_struc.ma_sta_ori_n[i] = 1 
            ma_struc.ma_sta_ori[i] = 1 
            ma_struc.ma_fra_p[i] = frac 
            ma_struc.ma_fra_n[i] = frac 
            ma_struc.ma_fra[i] = frac
        endif 

        ;; for the components 
        for k = 0, ( ncomp - 1 ), 1 do begin 

            ;; positive side 
            comps_box_p = $
                comps_rot[k].img[ ma_struc.ma_x0p[i]:ma_struc.ma_x1p[i], $
                ma_struc.ma_y0p[i]:ma_struc.ma_y2p[i] ] 
            ma_comp[k].ma_tot_p[i] = total( comps_box_p, /double ) 
            ma_comp[k].ma_mea_p[i] = ( ma_comp[k].ma_tot_p[i] / $
                ma_struc.ma_pix_p[i] )
            ma_comp[k].ma_sca_p[i] = robust_sigma( comps_box_p ) 
            ma_comp[k].ma_med_p[i] = median( comps_box_p, /double ) 
            ma_comp[k].ma_tmu_p[i] = -2.50 * alog10( ma_comp[k].ma_mea_p[i] / $ 
                ( exptime * pix_area ) ) + magzpt 
            ma_comp[k].ma_amu_p[i] = -2.50 * alog10( ma_comp[k].ma_med_p[i] / $ 
                ( exptime * pix_area ) ) + magzpt 

            ;; negative side 
            comps_box_n = $
                comps_rot[k].img[ ma_struc.ma_x0n[i]:ma_struc.ma_x1n[i], $
                ma_struc.ma_y0n[i]:ma_struc.ma_y2n[i] ] 
            ma_comp[k].ma_tot_n[i] = total( comps_box_n, /double ) 
            ma_comp[k].ma_mea_n[i] = ( ma_comp[k].ma_tot_n[i] / $
                ma_struc.ma_pix_n[i] )
            ma_comp[k].ma_sca_n[i] = robust_sigma( comps_box_n ) 
            ma_comp[k].ma_med_n[i] = median( comps_box_n, /double ) 
            ma_comp[k].ma_tmu_n[i] = -2.50 * alog10( ma_comp[k].ma_mea_n[i] / $ 
                ( exptime * pix_area ) ) + magzpt 
            ma_comp[k].ma_amu_n[i] = -2.50 * alog10( ma_comp[k].ma_med_n[i] / $ 
                ( exptime * pix_area ) ) + magzpt 

            index_p = where( comps_box_p gt -99999 ) 
            index_n = where( comps_box_n gt -99999 )
            comps_pix_p = comps_box_p[ index_p ]
            comps_pix_n = comps_box_n[ index_n ]
            ;;average 
            ma_comp[k].ma_tot[i] = ( ma_comp[k].ma_tot_p[i] + $
                ma_comp[k].ma_tot_n[i] ) 
            ma_comp[k].ma_mea[i] = ( ma_comp[k].ma_tot[i] / ( $ 
                ma_struc.ma_pix_p[i] + ma_struc.ma_pix_n[i] ) )
            ma_comp[k].ma_med[i] = median( [ comps_pix_p, comps_pix_n ], $
                /double )
            ma_comp[k].ma_sca[i] = robust_sigma( [ comps_pix_p, comps_pix_n ] ) 
            ma_comp[k].ma_tmu[i] = -2.50 * alog10( ma_comp[k].ma_mea[i] / $ 
                ( pix_area * exptime ) ) + magzpt
            ma_comp[k].ma_amu[i] = -2.50 * alog10( ma_comp[k].ma_med[i] / $ 
                ( pix_area * exptime ) ) + magzpt
        endfor

    ;; for non-central box 
    endif else begin 

        ;; the width of the box should increase stepwise 
        ;; now use R50 and R80 as standard 
        if ( ( ma_struc.ma_dis[ i-1 ] + ma_width ) gt ( r50 * 0.6 ) ) then begin 
            h = h + s1 
            w = w + s2 
        endif 
        if ( ( ma_struc.ma_dis[ i-1 ] + ma_width ) gt r50 )  then begin 
            h = h + s1 
            w = w + s2 
        endif 
        if ( ( ma_struc.ma_dis[ i-1 ] + ma_width ) gt $
            ( ( r50 + r80 ) / 2.0 ) ) then begin 
            h = h + s1 
            w = w + s2
        endif
        if ( ( ma_struc.ma_dis[ i-1 ] + ma_width ) gt r80 ) then begin 
            h = h + s1 
            w = w + s2
        endif
        if ( ( ma_struc.ma_dis[ i-1 ] + ma_width ) gt ( r80 * 1.4 ) ) then begin 
            h = h + s1 
            w = w + s2
        endif

        frac_p = 0.0
        frac_n = 0.0
        nround = 0
        flux_step = 0 

        while ( ( flux_step eq 0 ) and ( nround le 5 ) ) do begin 

            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ;; positive side  
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ma_struc.ma_x0p[i] = ( long( ma_struc.ma_x1p[i-1] ) ) < $
                ( naxis1 - 1 ) 
            ma_struc.ma_y0p[i] = ( long( new_ycen - h / 2 ) ) + o 
            ma_struc.ma_x1p[i] = ( long( ma_struc.ma_x1p[i-1] + w ) ) < $
                ( naxis1 - 1 )
            ma_struc.ma_y1p[i] = ( long( new_ycen - h / 2 ) ) + o
            ma_struc.ma_x2p[i] = ( long( ma_struc.ma_x1p[i-1] + w ) ) < $
                ( naxis1 - 1 )
            ma_struc.ma_y2p[i] = ( long( new_ycen + h / 2 ) ) + o
            ma_struc.ma_x3p[i] = ( long( ma_struc.ma_x1p[i-1] ) ) < $
                ( naxis1 - 1 )
            ma_struc.ma_y3p[i] = ( long( new_ycen + h / 2 ) ) + o
            ma_struc.ma_pix_p[i] = ( long( ( w + 1 ) * ( h + 1 ) ) ) 
            ma_struc.ma_dis_p[i] = abs( ( float( ma_struc.ma_x0p[i] + $
                ma_struc.ma_x1p[i] ) / 2.0 ) - new_xcen_2 ) 
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ;; negative side 
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ma_struc.ma_x0n[i] = ( long( ma_struc.ma_x1n[i-1] ) ) > 0 
            ma_struc.ma_y0n[i] = ( long( new_ycen - h / 2 ) ) + o
            ma_struc.ma_x1n[i] = ( long( ma_struc.ma_x1n[i-1] - w ) ) > 0
            ma_struc.ma_y1n[i] = ( long( new_ycen - h / 2 ) ) + o
            ma_struc.ma_x2n[i] = ( long( ma_struc.ma_x1n[i-1] - w ) ) > 0
            ma_struc.ma_y2n[i] = ( long( new_ycen + h / 2 ) ) + o
            ma_struc.ma_x3n[i] = ( long( ma_struc.ma_x1n[i-1] ) ) > 0 
            ma_struc.ma_y3n[i] = ( long( new_ycen + h / 2 ) ) + o 
            ma_struc.ma_pix_n[i] = ( long( ( w + 1 ) * ( h + 1 ) ) ) 
            ma_struc.ma_dis_n[i] = abs( new_xcen_2 - ( $
                float( ma_struc.ma_x0n[i] + ma_struc.ma_x1n[i] ) / 2.0 ) )

            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ;; save the region information 
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            xx = string( ( ma_struc.ma_x0p[i] + ma_struc.ma_x1p[i] ) / 2.0 )
            yy = string( ( ma_struc.ma_y0p[i] + ma_struc.ma_y2p[i] ) / 2.0 )
            xw = string( abs( ma_struc.ma_x1p[i] - ma_struc.ma_x0p[i] ) )
            yw = string( abs( ma_struc.ma_y2p[i] - ma_struc.ma_y0p[i] ) )
            ma_struc.ma_reg_p[i] = 'box( ' + xx + ' , ' + yy + ' , ' + $
                xw + ' , ' + yw + ' )'
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            xx = string( ( ma_struc.ma_x0n[i] + ma_struc.ma_x1n[i] ) / 2.0 )
            yy = string( ( ma_struc.ma_y0n[i] + ma_struc.ma_y2n[i] ) / 2.0 )
            xw = string( abs( ma_struc.ma_x1n[i] - ma_struc.ma_x0n[i] ) )
            yw = string( abs( ma_struc.ma_y2n[i] - ma_struc.ma_y0n[i] ) ) 
            ma_struc.ma_reg_n[i] = 'box( ' + xx + ' , ' + yy + ' , ' + $
                xw + ' , ' + yw + ' )'

            ori_box_p = ori_rot[ ma_struc.ma_x0p[i]:ma_struc.ma_x1p[i], $
                ma_struc.ma_y0p[i]:ma_struc.ma_y2p[i] ]  
            mod_box_p = mod_rot[ ma_struc.ma_x0p[i]:ma_struc.ma_x1p[i], $
                ma_struc.ma_y0p[i]:ma_struc.ma_y2p[i] ]  
            msk_box_p = msk_rot[ ma_struc.ma_x0p[i]:ma_struc.ma_x1p[i], $
                ma_struc.ma_y0p[i]:ma_struc.ma_y2p[i] ]  
            index_use_p = where( msk_box_p lt 1 ) 
            num_use_p = n_elements( index_use_p )

            ori_box_n = ori_rot[ ma_struc.ma_x1n[i]:ma_struc.ma_x0n[i], $
                ma_struc.ma_y0n[i]:ma_struc.ma_y2n[i] ]  
            mod_box_n = mod_rot[ ma_struc.ma_x1n[i]:ma_struc.ma_x0n[i], $
                ma_struc.ma_y0n[i]:ma_struc.ma_y2n[i] ]  
            msk_box_n = msk_rot[ ma_struc.ma_x1n[i]:ma_struc.ma_x0n[i], $
                ma_struc.ma_y0n[i]:ma_struc.ma_y2n[i] ]  
            index_use_n = where( msk_box_n lt 1 ) 
            num_use_n = n_elements( index_use_n ) 

            ma_struc.ma_mod_log_flux[i] = alog10( ( $
                total( mod_box_n, /double ) + $
                total( mod_box_p, /double ) ) / 2.0 )
            if ( ( ma_struc.ma_mod_log_flux[i-1] - $
                ma_struc.ma_mod_log_flux[i] ) gt 0.50 ) then begin 
                flux_step = 0 
            endif else begin 
                flux_step = 1 
            endelse

            ;; %%%%%%%%% "and" v.s "or"
            if ( index_use_p[0] eq -1 ) then begin 
                frac_p = 0.0 
            endif else begin 

                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ;; positive side 
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                frac_p = float( num_use_p ) / float( ma_struc.ma_pix_p[i] )
                ma_struc.ma_fra_p[i] = frac_p
                ma_struc.ma_use_p[i] = num_use_p 

                ;; ori 
                sum_ori = total( ori_box_p[ index_use_p ], /double ) 
                ma_struc.ma_tot_ori_p[i] = ( ( sum_ori / num_use_p ) * $
                    ma_struc.ma_pix_p[i] )
                ma_struc.ma_mea_ori_p[i] = ( sum_ori / num_use_p ) 
                ma_struc.ma_med_ori_p[i] = median( ori_box_p[ index_use_p ], $
                    /double )
                useful_ori_p = ori_box_p[ index_use_p ]
                ma_struc.ma_sca_ori_p[i] = robust_sigma( useful_ori_p ) 
                if ( ma_struc.ma_mea_ori_p[i] lt 0.0 ) then begin 
                    ma_struc.ma_tmu_ori_p[i] = 99.0
                endif else begin 
                    ma_struc.ma_tmu_ori_p[i] = -2.50 * alog10( $
                        ma_struc.ma_mea_ori_p[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse
                if ( ma_struc.ma_med_ori_p[i] lt 0.0 ) then begin 
                    ma_struc.ma_amu_ori_p[i] = 99.0 
                endif else begin 
                    ma_struc.ma_amu_ori_p[i] = -2.50 * alog10( $
                        ma_struc.ma_med_ori_p[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse

                ;; mod 
                sum_mod = total( mod_box_p[ index_use_p ] )
                ma_struc.ma_tot_mod_p[i] = ( ( sum_mod / num_use_p ) * $ 
                    ma_struc.ma_pix_p[i] )
                ma_struc.ma_mea_mod_p[i] = ( sum_mod / num_use_p ) 
                ma_struc.ma_med_mod_p[i] = median( mod_box_p[ index_use_p ], $
                    /double )
                useful_mod_p = mod_box_p[ index_use_p ]
                ma_struc.ma_sca_mod_p[i] = robust_sigma( useful_mod_p ) 
                if ( ma_struc.ma_mea_mod_p[i] lt 0.0 ) then begin 
                    ma_struc.ma_tmu_mod_p[i] = 99.0
                endif else begin 
                    ma_struc.ma_tmu_mod_p[i] = -2.50 * alog10( $
                        ma_struc.ma_mea_mod_p[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse
                if ( ma_struc.ma_med_mod_p[i] lt 0.0 ) then begin 
                    ma_struc.ma_amu_mod_p[i] = 99.0 
                endif else begin 
                    ma_struc.ma_amu_mod_p[i] = -2.50 * alog10( $
                        ma_struc.ma_med_mod_p[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse

                ;; res 
                ;ma_struc.ma_res_p[i] = ( ma_struc.ma_tmu_ori_p[i] - $
                ;    ma_struc.ma_tmu_mod_p[i] )
                ma_struc.ma_res_p[i] = ( ma_struc.ma_amu_ori_p[i] - $
                    ma_struc.ma_amu_mod_p[i] )
            endelse

            if ( index_use_n[0] eq -1 ) then begin 
                frac_n = 0.0 
            endif else begin 

                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ;; negative side 
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                frac_n = float( num_use_n ) / float( ma_struc.ma_pix_n[i] )
                ma_struc.ma_fra_n[i] = frac_n
                ma_struc.ma_use_n[i] = num_use_n 

                ;; ori 
                sum_ori = total( ori_box_n[ index_use_n ], /double ) 
                ma_struc.ma_tot_ori_n[i] = ( ( sum_ori / num_use_n ) * $
                    ma_struc.ma_pix_n[i] )
                ma_struc.ma_mea_ori_n[i] = ( sum_ori / num_use_n ) 
                ma_struc.ma_med_ori_n[i] = median( ori_box_n[ index_use_n ], $
                    /double )
                useful_ori_n = ori_box_n[ index_use_n ]
                ma_struc.ma_sca_ori_n[i] = robust_sigma( useful_ori_n ) 
                if ( ma_struc.ma_mea_ori_n[i] lt 0.0 ) then begin 
                    ma_struc.ma_tmu_ori_n[i] = 99.0
                endif else begin 
                    ma_struc.ma_tmu_ori_n[i] = -2.50 * alog10( $
                        ma_struc.ma_mea_ori_n[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse
                if ( ma_struc.ma_med_ori_n[i] lt 0.0 ) then begin 
                    ma_struc.ma_amu_ori_n[i] = 99.0 
                endif else begin 
                    ma_struc.ma_amu_ori_n[i] = -2.50 * alog10( $
                        ma_struc.ma_med_ori_n[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse

                ;; mod 
                sum_mod = total( mod_box_n[ index_use_n ] )
                ma_struc.ma_tot_mod_n[i] = ( ( sum_mod / num_use_n ) * $ 
                    ma_struc.ma_pix_n[i] )
                ma_struc.ma_mea_mod_n[i] = ( sum_mod / num_use_n ) 
                ma_struc.ma_med_mod_n[i] = median( mod_box_n[ index_use_n ], $
                    /double )
                useful_mod_n = mod_box_n[ index_use_n ]
                ma_struc.ma_sca_mod_n[i] = robust_sigma( useful_mod_n ) 
                if ( ma_struc.ma_mea_mod_n[i] lt 0.0 ) then begin 
                    ma_struc.ma_tmu_mod_n[i] = 99.0
                endif else begin 
                    ma_struc.ma_tmu_mod_n[i] = -2.50 * alog10( $
                        ma_struc.ma_mea_mod_n[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse
                if ( ma_struc.ma_med_mod_n[i] lt 0.0 ) then begin 
                    ma_struc.ma_amu_mod_n[i] = 99.0 
                endif else begin 
                    ma_struc.ma_amu_mod_n[i] = -2.50 * alog10( $
                        ma_struc.ma_med_mod_n[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse

                ;; res 
                ;ma_struc.ma_res_n[i] = ( ma_struc.ma_tmu_ori_n[i] - $
                ;    ma_struc.ma_tmu_mod_n[i] )
                ma_struc.ma_res_n[i] = ( ma_struc.ma_amu_ori_n[i] - $
                    ma_struc.ma_amu_mod_n[i] )
            endelse
            
            ;; status of the box 
            if ( frac_p ge bad_thre ) then begin 
                ma_struc.ma_sta_ori_p[i] = 1 
            endif else begin 
                ma_struc.ma_sta_ori_p[i] = 0 
            endelse

            if ( frac_n ge bad_thre ) then begin 
                ma_struc.ma_sta_ori_n[i] = 1 
            endif else begin 
                ma_struc.ma_sta_ori_n[i] = 0 
            endelse

            ;; if the good pixel fractions are low at both side 
            ;; %%%%%% "and" or "or" !!!!!!!!!!
            if ( ( frac_p lt bad_thre ) and ( frac_n lt bad_thre ) ) then begin 
                flux_step = 0 
            endif
            ;; 
            if ( flux_step eq 0 ) then begin 
                w = w + s2
            endif 

            ;; average     
            ;; four situation 
            ma_struc.ma_dis[i] = ma_struc.ma_dis_p[i]
            ;; 1
            if ( ( frac_p ge bad_thre ) and ( frac_n lt bad_thre ) ) then begin 
                ma_struc.ma_pix[i] = ma_struc.ma_pix_p[i] 
                ma_struc.ma_use[i] = ma_struc.ma_use_p[i] 
                ma_struc.ma_fra[i] = ma_struc.ma_fra_p[i] 
                ma_struc.ma_tot_ori[i] = ma_struc.ma_tot_ori_p[i] 
                ma_struc.ma_mea_ori[i] = ma_struc.ma_mea_ori_p[i]
                ma_struc.ma_med_ori[i] = ma_struc.ma_med_ori_p[i] 
                ma_struc.ma_sca_ori[i] = ma_struc.ma_sca_ori_p[i] 
                ma_struc.ma_tmu_ori[i] = ma_struc.ma_tmu_ori_p[i] 
                ma_struc.ma_amu_ori[i] = ma_struc.ma_tmu_ori_p[i] 
                ma_struc.ma_tot_mod[i] = ma_struc.ma_tot_mod_p[i] 
                ma_struc.ma_mea_mod[i] = ma_struc.ma_mea_mod_p[i]
                ma_struc.ma_med_mod[i] = ma_struc.ma_med_mod_p[i] 
                ma_struc.ma_sca_mod[i] = ma_struc.ma_sca_mod_p[i] 
                ma_struc.ma_tmu_mod[i] = ma_struc.ma_tmu_mod_p[i] 
                ma_struc.ma_amu_mod[i] = ma_struc.ma_tmu_mod_p[i] 
                ma_struc.ma_sta_ori[i] = 1 
            endif
            ;; 2
            if ( ( frac_p lt bad_thre ) and ( frac_n ge bad_thre ) ) then begin 
                ma_struc.ma_pix[i] = ma_struc.ma_pix_n[i] 
                ma_struc.ma_use[i] = ma_struc.ma_use_n[i] 
                ma_struc.ma_fra[i] = ma_struc.ma_fra_n[i] 
                ma_struc.ma_tot_ori[i] = ma_struc.ma_tot_ori_n[i] 
                ma_struc.ma_mea_ori[i] = ma_struc.ma_mea_ori_n[i]
                ma_struc.ma_med_ori[i] = ma_struc.ma_med_ori_n[i] 
                ma_struc.ma_sca_ori[i] = ma_struc.ma_sca_ori_n[i] 
                ma_struc.ma_tmu_ori[i] = ma_struc.ma_tmu_ori_n[i] 
                ma_struc.ma_amu_ori[i] = ma_struc.ma_tmu_ori_n[i] 
                ma_struc.ma_tot_mod[i] = ma_struc.ma_tot_mod_n[i] 
                ma_struc.ma_med_mod[i] = ma_struc.ma_med_mod_n[i] 
                ma_struc.ma_med_mod[i] = ma_struc.ma_med_mod_n[i] 
                ma_struc.ma_sca_mod[i] = ma_struc.ma_sca_mod_n[i] 
                ma_struc.ma_tmu_mod[i] = ma_struc.ma_tmu_mod_n[i] 
                ma_struc.ma_amu_mod[i] = ma_struc.ma_tmu_mod_n[i] 
                ma_struc.ma_sta_ori[i] = 1 
            endif
            ;; 3 
            if ( ( frac_p lt bad_thre ) and ( frac_n lt bad_thre ) ) then begin 
                ma_struc.ma_pix[i] = ma_struc.ma_pix_n[i] 
                ma_struc.ma_use[i] = 0 
                ma_struc.ma_fra[i] = 0.0
                ma_struc.ma_tot_ori[i] = 0.0
                ma_struc.ma_mea_ori[i] = 0.0 
                ma_struc.ma_med_ori[i] = 0.0
                ma_struc.ma_sca_ori[i] = 0.0
                ma_struc.ma_tmu_ori[i] = 99.0
                ma_struc.ma_amu_ori[i] = 99.0
                ma_struc.ma_tot_mod[i] = 0.0
                ma_struc.ma_med_mod[i] = 0.0
                ma_struc.ma_med_mod[i] = 0.0
                ma_struc.ma_sca_mod[i] = 0.0
                ma_struc.ma_tmu_mod[i] = 99.0
                ma_struc.ma_amu_mod[i] = 99.0
                ma_struc.ma_sta_ori[i] = 0 
            endif
            ;; 4 
            if ( ( frac_p ge bad_thre ) and ( frac_n ge bad_thre ) ) then begin 
                ma_struc.ma_pix[i] = ( ma_struc.ma_pix_p[i] + $
                    ma_struc.ma_pix_n[i] )
                ma_struc.ma_use[i] = ( ma_struc.ma_use_p[i] + $
                    ma_struc.ma_use_n[i] )
                ma_struc.ma_fra[i] = ( float( ma_struc.ma_use[i] ) / $ 
                    float( ma_struc.ma_pix[i] ) )

                ;; ori 
                ma_struc.ma_tot_ori[i] = ( ma_struc.ma_tot_ori_p[i] + $
                    ma_struc.ma_tot_ori_n[i] ) 
                useful_ori = [ useful_ori_p, useful_ori_n ] 
                sum_ori = total( useful_ori, /double ) 
                ma_struc.ma_mea_ori[i] = ( sum_ori / ma_struc.ma_use[i] ) 
                ma_struc.ma_med_ori[i] = median( useful_ori, /double ) 
                ma_struc.ma_sca_ori[i] = robust_sigma( useful_ori )
                ma_struc.ma_tmu_ori[i] = -2.50 * alog10( $
                    ma_struc.ma_mea_ori[i] / ( pix_area * exptime ) ) + magzpt  
                ma_struc.ma_amu_ori[i] = -2.50 * alog10( $
                    ma_struc.ma_med_ori[i] / ( pix_area * exptime ) ) + magzpt

                ;; mod 
                ma_struc.ma_tot_mod[i] = ( ma_struc.ma_tot_mod_p[i] + $
                    ma_struc.ma_tot_mod_n[i] ) 
                useful_mod = [ useful_mod_p, useful_mod_n ] 
                sum_mod = total( useful_mod, /double ) 
                ma_struc.ma_mea_mod[i] = ( sum_mod / ma_struc.ma_use[i] ) 
                ma_struc.ma_med_mod[i] = median( useful_mod, /double ) 
                ma_struc.ma_sca_mod[i] = robust_sigma( useful_mod )
                ma_struc.ma_tmu_mod[i] = -2.50 * alog10( $
                    ma_struc.ma_mea_mod[i] / ( pix_area * exptime ) ) + magzpt  
                ma_struc.ma_amu_mod[i] = -2.50 * alog10( $
                    ma_struc.ma_med_mod[i] / ( pix_area * exptime ) ) + magzpt
                ma_struc.ma_sta_ori[i] = 1 
            endif
            ;; res 
            ;ma_struc.ma_res[i] = ( ma_struc.ma_tmu_ori[i] - $
            ;    ma_struc.ma_tmu_mod[i] )
            ma_struc.ma_res[i] = ( ma_struc.ma_amu_ori[i] - $
                ma_struc.ma_amu_mod[i] )

            ;; for the components 
            for k = 0, ( ncomp - 1 ), 1 do begin 

                ;; positive side 
                comps_box_p = $
                    comps_rot[k].img[ ma_struc.ma_x0p[i]:ma_struc.ma_x1p[i], $
                    ma_struc.ma_y0p[i]:ma_struc.ma_y2p[i] ] 
                ma_comp[k].ma_tot_p[i] = total( comps_box_p, /double ) 
                ma_comp[k].ma_mea_p[i] = ( ma_comp[k].ma_tot_p[i] / $
                    ma_struc.ma_pix_p[i] )
                ma_comp[k].ma_sca_p[i] = robust_sigma( comps_box_p ) 
                ma_comp[k].ma_med_p[i] = median( comps_box_p, /double ) 
                ma_comp[k].ma_tmu_p[i] = -2.50 * alog10( $
                    ma_comp[k].ma_mea_p[i] / ( exptime * pix_area ) ) + magzpt 
                ma_comp[k].ma_amu_p[i] = -2.50 * alog10( $
                    ma_comp[k].ma_med_p[i] / ( exptime * pix_area ) ) + magzpt 

                ;; negative side 
                comps_box_n = $
                    comps_rot[k].img[ ma_struc.ma_x1n[i]:ma_struc.ma_x0n[i], $
                    ma_struc.ma_y0n[i]:ma_struc.ma_y2n[i] ] 
                ma_comp[k].ma_tot_n[i] = total( comps_box_n, /double ) 
                ma_comp[k].ma_mea_n[i] = ( ma_comp[k].ma_tot_n[i] / $
                    ma_struc.ma_pix_n[i] )
                ma_comp[k].ma_sca_n[i] = robust_sigma( comps_box_n ) 
                ma_comp[k].ma_med_n[i] = median( comps_box_n, /double ) 
                ma_comp[k].ma_tmu_n[i] = -2.50 * alog10( $
                    ma_comp[k].ma_mea_n[i] / ( exptime * pix_area ) ) + magzpt 
                ma_comp[k].ma_amu_n[i] = -2.50 * alog10( $
                    ma_comp[k].ma_med_n[i] / ( exptime * pix_area ) ) + magzpt 

                index_p = where( comps_box_p gt -99999 ) 
                index_n = where( comps_box_n gt -99999 )
                comps_pix_p = comps_box_p[ index_p ]
                comps_pix_n = comps_box_n[ index_n ]
                ;;average 
                ma_comp[k].ma_tot[i] = ( ma_comp[k].ma_tot_p[i] + $
                    ma_comp[k].ma_tot_n[i] ) 
                ma_comp[k].ma_mea[i] = ( ma_comp[k].ma_tot[i] / ( $ 
                    ma_struc.ma_pix_p[i] + ma_struc.ma_pix_n[i] ) )
                ma_comp[k].ma_med[i] = median( [ comps_pix_p, comps_pix_n ], $
                    /double )
                ma_comp[k].ma_sca[i] = robust_sigma( [ comps_pix_p, $ 
                    comps_pix_n ] ) 
                ma_comp[k].ma_tmu[i] = -2.50 * alog10( ma_comp[k].ma_mea[i] / $ 
                    ( pix_area * exptime ) ) + magzpt
                ma_comp[k].ma_amu[i] = -2.50 * alog10( ma_comp[k].ma_med[i] / $ 
                    ( pix_area * exptime ) ) + magzpt
            endfor

            ;; decide whether or not go another round 
            nround = nround + 1

        endwhile 

    endelse

    ;; print 
    print, '##################################################################' 
    ;print, ma_struc.ma_pix_p[i], ma_struc.ma_use_p[i], ma_struc.ma_pix_n[i], $
    ;    ma_struc.ma_use_n[i]     
    ;print, ma_struc.ma_dis_p[i]
    ;print, ma_tot_ori_p[i], ma_tot_mod_p[i], ma_tot_ori_n[i], ma_tot_mod_n[i]
    ;print, ma_med_ori_p[i], ma_med_mod_p[i], ma_med_ori_n[i], ma_med_mod_n[i] 
    ;print, ma_sca_ori_p[i], ma_sca_mod_p[i], ma_sca_ori_n[i], ma_sca_mod_n[i] 
    ;print, ma_tmu_ori_p[i], ma_tmu_mod_p[i], ma_tmu_ori_n[i], ma_tmu_mod_n[i] 
    ;print, ma_struc.ma_amu_ori_p[i], ma_struc.ma_amu_mod_p[i], $
    ;    ma_struc.ma_amu_ori_n[i], ma_struc.ma_amu_ori_n[i] 
    ;print, '###############       Average         ################'
    ;print, ma_struc.ma_pix[i], ma_struc.ma_use[i] 
    ;print, ma_tot_ori[i], ma_tot_mod[i] 
    ;print, ma_med_ori[i], ma_med_mod[i] 
    ;print, ma_sca_ori[i], ma_sca_mod[i] 
    ;print, ma_tmu_ori[i], ma_tmu_mod[i] 
    ;print, ma_struc.ma_amu_ori[i], ma_struc.ma_amu_mod[i] 

endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define parameters for minor axis cut 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
mi_nbox = long( ( new_ycen - ( mi_height ) ) / mi_height ) 

mi_struc = { mi_x0p:lonarr( mi_nbox ), mi_y0p:lonarr( mi_nbox ), $ 
    mi_x0n:lonarr( mi_nbox ), mi_y0n:lonarr( mi_nbox ), $  
    mi_x1p:lonarr( mi_nbox ), mi_y1p:lonarr( mi_nbox ), $  
    mi_x1n:lonarr( mi_nbox ), mi_y1n:lonarr( mi_nbox ), $  
    mi_x2p:lonarr( mi_nbox ), mi_y2p:lonarr( mi_nbox ), $  
    mi_x2n:lonarr( mi_nbox ), mi_y2n:lonarr( mi_nbox ), $ 
    mi_x3p:lonarr( mi_nbox ), mi_y3p:lonarr( mi_nbox ), $  
    mi_x3n:lonarr( mi_nbox ), mi_y3n:lonarr( mi_nbox ), $
    mi_pix_p:lonarr( mi_nbox ), mi_use_p:lonarr( mi_nbox ), $ 
    mi_dis_p:fltarr( mi_nbox ), mi_fra_p:lonarr( mi_nbox ), $ 
    mi_pix_n:lonarr( mi_nbox ), mi_use_n:lonarr( mi_nbox ), $ 
    mi_dis_n:fltarr( mi_nbox ), mi_fra_n:lonarr( mi_nbox ), $
    mi_mea_ori_p:fltarr( mi_nbox ), mi_mea_ori_n:fltarr( mi_nbox ), $
    mi_tot_ori_p:fltarr( mi_nbox ), mi_tot_ori_n:fltarr( mi_nbox ), $
    mi_sca_ori_p:fltarr( mi_nbox ), mi_sca_ori_n:fltarr( mi_nbox ), $
    mi_med_ori_p:fltarr( mi_nbox ), mi_med_ori_n:fltarr( mi_nbox ), $
    mi_tmu_ori_p:fltarr( mi_nbox ), mi_tmu_ori_n:fltarr( mi_nbox ), $
    mi_amu_ori_p:fltarr( mi_nbox ), mi_amu_ori_n:fltarr( mi_nbox ), $
    mi_sta_ori_p:fltarr( mi_nbox ), mi_sta_ori_n:fltarr( mi_nbox ), $
    mi_mea_mod_p:fltarr( mi_nbox ), mi_mea_mod_n:fltarr( mi_nbox ), $
    mi_tot_mod_p:fltarr( mi_nbox ), mi_tot_mod_n:fltarr( mi_nbox ), $
    mi_sca_mod_p:fltarr( mi_nbox ), mi_sca_mod_n:fltarr( mi_nbox ), $
    mi_med_mod_p:fltarr( mi_nbox ), mi_med_mod_n:fltarr( mi_nbox ), $
    mi_tmu_mod_p:fltarr( mi_nbox ), mi_tmu_mod_n:fltarr( mi_nbox ), $
    mi_amu_mod_p:fltarr( mi_nbox ), mi_amu_mod_n:fltarr( mi_nbox ), $
    mi_res_p:fltarr( mi_nbox ), mi_res_n:fltarr( mi_nbox ), $ 
    mi_pix:lonarr( mi_nbox ), mi_use:lonarr( mi_nbox ), $
    mi_dis:lonarr( mi_nbox ), mi_fra:lonarr( mi_nbox ), $
    mi_mea_ori:fltarr( mi_nbox ), mi_mea_mod:fltarr( mi_nbox ), $ 
    mi_tot_ori:fltarr( mi_nbox ), mi_tot_mod:fltarr( mi_nbox ), $ 
    mi_med_ori:fltarr( mi_nbox ), mi_med_mod:fltarr( mi_nbox ), $ 
    mi_sca_ori:fltarr( mi_nbox ), mi_sca_mod:fltarr( mi_nbox ), $ 
    mi_tmu_ori:fltarr( mi_nbox ), mi_tmu_mod:fltarr( mi_nbox ), $ 
    mi_amu_ori:fltarr( mi_nbox ), mi_amu_mod:fltarr( mi_nbox ), $ 
    mi_sta_ori:fltarr( mi_nbox ), mi_sta_mod:fltarr( mi_nbox ), $ 
    mi_res: fltarr( mi_nbox ), mi_mod_log_flux:fltarr( mi_nbox ), $
    mi_reg_p:strarr( mi_nbox ), mi_reg_n:strarr( mi_nbox ) }
mi_comp = { mi_mea_p:fltarr( mi_nbox ), mi_tot_p:fltarr( mi_nbox ), $ 
    mi_med_p:fltarr( mi_nbox ), mi_sca_p:fltarr( mi_nbox ), $
    mi_tmu_p:fltarr( mi_nbox ), mi_amu_p:fltarr( mi_nbox ), $
    mi_mea_n:fltarr( mi_nbox ), mi_tot_n:fltarr( mi_nbox ), $ 
    mi_med_n:fltarr( mi_nbox ), mi_sca_n:fltarr( mi_nbox ), $
    mi_tmu_n:fltarr( mi_nbox ), mi_amu_n:fltarr( mi_nbox ), $
    mi_mea:fltarr( mi_nbox ), mi_tot:fltarr( mi_nbox ), $ 
    mi_med:fltarr( mi_nbox ), mi_sca:fltarr( mi_nbox ), $
    mi_tmu:fltarr( mi_nbox ), mi_amu:fltarr( mi_nbox )  }
mi_comp = replicate( mi_comp, ncomp )
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
for i = 0, ( mi_nbox - 1 ), 1 do begin  

    h = mi_height 
    c = mi_cenbox 
    w = mi_width 
    o = mi_offset
   s1 = mi_step1 
   s2 = mi_step2 

    print, '#########################################################'
    print, '#########' + string( i ) + '         ##########'
    print, '#########################################################'

    if ( i eq 0 ) then begin
        frac = 0.0 
        ;; central box 
        mi_struc.mi_x0p[i] = long( new_xcen - c / 2 )
        mi_struc.mi_y0p[i] = long( new_ycen - h / 2 )
        mi_struc.mi_x1p[i] = long( new_xcen + c / 2 ) 
        mi_struc.mi_y1p[i] = long( new_ycen - h / 2 )
        mi_struc.mi_x2p[i] = long( new_xcen + c / 2 )
        mi_struc.mi_y2p[i] = long( new_ycen + h / 2 )
        mi_struc.mi_x3p[i] = long( new_xcen - c / 2 )
        mi_struc.mi_y3p[i] = long( new_ycen + h / 2 )

        mi_struc.mi_x0n[i] = long( new_xcen - c / 2 )
        mi_struc.mi_y0n[i] = long( new_ycen - h / 2 )
        mi_struc.mi_x1n[i] = long( new_xcen + c / 2 )
        mi_struc.mi_y1n[i] = long( new_ycen - h / 2 )
        mi_struc.mi_x2n[i] = long( new_xcen + c / 2 )
        mi_struc.mi_y2n[i] = long( new_ycen + h / 2 )
        mi_struc.mi_x3n[i] = long( new_xcen - c / 2 )
        mi_struc.mi_y3n[i] = long( new_ycen + h / 2 )

        xx = string( ( mi_struc.mi_x0p[i] + mi_struc.mi_x1p[i] ) / 2.0 )
        yy = string( ( mi_struc.mi_y0p[i] + mi_struc.mi_y2p[i] ) / 2.0 )
        xw = string( abs( mi_struc.mi_x1p[i] - mi_struc.mi_x0p[i] ) )
        yw = string( abs( mi_struc.mi_y2p[i] - mi_struc.mi_y0p[i] ) )
        mi_struc.mi_reg_p[i] = 'box( ' + xx + ' , ' + yy + ' , ' + $
            xw + ' , ' + yw + ' )'

        xx = string( ( mi_struc.mi_x0n[i] + mi_struc.mi_x1n[i] ) / 2.0 )
        yy = string( ( mi_struc.mi_y0n[i] + mi_struc.mi_y2n[i] ) / 2.0 )
        xw = string( abs( mi_struc.mi_x1n[i] - mi_struc.mi_x0n[i] ) )
        yw = string( abs( mi_struc.mi_y2n[i] - mi_struc.mi_y0n[i] ) ) 
        mi_struc.mi_reg_n[i] = 'box( ' + xx + ' , ' + yy + ' , ' + $
            xw + ' , ' + yw + ' )'

        new_ycen_2 = ( float( mi_struc.mi_y0p[i] + mi_struc.mi_y2p[i] ) / 2.0 )

        mi_struc.mi_pix[i] = long( ( c + 1 ) * ( h + 1 ) ) 
        mi_struc.mi_pix_p[i] = mi_struc.mi_pix[i]
        mi_struc.mi_pix_n[i] = mi_struc.mi_pix[i] 

        mi_struc.mi_dis[i] = 0 
        mi_struc.mi_dis_p[i] = mi_struc.mi_dis[i] 
        mi_struc.mi_dis_n[i] = mi_struc.mi_dis[i] 

        if ( mi_struc.mi_x0p[i] ne mi_struc.mi_x1p[i] ) then begin 
            ori_box = ori_rot[ mi_struc.mi_x0p[i]:mi_struc.mi_x1p[i], $
                mi_struc.mi_y0p[i]:mi_struc.mi_y2p[i] ]  
            mod_box = mod_rot[ mi_struc.mi_x0p[i]:mi_struc.mi_x1p[i], $
                mi_struc.mi_y0p[i]:mi_struc.mi_y2p[i] ]  
            msk_box = msk_rot[ mi_struc.mi_x0p[i]:mi_struc.mi_x1p[i], $
                mi_struc.mi_y0p[i]:mi_struc.mi_y2p[i] ]  
        endif else begin 
            ori_box = ori_rot[ mi_struc.mi_x0p[i], $
                mi_struc.mi_y0p[i]:mi_struc.mi_y2p[i] ]  
            mod_box = mod_rot[ mi_struc.mi_x0p[i], $
                mi_struc.mi_y0p[i]:mi_struc.mi_y2p[i] ]  
            msk_box = msk_rot[ mi_struc.mi_x0p[i], $
                mi_struc.mi_y0p[i]:mi_struc.mi_y2p[i] ]  
        endelse

        index_use = where( msk_box lt 1 ) 

        if ( index_use[0] eq -1 ) then begin 
            frac = 0.0 
        endif else begin 
            num_use = n_elements( index_use ) 
            frac = float( num_use ) / float( mi_struc.mi_pix[i] ) 
            mi_struc.mi_use[i] = num_use 
            mi_struc.mi_use_p[i] = num_use
            mi_struc.mi_use_n[i] = num_use 

            ;; ori 
            sum_ori = total( ori_box[ index_use ], /double )
            mi_struc.mi_tot_ori[i] = ( sum_ori / num_use ) * $
                mi_struc.mi_pix[i]
            mi_struc.mi_mea_ori[i] = ( sum_ori / num_use ) 
            mi_struc.mi_med_ori[i] = median( ori_box[ index_use ], /double )
            useful_ori = ori_box[ index_use ]
            mi_struc.mi_sca_ori[i] = robust_sigma( useful_ori ) 
            mi_struc.mi_tmu_ori[i] = -2.50 * alog10( mi_struc.mi_mea_ori[i] / $
                ( pix_area * exptime ) ) + magzpt
            mi_struc.mi_amu_ori[i] = -2.50 * alog10( mi_struc.mi_med_ori[i] / $
                ( pix_area * exptime ) ) + magzpt

            ;; mod 
            sum_mod = total( mod_box[ index_use ], /double )
            mi_struc.mi_tot_mod[i] = ( sum_mod / num_use ) * $
                mi_struc.mi_pix[i]
            mi_struc.mi_mea_mod[i] = ( sum_mod / num_use ) 
            mi_struc.mi_med_mod[i] = median( mod_box[ index_use ], /double )
            useful = mod_box[ index_use ]
            mi_struc.mi_sca_mod[i] = robust_sigma( useful ) 
            mi_struc.mi_tmu_mod[i] = -2.50 * alog10( mi_struc.mi_mea_mod[i] / $
                ( pix_area * exptime ) ) + magzpt
            mi_struc.mi_amu_mod[i] = -2.50 * alog10( mi_struc.mi_med_mod[i] / $
                ( pix_area * exptime ) ) + magzpt

            ;; res
            ;mi_struc.mi_res[i] = ( mi_struc.mi_tmu_ori[i] - $
            ;    mi_struc.mi_tmu_mod[i] )
            mi_struc.mi_res[i] = ( mi_struc.mi_amu_ori[i] - $
                mi_struc.mi_amu_mod[i] )

            mi_struc.mi_mod_log_flux[i] = alog10( total( mod_box ) )

            mi_struc.mi_tot_ori_p[i] = mi_struc.mi_tot_ori[i]  
            mi_struc.mi_tot_ori_n[i] = mi_struc.mi_tot_ori[i] 
            mi_struc.mi_mea_ori_p[i] = mi_struc.mi_mea_ori[i]
            mi_struc.mi_mea_ori_n[i] = mi_struc.mi_mea_ori[i] 
            mi_struc.mi_med_ori_p[i] = mi_struc.mi_med_ori[i]  
            mi_struc.mi_med_ori_n[i] = mi_struc.mi_med_ori[i] 
            mi_struc.mi_sca_ori_p[i] = mi_struc.mi_sca_ori[i] 
            mi_struc.mi_sca_ori_n[i] = mi_struc.mi_sca_ori[i]
            mi_struc.mi_tmu_ori_p[i] = mi_struc.mi_tmu_ori[i] 
            mi_struc.mi_tmu_ori_n[i] = mi_struc.mi_tmu_ori[i] 
            mi_struc.mi_amu_ori_p[i] = mi_struc.mi_amu_ori[i]
            mi_struc.mi_amu_ori_n[i] = mi_struc.mi_amu_ori[i]

            mi_struc.mi_tot_mod_p[i] = mi_struc.mi_tot_mod[i]  
            mi_struc.mi_tot_mod_n[i] = mi_struc.mi_tot_mod[i] 
            mi_struc.mi_mea_mod_p[i] = mi_struc.mi_mea_mod[i]
            mi_struc.mi_mea_mod_n[i] = mi_struc.mi_mea_mod[i] 
            mi_struc.mi_med_mod_p[i] = mi_struc.mi_med_mod[i]  
            mi_struc.mi_med_mod_n[i] = mi_struc.mi_med_mod[i] 
            mi_struc.mi_sca_mod_p[i] = mi_struc.mi_sca_mod[i] 
            mi_struc.mi_sca_mod_n[i] = mi_struc.mi_sca_mod[i]
            mi_struc.mi_tmu_mod_p[i] = mi_struc.mi_tmu_mod[i] 
            mi_struc.mi_tmu_mod_n[i] = mi_struc.mi_tmu_mod[i] 
            mi_struc.mi_amu_mod_p[i] = mi_struc.mi_amu_mod[i]
            mi_struc.mi_amu_mod_n[i] = mi_struc.mi_amu_mod[i]

        endelse

        if ( frac gt bad_thre ) then begin 
            mi_struc.mi_sta_ori_p[i] = 1 
            mi_struc.mi_sta_ori_n[i] = 1 
            mi_struc.mi_sta_ori[i] = 1 
            mi_struc.mi_fra_p[i] = frac 
            mi_struc.mi_fra_n[i] = frac 
            mi_struc.mi_fra[i] = frac
        endif 

        ;; for the components 
        for k = 0, ( ncomp - 1 ), 1 do begin 

            ;; positive side 
            comps_box_p = $
                comps_rot[k].img[ mi_struc.mi_x0p[i]:mi_struc.mi_x1p[i], $
                mi_struc.mi_y0p[i]:mi_struc.mi_y2p[i] ] 
            mi_comp[k].mi_tot_p[i] = total( comps_box_p, /double ) 
            mi_comp[k].mi_mea_p[i] = ( mi_comp[k].mi_tot_p[i] / $
                mi_struc.mi_pix_p[i] )
            mi_comp[k].mi_sca_p[i] = robust_sigma( comps_box_p ) 
            mi_comp[k].mi_med_p[i] = median( comps_box_p, /double ) 
            mi_comp[k].mi_tmu_p[i] = -2.50 * alog10( mi_comp[k].mi_mea_p[i] / $ 
                ( exptime * pix_area ) ) + magzpt 
            mi_comp[k].mi_amu_p[i] = -2.50 * alog10( mi_comp[k].mi_med_p[i] / $ 
                ( exptime * pix_area ) ) + magzpt 

            ;; negative side 
            comps_box_n = $
                comps_rot[k].img[ mi_struc.mi_x0n[i]:mi_struc.mi_x1n[i], $
                mi_struc.mi_y0n[i]:mi_struc.mi_y2n[i] ] 
            mi_comp[k].mi_tot_n[i] = total( comps_box_n, /double ) 
            mi_comp[k].mi_mea_n[i] = ( mi_comp[k].mi_tot_n[i] / $
                mi_struc.mi_pix_n[i] )
            mi_comp[k].mi_sca_n[i] = robust_sigma( comps_box_n ) 
            mi_comp[k].mi_med_n[i] = median( comps_box_n, /double ) 
            mi_comp[k].mi_tmu_n[i] = -2.50 * alog10( mi_comp[k].mi_mea_n[i] / $ 
                ( exptime * pix_area ) ) + magzpt 
            mi_comp[k].mi_amu_n[i] = -2.50 * alog10( mi_comp[k].mi_med_n[i] / $ 
                ( exptime * pix_area ) ) + magzpt 

           index_p = where( comps_box_p gt -99999 ) 
           index_n = where( comps_box_n gt -99999 )
           comps_pix_p = comps_box_p[ index_p ]
           comps_pix_n = comps_box_n[ index_n ]
            ;;average 
            mi_comp[k].mi_tot[i] = ( mi_comp[k].mi_tot_p[i] + $
                mi_comp[k].mi_tot_n[i] ) 
            mi_comp[k].mi_mea[i] = ( mi_comp[k].mi_tot[i] / ( $ 
                mi_struc.mi_pix_p[i] + mi_struc.mi_pix_n[i] ) )
            mi_comp[k].mi_med[i] = median( [ comps_pix_p, comps_pix_n ], $
                /double )
            mi_comp[k].mi_sca[i] = robust_sigma( [ comps_pix_p, comps_pix_n ] ) 
            mi_comp[k].mi_tmu[i] = -2.50 * alog10( mi_comp[k].mi_mea[i] / $ 
                ( pix_area * exptime ) ) + magzpt
            mi_comp[k].mi_amu[i] = -2.50 * alog10( mi_comp[k].mi_med[i] / $ 
                ( pix_area * exptime ) ) + magzpt
        endfor

    ;; for non-central box 
    endif else begin 

        ;; the width of the box should increase stepwise 
        ;; now use R50 and R80 as standard 
        if ( ( mi_struc.mi_dis[ i-1 ] + mi_width ) gt ( r50 * 0.3 * ell ) ) $
            then begin 
            h = h + s1 
            w = w + s2 
        endif 
        if ( ( mi_struc.mi_dis[ i-1 ] + mi_width ) gt ( r50 * 0.6 * ell ) ) $
            then begin 
            h = h + s1 
            w = w + s2 
        endif 
        if ( ( mi_struc.mi_dis[ i-1 ] + mi_width ) gt ( r50 * ell ) ) then begin 
            h = h + s1 
            w = w + s2
        endif
        if ( ( mi_struc.mi_dis[ i-1 ] + mi_width ) gt $
            ( ( r80 + r50 ) / 2.0 * ell ) ) then begin 
            h = h + s1 
            w = w + s2
        endif
        if ( ( mi_struc.mi_dis[ i-1 ] + mi_width ) gt ( r80 * ell ) ) $
            then begin 
            h = h + s1 
            w = w + s2
        endif

        frac_p = 0.0
        frac_n = 0.0
        nround = 0
        flux_step = 0 

        while ( ( flux_step eq 0 ) and ( nround le 5 ) ) do begin 

            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ;; positive side  
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            mi_struc.mi_x0p[i] = ( long( new_xcen - w / 2 ) ) 
            mi_struc.mi_y0p[i] = ( long( mi_struc.mi_y2p[i-1] ) ) < $
                ( naxis2 - 1 ) 
            mi_struc.mi_x1p[i] = ( long( new_xcen + w / 2 ) )
            mi_struc.mi_y1p[i] = ( long( mi_struc.mi_y2p[i-1] ) ) < $
                ( naxis2 - 1 ) 
            mi_struc.mi_x2p[i] = ( long( new_xcen + w / 2 ) ) 
            mi_struc.mi_y2p[i] = ( long( mi_struc.mi_y2p[i-1] + h ) ) < $
                ( naxis2 - 1 ) 
            mi_struc.mi_x3p[i] = ( long( new_xcen - w / 2 ) ) 
            mi_struc.mi_y3p[i] = ( long( mi_struc.mi_y2p[i-1] + h ) ) < $
                ( naxis2 - 1 ) 
            mi_struc.mi_pix_p[i] = ( long( ( w + 1 ) * ( h + 1 ) ) ) 
            mi_struc.mi_dis_p[i] = abs( ( float( mi_struc.mi_y2p[i] + $
                mi_struc.mi_y1p[i] ) / 2.0 ) - new_ycen_2 - 1 ) 

            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ;; negative side 
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            mi_struc.mi_x0n[i] = ( long( new_xcen - w / 2 ) ) 
            mi_struc.mi_y0n[i] = ( long( mi_struc.mi_y2n[i-1] ) ) > 0.0
            mi_struc.mi_x1n[i] = ( long( new_xcen + w / 2 ) ) 
            mi_struc.mi_y1n[i] = ( long( mi_struc.mi_y2n[i-1] ) ) > 0.0
            mi_struc.mi_x2n[i] = ( long( new_xcen + w / 2 ) ) 
            mi_struc.mi_y2n[i] = ( long( mi_struc.mi_y2n[i-1] - h ) ) > 0.0
            mi_struc.mi_x3n[i] = ( long( new_xcen - w / 2 ) ) 
            mi_struc.mi_y3n[i] = ( long( mi_struc.mi_y2n[i-1] - h ) ) > 0.0
            mi_struc.mi_pix_n[i] = ( long( ( w + 1 ) * ( h + 1 ) ) ) 
            mi_struc.mi_dis_n[i] = abs( new_ycen_2 + 1 - ( $
                float( mi_struc.mi_y0n[i] + mi_struc.mi_y2n[i] ) / 2.0 ) )

            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            ;; save the region information 
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            xx = string( ( mi_struc.mi_x0p[i] + mi_struc.mi_x1p[i] ) / 2.0 )
            yy = string( ( mi_struc.mi_y0p[i] + mi_struc.mi_y2p[i] ) / 2.0 )
            xw = string( abs( mi_struc.mi_x1p[i] - mi_struc.mi_x0p[i] ) )
            yw = string( abs( mi_struc.mi_y2p[i] - mi_struc.mi_y0p[i] ) )
            mi_struc.mi_reg_p[i] = 'box( ' + xx + ' , ' + yy + ' , ' + $
                xw + ' , ' + yw + ' )'
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
            xx = string( ( mi_struc.mi_x0n[i] + mi_struc.mi_x1n[i] ) / 2.0 )
            yy = string( ( mi_struc.mi_y0n[i] + mi_struc.mi_y2n[i] ) / 2.0 )
            xw = string( abs( mi_struc.mi_x1n[i] - mi_struc.mi_x0n[i] ) )
            yw = string( abs( mi_struc.mi_y0n[i] - mi_struc.mi_y2n[i] ) ) 
            mi_struc.mi_reg_n[i] = 'box( ' + xx + ' , ' + yy + ' , ' + $
                xw + ' , ' + yw + ' )'

            ori_box_p = ori_rot[ mi_struc.mi_x0p[i]:mi_struc.mi_x1p[i], $
                mi_struc.mi_y0p[i]:mi_struc.mi_y2p[i] ]  
            mod_box_p = mod_rot[ mi_struc.mi_x0p[i]:mi_struc.mi_x1p[i], $
                mi_struc.mi_y0p[i]:mi_struc.mi_y2p[i] ]  
            msk_box_p = msk_rot[ mi_struc.mi_x0p[i]:mi_struc.mi_x1p[i], $
                mi_struc.mi_y0p[i]:mi_struc.mi_y2p[i] ]  
            index_use_p = where( msk_box_p lt 1 ) 
            num_use_p = n_elements( index_use_p )

            ori_box_n = ori_rot[ mi_struc.mi_x0n[i]:mi_struc.mi_x1n[i], $
                mi_struc.mi_y2n[i]:mi_struc.mi_y0n[i] ]  
            mod_box_n = mod_rot[ mi_struc.mi_x0n[i]:mi_struc.mi_x1n[i], $
                mi_struc.mi_y2n[i]:mi_struc.mi_y0n[i] ]  
            msk_box_n = msk_rot[ mi_struc.mi_x0n[i]:mi_struc.mi_x1n[i], $
                mi_struc.mi_y2n[i]:mi_struc.mi_y0n[i] ]  
            index_use_n = where( msk_box_n lt 1 ) 
            num_use_n = n_elements( index_use_n ) 

            mi_struc.mi_mod_log_flux[i] = alog10( ( $
                total( mod_box_n, /double ) + $
                total( mod_box_p, /double ) ) / 2.0 )
            if ( ( mi_struc.mi_mod_log_flux[i-1] - $
                mi_struc.mi_mod_log_flux[i] ) gt 0.50 ) then begin 
                flux_step = 0 
            endif else begin 
                flux_step = 1 
            endelse

            ;; %%%%%%%%% "and" v.s "or"
            if ( index_use_p[0] eq -1 ) then begin 
                frac_p = 0.0 
            endif else begin 
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ;; positive side 
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                frac_p = float( num_use_p ) / float( mi_struc.mi_pix_p[i] )
                mi_struc.mi_fra_p[i] = frac_p
                mi_struc.mi_use_p[i] = num_use_p 
                ;; ori 
                sum_ori = total( ori_box_p[ index_use_p ], /double ) 
                mi_struc.mi_tot_ori_p[i] = ( ( sum_ori / num_use_p ) * $
                    mi_struc.mi_pix_p[i] )
                mi_struc.mi_mea_ori_p[i] = ( sum_ori / num_use_p ) 
                mi_struc.mi_med_ori_p[i] = median( ori_box_p[ index_use_p ], $
                    /double )
                useful_ori_p = ori_box_p[ index_use_p ]
                mi_struc.mi_sca_ori_p[i] = robust_sigma( useful_ori_p ) 
                if ( mi_struc.mi_mea_ori_p[i] lt 0.0 ) then begin 
                    mi_struc.mi_tmu_ori_p[i] = 99.0
                endif else begin 
                    mi_struc.mi_tmu_ori_p[i] = -2.50 * alog10( $
                        mi_struc.mi_mea_ori_p[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse
                if ( mi_struc.mi_med_ori_p[i] lt 0.0 ) then begin 
                    mi_struc.mi_amu_ori_p[i] = 99.0 
                endif else begin 
                    mi_struc.mi_amu_ori_p[i] = -2.50 * alog10( $
                        mi_struc.mi_med_ori_p[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse

                ;; mod 
                sum_mod = total( mod_box_p[ index_use_p ] )
                mi_struc.mi_tot_mod_p[i] = ( ( sum_mod / num_use_p ) * $ 
                    mi_struc.mi_pix_p[i] )
                mi_struc.mi_mea_mod_p[i] = ( sum_mod / num_use_p ) 
                mi_struc.mi_med_mod_p[i] = median( mod_box_p[ index_use_p ], $
                    /double )
                useful_mod_p = mod_box_p[ index_use_p ]
                mi_struc.mi_sca_mod_p[i] = robust_sigma( useful_mod_p ) 
                if ( mi_struc.mi_mea_mod_p[i] lt 0.0 ) then begin 
                    mi_struc.mi_tmu_mod_p[i] = 99.0
                endif else begin 
                    mi_struc.mi_tmu_mod_p[i] = -2.50 * alog10( $
                        mi_struc.mi_mea_mod_p[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse
                if ( mi_struc.mi_med_mod_p[i] lt 0.0 ) then begin 
                    mi_struc.mi_amu_mod_p[i] = 99.0 
                endif else begin 
                    mi_struc.mi_amu_mod_p[i] = -2.50 * alog10( $
                        mi_struc.mi_med_mod_p[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse
                ;; res 
                ;mi_struc.mi_res_p[i] = ( mi_struc.mi_tmu_ori_p[i] - $
                ;    mi_struc.mi_tmu_mod_p[i] )
                mi_struc.mi_res_p[i] = ( mi_struc.mi_amu_ori_p[i] - $
                    mi_struc.mi_amu_mod_p[i] )
            endelse

            if ( index_use_n[0] eq -1 ) then begin 
                frac_n = 0.0 
            endif else begin 
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                ;; negative side 
                ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                frac_n = float( num_use_n ) / float( mi_struc.mi_pix_n[i] )
                mi_struc.mi_fra_n[i] = frac_n
                mi_struc.mi_use_n[i] = num_use_n 
                ;; ori 
                sum_ori = total( ori_box_n[ index_use_n ], /double ) 
                mi_struc.mi_tot_ori_n[i] = ( ( sum_ori / num_use_n ) * $
                    mi_struc.mi_pix_n[i] )
                mi_struc.mi_mea_ori_n[i] = ( sum_ori / num_use_n ) 
                mi_struc.mi_med_ori_n[i] = median( ori_box_n[ index_use_n ], $
                    /double )
                useful_ori_n = ori_box_n[ index_use_n ]
                mi_struc.mi_sca_ori_n[i] = robust_sigma( useful_ori_n ) 
                if ( mi_struc.mi_mea_ori_n[i] lt 0.0 ) then begin 
                    mi_struc.mi_tmu_ori_n[i] = 99.0
                endif else begin 
                    mi_struc.mi_tmu_ori_n[i] = -2.50 * alog10( $
                        mi_struc.mi_mea_ori_n[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse
                if ( mi_struc.mi_med_ori_n[i] lt 0.0 ) then begin 
                    mi_struc.mi_amu_ori_n[i] = 99.0 
                endif else begin 
                    mi_struc.mi_amu_ori_n[i] = -2.50 * alog10( $
                        mi_struc.mi_med_ori_n[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse

                ;; mod 
                sum_mod = total( mod_box_n[ index_use_n ] )
                mi_struc.mi_tot_mod_n[i] = ( ( sum_mod / num_use_n ) * $ 
                    mi_struc.mi_pix_n[i] )
                mi_struc.mi_mea_mod_n[i] = ( sum_mod / num_use_n ) 
                mi_struc.mi_med_mod_n[i] = median( mod_box_n[ index_use_n ], $
                    /double )
                useful_mod_n = mod_box_n[ index_use_n ]
                mi_struc.mi_sca_mod_n[i] = robust_sigma( useful_mod_n ) 
                if ( mi_struc.mi_mea_mod_n[i] lt 0.0 ) then begin 
                    mi_struc.mi_tmu_mod_n[i] = 99.0
                endif else begin 
                    mi_struc.mi_tmu_mod_n[i] = -2.50 * alog10( $
                        mi_struc.mi_mea_mod_n[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse
                if ( mi_struc.mi_med_mod_n[i] lt 0.0 ) then begin 
                    mi_struc.mi_amu_mod_n[i] = 99.0 
                endif else begin 
                    mi_struc.mi_amu_mod_n[i] = -2.50 * alog10( $
                        mi_struc.mi_med_mod_n[i] / ( pix_area * exptime ) ) $
                        + magzpt
                endelse

                ;; res 
                ;mi_struc.mi_res_n[i] = ( mi_struc.mi_tmu_ori_n[i] - $
                ;    mi_struc.mi_tmu_mod_n[i] )
                mi_struc.mi_res_n[i] = ( mi_struc.mi_amu_ori_n[i] - $
                    mi_struc.mi_amu_mod_n[i] )
            endelse
            
            ;; status of the box 
            if ( frac_p ge bad_thre ) then begin 
                mi_struc.mi_sta_ori_p[i] = 1 
            endif else begin 
                mi_struc.mi_sta_ori_p[i] = 0 
            endelse

            if ( frac_n ge bad_thre ) then begin 
                mi_struc.mi_sta_ori_n[i] = 1 
            endif else begin 
                mi_struc.mi_sta_ori_n[i] = 0 
            endelse

            ;; if the good pixel fractions are low at both side 
            ;; %%%%%% "and" or "or" !!!!!!!!!!
            if ( ( frac_p lt bad_thre ) and ( frac_n lt bad_thre ) ) then begin 
                flux_step = 0 
            endif
            ;; 
            if ( flux_step eq 0 ) then begin 
                h = h + s1
            endif 

            ;; average     
            ;; four situation 
            mi_struc.mi_dis[i] = mi_struc.mi_dis_p[i]
            ;; 1
            if ( ( frac_p ge bad_thre ) and ( frac_n lt bad_thre ) ) then begin 
                mi_struc.mi_pix[i] = mi_struc.mi_pix_p[i] 
                mi_struc.mi_use[i] = mi_struc.mi_use_p[i] 
                mi_struc.mi_fra[i] = mi_struc.mi_fra_p[i] 
                mi_struc.mi_tot_ori[i] = mi_struc.mi_tot_ori_p[i] 
                mi_struc.mi_mea_ori[i] = mi_struc.mi_mea_ori_p[i]
                mi_struc.mi_med_ori[i] = mi_struc.mi_med_ori_p[i] 
                mi_struc.mi_sca_ori[i] = mi_struc.mi_sca_ori_p[i] 
                mi_struc.mi_tmu_ori[i] = mi_struc.mi_tmu_ori_p[i] 
                mi_struc.mi_amu_ori[i] = mi_struc.mi_tmu_ori_p[i] 
                mi_struc.mi_tot_mod[i] = mi_struc.mi_tot_mod_p[i] 
                mi_struc.mi_mea_mod[i] = mi_struc.mi_mea_mod_p[i]
                mi_struc.mi_med_mod[i] = mi_struc.mi_med_mod_p[i] 
                mi_struc.mi_sca_mod[i] = mi_struc.mi_sca_mod_p[i] 
                mi_struc.mi_tmu_mod[i] = mi_struc.mi_tmu_mod_p[i] 
                mi_struc.mi_amu_mod[i] = mi_struc.mi_tmu_mod_p[i] 
                mi_struc.mi_sta_ori[i] = 1 
            endif
            ;; 2
            if ( ( frac_p lt bad_thre ) and ( frac_n ge bad_thre ) ) then begin 
                mi_struc.mi_pix[i] = mi_struc.mi_pix_n[i] 
                mi_struc.mi_use[i] = mi_struc.mi_use_n[i] 
                mi_struc.mi_fra[i] = mi_struc.mi_fra_n[i] 
                mi_struc.mi_tot_ori[i] = mi_struc.mi_tot_ori_n[i] 
                mi_struc.mi_mea_ori[i] = mi_struc.mi_mea_ori_n[i]
                mi_struc.mi_med_ori[i] = mi_struc.mi_med_ori_n[i] 
                mi_struc.mi_sca_ori[i] = mi_struc.mi_sca_ori_n[i] 
                mi_struc.mi_tmu_ori[i] = mi_struc.mi_tmu_ori_n[i] 
                mi_struc.mi_amu_ori[i] = mi_struc.mi_tmu_ori_n[i] 
                mi_struc.mi_tot_mod[i] = mi_struc.mi_tot_mod_n[i] 
                mi_struc.mi_med_mod[i] = mi_struc.mi_med_mod_n[i] 
                mi_struc.mi_med_mod[i] = mi_struc.mi_med_mod_n[i] 
                mi_struc.mi_sca_mod[i] = mi_struc.mi_sca_mod_n[i] 
                mi_struc.mi_tmu_mod[i] = mi_struc.mi_tmu_mod_n[i] 
                mi_struc.mi_amu_mod[i] = mi_struc.mi_tmu_mod_n[i] 
                mi_struc.mi_sta_ori[i] = 1 
            endif
            ;; 3 
            if ( ( frac_p lt bad_thre ) and ( frac_n lt bad_thre ) ) then begin 
                mi_struc.mi_pix[i] = mi_struc.mi_pix_n[i] 
                mi_struc.mi_use[i] = 0 
                mi_struc.mi_fra[i] = 0.0
                mi_struc.mi_tot_ori[i] = 0.0
                mi_struc.mi_mea_ori[i] = 0.0 
                mi_struc.mi_med_ori[i] = 0.0
                mi_struc.mi_sca_ori[i] = 0.0
                mi_struc.mi_tmu_ori[i] = 99.0
                mi_struc.mi_amu_ori[i] = 99.0
                mi_struc.mi_tot_mod[i] = 0.0
                mi_struc.mi_med_mod[i] = 0.0
                mi_struc.mi_med_mod[i] = 0.0
                mi_struc.mi_sca_mod[i] = 0.0
                mi_struc.mi_tmu_mod[i] = 99.0
                mi_struc.mi_amu_mod[i] = 99.0
                mi_struc.mi_sta_ori[i] = 0 
            endif
            ;; 4 
            if ( ( frac_p ge bad_thre ) and ( frac_n ge bad_thre ) ) then begin 
                mi_struc.mi_pix[i] = ( mi_struc.mi_pix_p[i] + $
                    mi_struc.mi_pix_n[i] )
                mi_struc.mi_use[i] = ( mi_struc.mi_use_p[i] + $
                    mi_struc.mi_use_n[i] )
                mi_struc.mi_fra[i] = ( float( mi_struc.mi_use[i] ) / $ 
                    float( mi_struc.mi_pix[i] ) )

                ;; ori 
                mi_struc.mi_tot_ori[i] = ( mi_struc.mi_tot_ori_p[i] + $
                    mi_struc.mi_tot_ori_n[i] ) 
                useful_ori = [ useful_ori_p, useful_ori_n ] 
                sum_ori = total( useful_ori, /double ) 
                mi_struc.mi_mea_ori[i] = ( sum_ori / mi_struc.mi_use[i] ) 
                mi_struc.mi_med_ori[i] = median( useful_ori, /double ) 
                mi_struc.mi_sca_ori[i] = robust_sigma( useful_ori )
                mi_struc.mi_tmu_ori[i] = -2.50 * alog10( $
                    mi_struc.mi_mea_ori[i] / ( pix_area * exptime ) ) + magzpt  
                mi_struc.mi_amu_ori[i] = -2.50 * alog10( $
                    mi_struc.mi_med_ori[i] / ( pix_area * exptime ) ) + magzpt

                ;; mod 
                mi_struc.mi_tot_mod[i] = ( mi_struc.mi_tot_mod_p[i] + $
                    mi_struc.mi_tot_mod_n[i] ) 
                useful_mod = [ useful_mod_p, useful_mod_n ] 
                sum_mod = total( useful_mod, /double ) 
                mi_struc.mi_mea_mod[i] = ( sum_mod / mi_struc.mi_use[i] ) 
                mi_struc.mi_med_mod[i] = median( useful_mod, /double ) 
                mi_struc.mi_sca_mod[i] = robust_sigma( useful_mod )
                mi_struc.mi_tmu_mod[i] = -2.50 * alog10( $
                    mi_struc.mi_mea_mod[i] / ( pix_area * exptime ) ) + magzpt  
                mi_struc.mi_amu_mod[i] = -2.50 * alog10( $
                    mi_struc.mi_med_mod[i] / ( pix_area * exptime ) ) + magzpt
                mi_struc.mi_sta_ori[i] = 1 
            endif
            ;; res 
            ;mi_struc.mi_res[i] = ( mi_struc.mi_tmu_ori[i] - $
            ;    mi_struc.mi_tmu_mod[i] )
            mi_struc.mi_res[i] = ( mi_struc.mi_amu_ori[i] - $
                mi_struc.mi_amu_mod[i] )

            ;; for the components 
            for k = 0, ( ncomp - 1 ), 1 do begin 

                ;; positive side 
                comps_box_p = $
                    comps_rot[k].img[ mi_struc.mi_x0p[i]:mi_struc.mi_x1p[i], $
                    mi_struc.mi_y0p[i]:mi_struc.mi_y2p[i] ] 
                mi_comp[k].mi_tot_p[i] = total( comps_box_p, /double ) 
                mi_comp[k].mi_mea_p[i] = ( mi_comp[k].mi_tot_p[i] / $
                    mi_struc.mi_pix_p[i] )
                mi_comp[k].mi_sca_p[i] = robust_sigma( comps_box_p ) 
                mi_comp[k].mi_med_p[i] = median( comps_box_p, /double ) 
                mi_comp[k].mi_tmu_p[i] = -2.50 * alog10( $
                    mi_comp[k].mi_mea_p[i] / ( exptime * pix_area ) ) + magzpt 
                mi_comp[k].mi_amu_p[i] = -2.50 * alog10( $
                    mi_comp[k].mi_med_p[i] / ( exptime * pix_area ) ) + magzpt 

                ;; negative side 
                comps_box_n = $
                    comps_rot[k].img[ mi_struc.mi_x0n[i]:mi_struc.mi_x1n[i], $
                    mi_struc.mi_y2n[i]:mi_struc.mi_y0n[i] ] 
                mi_comp[k].mi_tot_n[i] = total( comps_box_n, /double ) 
                mi_comp[k].mi_mea_n[i] = ( mi_comp[k].mi_tot_n[i] / $
                    mi_struc.mi_pix_n[i] )
                mi_comp[k].mi_sca_n[i] = robust_sigma( comps_box_n ) 
                mi_comp[k].mi_med_n[i] = median( comps_box_n, /double ) 
                mi_comp[k].mi_tmu_n[i] = -2.50 * alog10( $
                    mi_comp[k].mi_mea_n[i] / ( exptime * pix_area ) ) + magzpt 
                mi_comp[k].mi_amu_n[i] = -2.50 * alog10( $
                    mi_comp[k].mi_med_n[i] / ( exptime * pix_area ) ) + magzpt 

                index_p = where( comps_box_p gt -99999 )
                index_n = where( comps_box_n gt -99999 )
                comps_pix_p = comps_box_p[ index_p ] 
                comps_pix_n = comps_box_n[ index_n ]
                ;;average 
                mi_comp[k].mi_tot[i] = ( mi_comp[k].mi_tot_p[i] + $
                    mi_comp[k].mi_tot_n[i] ) 
                mi_comp[k].mi_mea[i] = ( mi_comp[k].mi_tot[i] / ( $ 
                    mi_struc.mi_pix_p[i] + mi_struc.mi_pix_n[i] ) )
                mi_comp[k].mi_med[i] = median( [ comps_pix_p, comps_pix_n ], $
                    /double )
                mi_comp[k].mi_sca[i] = robust_sigma( [ comps_pix_p, $ 
                    comps_pix_n ] ) 
                mi_comp[k].mi_tmu[i] = -2.50 * alog10( mi_comp[k].mi_mea[i] / $ 
                    ( pix_area * exptime ) ) + magzpt
                mi_comp[k].mi_amu[i] = -2.50 * alog10( mi_comp[k].mi_med[i] / $ 
                    ( pix_area * exptime ) ) + magzpt
            endfor

            ;; decide whether or not go another round 
            nround = nround + 1

        endwhile 

    endelse

    ;; print 
    print, '##################################################################' 
    ;print, mi_struc.mi_pix_p[i], mi_struc.mi_use_p[i], mi_struc.mi_pix_n[i], $
    ;    mi_struc.mi_use_n[i]     
    ;print, mi_struc.mi_dis_p[i]
    ;print, mi_tot_ori_p[i], mi_tot_mod_p[i], mi_tot_ori_n[i], mi_tot_mod_n[i]
    ;print, mi_med_ori_p[i], mi_med_mod_p[i], mi_med_ori_n[i], mi_med_mod_n[i] 
    ;print, mi_sca_ori_p[i], mi_sca_mod_p[i], mi_sca_ori_n[i], mi_sca_mod_n[i] 
    ;print, mi_tmu_ori_p[i], mi_tmu_mod_p[i], mi_tmu_ori_n[i], mi_tmu_mod_n[i] 
    ;print, mi_struc.mi_amu_ori_p[i], mi_struc.mi_amu_mod_p[i], $
    ;    mi_struc.mi_amu_ori_n[i], mi_struc.mi_amu_ori_n[i] 
    ;print, '###############       Average         ################'
    ;print, mi_struc.mi_pix[i], mi_struc.mi_use[i] 
    ;print, mi_tot_ori[i], mi_tot_mod[i] 
    ;print, mi_med_ori[i], mi_med_mod[i] 
    ;print, mi_sca_ori[i], mi_sca_mod[i] 
    ;print, mi_tmu_ori[i], mi_tmu_mod[i] 
    ;print, mi_struc.mi_amu_ori[i], mi_struc.mi_amu_mod[i] 

endfor


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;make plots 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

position1 = [ 0.10, 0.35, 0.54, 0.99 ] 
;position1 = [ 0.10, 0.35, 0.94, 0.99 ] 
position2 = [ 0.55, 0.35, 0.99, 0.99 ] 
position3 = [ 0.10, 0.15, 0.54, 0.345 ]
;position3 = [ 0.10, 0.15, 0.94, 0.345 ]
position4 = [ 0.55, 0.15, 0.99, 0.345 ]

;xtitle = 'RSMA (arcsec^0.25)'
xtitle = 'Radius (arcsec)'
ytitle1 = 'Surface Brightness (mag/arcsec^2)' 
ytitle2 = 'Residual' 

;; radius range 
rad_range1 = ( r80 * 3.5 * 0.259 ) < ( ( new_xcen - 2 * ma_width ) * 0.259 )  
xrange1 = [ 0.01, rad_range1 ]
rad_range2 = ( r80 * ell * 1.3 * 0.259 ) < $
    ( ( new_ycen - 2 * ma_height ) * 0.259 )
xrange2 = [ 0.01, rad_range2 ] 

width1 = ( xrange1[1] - xrange1[0] )
width2 = ( xrange2[1] - xrange2[0] )

;; tick size
;; major axis
if ( width1 le 50 ) then begin 
    xtick1 = 10 
endif 
if ( width1 gt 50 ) then begin 
    xtick1 = 20 
endif 
if ( width1 gt 120 ) then begin 
    xtick1 = 30 
endif 
if ( width1 gt 160 ) then begin 
    xtick1 = 40 
endif 
if ( width1 gt 220 ) then begin 
    xtick1 = 50 
endif 
if ( width1 gt 280 ) then begin 
    xtick1 = 60 
endif 
;; minor axis
if ( width2 le 50 ) then begin 
    xtick2 = 10 
endif 
if ( width2 gt 50 ) then begin 
    xtick2 = 20 
endif 
if ( width2 gt 120 ) then begin 
    xtick2 = 30 
endif 
if ( width2 gt 160 ) then begin 
    xtick2 = 40 
endif 
if ( width2 gt 220 ) then begin 
    xtick2 = 50 
endif 
if ( width2 gt 280 ) then begin 
    xtick2 = 60 
endif 

;; surface brightness range
min_sbp = ma_struc.ma_tmu_mod[0] < mi_struc.mi_tmu_mod[0]
sbp_range1 = ( min_sbp - 0.3 )  
sbp_range2 = 27.40  
yrange1 = [ sbp_range2, sbp_range1 ] 
print, 'yrange : ', sbp_range1, sbp_range2
yrange2 = [ -0.399, 0.399 ] 

psxsize = 26 
psysize = 18

case ncomp of 
    2: colors = [ 'RED', 'BLUE' ]
    3: colors = [ 'RED', 'GREEN', 'BLUE' ] 
    4: colors = [ 'RED', 'GREEN', 'ORANGE', 'BLUE' ] 
    5: colors = [ 'RED', 'GREEN', 'ORANGE', 'BROWN', 'BLUE' ] 
    6: colors = [ 'RED', 'GREEN', 'ORANGE', 'BROWN', 'SEA GREEN', 'BLUE' ]
endcase

;;
plot_1 = model_string + '_ma_mi.eps'

!P.Font=0
set_plot, 'PS'
device, filename=plot_1, font_size=7.5, /encapsulated, $
    /color, set_font='HELVETICA BOLD', /tt_font, /bold, xsize=psxsize, $
    ysize=psysize

ma_struc.ma_dis_p = ( ma_struc.ma_dis_p + 0.6 )
ma_struc.ma_dis_n = ( ma_struc.ma_dis_n + 0.6 )
ma_struc.ma_dis = ( ma_struc.ma_dis + 0.6 )

mi_struc.mi_dis_p = ( mi_struc.mi_dis_p + 0.6 )
mi_struc.mi_dis_n = ( mi_struc.mi_dis_n + 0.6 )
mi_struc.mi_dis = ( mi_struc.mi_dis + 0.6 )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; For major axis 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
sta_p = where( ma_struc.ma_sta_ori_p eq 1 )
sta_n = where( ma_struc.ma_sta_ori_n eq 1 ) 
; sta = where( ma_struc.ma_sta_ori eq 1 ) 
  sta = where( ( ma_struc.ma_sta_ori_p eq 1 ) or ( ma_struc.ma_sta_ori_n eq 1 ) )

;; surface brightness  
cgPlot, ( ma_struc.ma_dis_p[sta_p] * 0.259 ), ma_struc.ma_tmu_ori_p[sta_p], $
    xstyle=1, ystyle=1, xrange=xrange1, yrange=yrange1, /nodata, $
    ytitle=ytitle1, position=position1, /noerase, xtickinterval=xtick1, $
    charsize=2.0, charthick=3.5, ytickinterval=2.0, xthick=4.0, ythick=4.0, $ 
    xtickformat="(A1)" 

cgplot, ( ma_struc.ma_dis_p[sta_p] * 0.259 ), ma_struc.ma_tmu_ori_p[sta_p], $
    psym=9, symsize=1.2, thick=3.5, /overplot, $
    color=cgColor( 'MAGENTA', !D.Table_Size )  
cgplot, ( ma_struc.ma_dis_p[sta_p] * 0.259 ), ma_struc.ma_tmu_mod_p[sta_p], $
    psym=0, linestyle=1, thick=3.5, /overplot, $
    color=cgColor( 'MAGENTA', !D.Table_Size ) 

cgplot, ( ma_struc.ma_dis_n[sta_n] * 0.259 ), ma_struc.ma_tmu_ori_n[sta_n], $
    psym=4, symsize=1.2, thick=3.5, /overplot, $
    color=cgColor( 'GOLD', !D.Table_Size )  
cgplot, ( ma_struc.ma_dis_n[sta_n] * 0.259 ), ma_struc.ma_tmu_mod_n[sta_n], $
    psym=0, linestyle=1, thick=3.5, /overplot, $
    color=cgColor( 'GOLD', !D.Table_Size ) 

cgplot, ( ma_struc.ma_dis[sta] * 0.259 ), ma_struc.ma_tmu_ori[sta], psym=16, $
    symsize=1.1, thick=2.5, /overplot, $
    color=cgColor( 'DARK GRAY', !D.Table_Size )  
cgplot, ( ma_struc.ma_dis[sta] * 0.259 ), ma_struc.ma_tmu_mod[sta], psym=0, $
    linestyle=0, thick=2.4, /overplot, $
    color=cgColor( 'DARK GRAY', !D.Table_Size ) 

;; plot for components 
if ( ncomp ge 2 ) then begin 
    for k = 0, ( ncomp - 1 ), 1 do begin 
        cgplot, ( ma_struc.ma_dis * 0.259 ), ma_comp[k].ma_tmu, psym=0, $ 
            linestyle = ( k + 2 ), thick=3.5, /overplot, $
            color=cgColor( colors[k], !D.Table_Size ) 
    endfor
endif

aa = position1[2] - 0.20 
bb = position1[3] - 0.10 
cgText, aa, bb, 'Major Axis', /normal, charsize=2.5, charthick=3.5, $
    color=cgColor( 'BLACK', !D.Table_Size )

;; residual 
cgPlot, ( ma_struc.ma_dis_p[sta_p] * 0.259 ), ma_struc.ma_res_p[sta_p], $
    xstyle=1, ystyle=1, xrange=xrange1, yrange=yrange2, /nodata, $
    xtitle=xtitle, position=position3, /noerase, xtickinterval=xtick1, $
    charsize=2.0, charthick=3.5, ytickinterval=0.2, xthick=4.0, ythick=4.0, $ 
    ytitle=ytitle2

cgplot, ( ma_struc.ma_dis_p[sta_p] * 0.259 ), ma_struc.ma_res_p[sta_p], $
    psym=9, symsize=1.3, $
    thick=3.5, /overplot, color=cgColor( 'MAGENTA', !D.Table_Size )  
cgplot, ( ma_struc.ma_dis_n[sta_n] * 0.259 ), ma_struc.ma_res_n[sta_n], $
    psym=4, symsize=1.1, $
    thick=3.5, /overplot, color=cgColor( 'GOLD', !D.Table_Size ) 
cgplot, ( ma_struc.ma_dis[sta] * 0.259 ), ma_struc.ma_res[sta], $
    psym=16, symsize=1.3, $
    thick=3.5, /overplot, color=cgColor( 'DARK GRAY', !D.Table_Size ) 
cgPlot, !X.Crange, [ 0.0, 0.0], psym=0, linestyle=0, thick=3.0, $ 
    color=cgColor( 'GRAY', !D.Table_Size ), /overplot

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; For minor axis 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
sta_p = where( mi_struc.mi_sta_ori_p eq 1 )
sta_n = where( mi_struc.mi_sta_ori_n eq 1 ) 
; sta = where( mi_struc.mi_sta_ori eq 1 ) 
  sta = where( ( mi_struc.mi_sta_ori_p eq 1 ) or ( mi_struc.mi_sta_ori_n eq 1 ) )

;; surface brightness  
cgPlot, ( mi_struc.mi_dis_p[sta_p] * 0.259 ), mi_struc.mi_tmu_ori_p[sta_p], $
    xstyle=1, ystyle=1, xrange=xrange2, yrange=yrange1, /nodata, $
    ytickformat="(A1)", position=position2, /noerase, xtickinterval=xtick2, $
    charsize=2.0, charthick=3.5, ytickinterval=2.0, xthick=4.0, ythick=4.0, $ 
    xtickformat="(A1)" 

cgplot, ( mi_struc.mi_dis_p[sta_p] * 0.259 ), mi_struc.mi_tmu_ori_p[sta_p], $
    psym=9, symsize=1.2, thick=3.5, /overplot, $
    color=cgColor( 'MAGENTA', !D.Table_Size )  
cgplot, ( mi_struc.mi_dis_p[sta_p] * 0.259 ), mi_struc.mi_tmu_mod_p[sta_p], $
    psym=0, linestyle=1, thick=3.5, /overplot, $
    color=cgColor( 'MAGENTA', !D.Table_Size ) 

cgplot, ( mi_struc.mi_dis_n[sta_n] * 0.259 ), mi_struc.mi_tmu_ori_n[sta_n], $
    psym=4, symsize=1.2, thick=3.5, /overplot, $
    color=cgColor( 'GOLD', !D.Table_Size )  
cgplot, ( mi_struc.mi_dis_n[sta_n] * 0.259 ), mi_struc.mi_tmu_mod_n[sta_n], $
    psym=0, linestyle=1, thick=3.5, /overplot, $
    color=cgColor( 'GOLD', !D.Table_Size ) 

cgplot, ( mi_struc.mi_dis[sta] * 0.259 ), mi_struc.mi_tmu_ori[sta], psym=16, $
    symsize=1.1, thick=2.5, /overplot, $
    color=cgColor( 'DARK GRAY', !D.Table_Size )  
cgplot, ( mi_struc.mi_dis[sta] * 0.259 ), mi_struc.mi_tmu_mod[sta], psym=0, $
    linestyle=0, thick=2.4, /overplot, $
    color=cgColor( 'DARK GRAY', !D.Table_Size ) 

;; plot for components 
if ( ncomp ge 2 ) then begin 
    for k = 0, ( ncomp - 1 ), 1 do begin 
        cgplot, ( mi_struc.mi_dis * 0.259 ), mi_comp[k].mi_tmu, psym=0, $ 
            linestyle = ( k + 2 ), thick=3.5, /overplot, $
            color=cgColor( colors[k], !D.Table_Size ) 
    endfor
endif

aa = position2[2] - 0.20 
bb = position2[3] - 0.10 
cgText, aa, bb, 'Minor Axis', /normal, charsize=2.5, charthick=3.5, $
    color=cgColor( 'BLACK', !D.Table_Size )

;; residual 
cgPlot, ( mi_struc.mi_dis_p[sta_p] * 0.259 ), mi_struc.mi_res_p[sta_p], $
    xstyle=1, ystyle=1, xrange=xrange2, yrange=yrange2, /nodata, $
    xtitle=xtitle, position=position4, /noerase, xtickinterval=xtick2, $
    charsize=2.0, charthick=3.5, ytickinterval=0.2, xthick=4.0, ythick=4.0, $ 
    ytickformat="(A1)"

cgplot, ( mi_struc.mi_dis_p[sta_p] * 0.259 ), mi_struc.mi_res_p[sta_p], $
    psym=9, symsize=1.3, $
    thick=3.5, /overplot, color=cgColor( 'MAGENTA', !D.Table_Size )  
cgplot, ( mi_struc.mi_dis_n[sta_n] * 0.259 ), mi_struc.mi_res_n[sta_n], $
    psym=4, symsize=1.3, $
    thick=3.5, /overplot, color=cgColor( 'GOLD', !D.Table_Size ) 
cgplot, ( mi_struc.mi_dis[sta] * 0.259 ), mi_struc.mi_res[sta], $
    psym=16, symsize=1.1, $
    thick=3.5, /overplot, color=cgColor( 'DARK GRAY', !D.Table_Size ) 
cgPlot, !X.Crange, [ 0.0, 0.0], psym=0, linestyle=0, thick=3.0, $ 
    color=cgColor( 'GRAY', !D.Table_Size ), /overplot

device, /close

;; save the region file 
;;;;;;;;;;;;;;;;;;;;;;;;
;; for major cut 
;;;;;;;;;;;;;;;;;;;;;;;;
ma_reg_file = model_string + '_ma.reg' 
openw, 10, ma_reg_file, width=400 
printf, 10, '# Region file format: DS9 version 4.1'
printf, 10, '# Filename: NGC4486_g_over.fits'
printf, 10, 'global color=blue dashlist=8 3 width=3 font="helvetica 10 ' + $ 
    ' normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 ' + $ 
    ' delete=1 include=1 source=1'
printf, 10, 'image'
for j = 0, ( ma_nbox - 1 ), 1 do begin 
    if ( ma_struc.ma_sta_ori_p[j] eq 1 ) then begin 
        printf, 10, ma_struc.ma_reg_p[j] 
    endif 
    if ( ma_struc.ma_sta_ori_n[j] eq 1 ) then begin 
        printf, 10, ma_struc.ma_reg_n[j] 
    endif 
endfor
close, 10

;;;;;;;;;;;;;;;;;;;;;;;;
;; for minor cut 
;;;;;;;;;;;;;;;;;;;;;;;;
mi_reg_file = model_string + '_mi.reg' 
openw, 10, mi_reg_file, width=400 
printf, 10, '# Region file formit: DS9 version 4.1'
printf, 10, '# Filename: NGC4486_g_over.fits'
printf, 10, 'global color=yellow dashlist=8 3 width=3 font="helvetica 10 ' + $ 
    ' normil romin" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 ' + $ 
    ' delete=1 include=1 source=1'
printf, 10, 'image'
for j = 0, ( mi_nbox - 1 ), 1 do begin 
    if ( mi_struc.mi_sta_ori_p[j] eq 1 ) then begin 
        printf, 10, mi_struc.mi_reg_p[j] 
    endif 
    if ( mi_struc.mi_sta_ori_n[j] eq 1 ) then begin 
        printf, 10, mi_struc.mi_reg_n[j] 
    endif 
endfor
close, 10

print, mi_struc.mi_tmu_ori_p 
print, '#############'
print, mi_struc.mi_tmu_ori_n 
print, '#############' 
print, mi_struc.mi_tmu_ori 
print, '#############' 
print, mi_struc.mi_dis_p 
print, '#############' 
print, mi_struc.mi_dis_n
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
