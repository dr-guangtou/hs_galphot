pro object_psf_building, image, snr, min_fwhm = min_fwhm, max_fwhm = max_fwhm, $
    vgt_size = vgt_size, psf_size = psf_size, magzpt = magzpt, pixel = pixel, $
    satu_level = satu_level, bg_size = bg_size, seeing_fwhm = seeing_fwhm, $
    ellipse = ellipse, psf_sum=psf_sum 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; using SExtractor and PSFExtractor to build the empirical PSF image 
;;; written by Song Huang 24.05.2012 
;;; modified by Song Huang 30.05.2012
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
on_error, 2
compile_opt idl2

if N_params() lt 2 then begin 
    print,  'Syntax - Object_PSF_Building, Image, SNR, Min_FWHM, Max_FWHM, '
    print,  '         VGT_SIZE, PSF_SIZE, MagZpt, Pixel, Satu_Level, BG_Size, '
    print,  '         Seeing_FWHM, /Ellipse  '
    return
endif
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
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
snr = float( SNR )
print, 'SNR Threshold  : ' + string( snr )
temp = size( img )
naxis1 = temp[1]
naxis2 = temp[2]
print, 'IMAGE SIZE : ' + string( naxis1 ) +' x ' + string( naxis2 )
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Parameter List
if NOT keyword_set( bg_size ) then begin 
    bg_size = 512
endif else begin 
    bg_size = float( bg_size )
endelse
if NOT keyword_set( magzpt ) then begin 
    magzpt = 23.5
endif else begin 
    magzpt = float( magzpt ) 
endelse
if NOT keyword_set( pixel ) then begin 
    pixel = 0.259
endif else begin 
    pixel = float( pixel ) 
endelse
if NOT keyword_set( satu_level ) then begin 
    satu_level = 50000
endif else begin 
    satu_level = long( satu_level ) 
endelse
if NOT keyword_set( seeing_fwhm ) then begin 
    seeing_fwhm = 1.00
endif else begin 
    seeing_fwhm = float( seeing_fwhm ) 
endelse
if NOT keyword_set( psf_size ) then begin 
    psf_size = 60 
endif else begin 
    psf_size = float( psf_size ) 
endelse 
if NOT keyword_set( vgt_size ) then begin 
    vgt_size = 90 
endif else begin 
    vgt_size = long( vgt_size ) 
endelse
if NOT keyword_set( min_fwhm ) then begin 
    min_fwhm = 2.0
endif else begin 
    min_fwhm = float( min_fwhm ) 
endelse
if NOT keyword_set( max_fwhm ) then begin 
    max_fwhm = 26.0
endif else begin 
    max_fwhm = float( max_fwhm ) 
endelse
;; constant
thresh      = 1.5
bg_thick    = 24
maxellip    = 0.30
poly_degree = 0
nsnap       = 1
filter_name = 'default.conv'
;filter_name = 'tophat_1.5_3x3.conv'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 0. Modify the parameter file for SExtractor 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
template = 'shuang_psf.param'
if NOT file_test( template ) then begin 
    message, 'Can not find shuang_psf.param file ! '
endif else begin 
    num_lines = file_lines( template )
    para = strarr( num_lines ) 
    openr, 10, template
    readf, 10, para 
    close, 10 
endelse
;; Change the vignet size 
index = where( strpos( para, 'VIGNET' ) ne -1 ) 
vgt_string = strcompress( string( vgt_size ), /remove_all )
para[ index ] = 'VIGNET(' + vgt_string + ',' + vgt_string + ')' 
;; Save the new parameter file 
openw, 20, template, width = 600 
for i = 0, ( num_lines - 1 ), 1 do begin 
    printf, 20, para[i] 
endfor 
close, 20 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 1. Run SExtractor, get the background subtraction and source detection 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
template = 'shuang_psf.sex' 
psf_sex  = name_string + '_psf.sex' 
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
cat_name = name_string + '_psf_cat.fits'
para[ index ] = 'CATALOG_NAME ' + cat_name 
;; Change the background size 
index = where( strpos( para, 'BACK_SIZE' ) ne -1 ) 
para[ index ] = 'BACK_SIZE   ' + string( bg_size ) 
;; Change the background thick 
index = where( strpos( para, 'BACKPHOTO_THICK' ) ne -1 ) 
para[ index ] = 'BACKPHOTO_THICK ' + string( bg_thick )
;; Change the background image 
index = where( strpos( para, 'CHECKIMAGE_TYPE' ) ne -1 ) 
para[ index ] = 'CHECKIMAGE_TYPE NONE'
index = where( strpos( para, 'CHECKIMAGE_NAME' ) ne -1 ) 
para[ index ] = 'CHECKIMAGE_NAME        ' 
;; Change the magnitude zeropoint
index = where( strpos( para, 'MAG_ZEROPOINT' ) ne -1 ) 
para[ index ] = 'MAG_ZEROPOINT ' + string( magzpt )
;; Change the pixel scale
index = where( strpos( para, 'PIXEL_SCALE' ) ne -1 ) 
para[ index ] = 'PIXEL_SCALE ' + string( pixel )
;; Change the saturation level
index = where( strpos( para, 'SATUR_LEVEL' ) ne -1 ) 
para[ index ] = 'SATUR_LEVEL ' + string( satu_level )
;; Change the seeing fwhm
index = where( strpos( para, 'SEEING_FWHM' ) ne -1 ) 
para[ index ] = 'SEEING_FWHM ' + string( seeing_fwhm )
;; Save the new Sex file 
openw, 20, psf_sex, width = 600 
for i = 0, ( num_lines - 1 ), 1 do begin 
    printf, 20, para[i] 
endfor 
close, 20 
;; Run SExtractor to get the background image 
spawn, 'sex  ' + image + ' -c ' + psf_sex 
if NOT file_test( cat_name ) then begin 
    message, 'Something wrong with the SExtractor !!' 
endif else begin 
    new_cat = 'psf.fits'
    spawn, 'cp  ' + cat_name + '  ' + new_cat 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 2. Make the diagnostic figure 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if NOT file_test( cat_name ) then begin 
    message, 'Can not find the Sextractor output catalog !!!' 
endif else begin 
    ;; use the half light radius and the magnitude 
    sex_cat = mrdfits( cat_name, 2 )
    rhalf    = sex_cat.FLUX_RADIUS 
    mag_best = sex_cat.MAG_BEST
    min_x = min( rhalf )
    max_x = max( rhalf )
    min_y = min( mag_best )
    max_y = max( mag_best )
    r1 = ( min_fwhm / 2.0 )
    r2 = ( max_fwhm / 2.0 )
    min_x = ( min_x > ( -0.5 ) )
    max_x = ( r2 * 4 ) < max_x
    ;; select the useful sources with enough S/N 
    snr_object = ( sex_cat.FLUX_BEST / sex_cat.FLUXERR_BEST )
    useful = where( snr_object gt snr ) 
    rhalf_useful = rhalf[ useful ]
    mag_best_useful = mag_best[ useful ]
    ;; make the plot
    psxsize = 20 
    psysize = 20
    plot0 = name_string + '_diag.ps'  
    mydevice = !D.NAME
    set_plot, 'PS'
    device, filename=plot0, font_size=8, /encapsul, $
        /color, set_font='HELVETICA BOLD', /tt_font, xsize=psxsize,$
        ysize=psysize
    xrange = [ min_x - 0.2 , max_x + 0.2 ] 
    yrange = [ max_y + 0.2 , min_y - 0.2 ]
    xtitle=textoidl('r_{half} (pixel)')
    ytitle=textoidl('magnitude best')
    ;; scatter plot of rhalf and mag_best  
    cgPlot, rhalf, mag_best, xstyle=1, ystyle=1, xrange=xrange, $
        yrange=yrange, /nodata, ytitle=ytitle, $
        /noerase, charsize=2.5, charthick=3.2, xthick=3.2, $
        ythick=3.2,  yminor=-1, xtitle=xtitle
    cgPlot, rhalf, mag_best, psym=9, symsize=1.4, $
        color=cgColor('BLACK', !D.Table_Size), /overplot
    cgPlot, rhalf_useful, mag_best_useful, psym=16, symsize=1.2, $
        color=cgColor('NAVY', !D.Table_Size), /overplot
    ;; draw the selected FWHM range
    cgPlot, [ r1, r1 ], !Y.Crange, color=cgColor('RED', !D.Table_Size), $
        linestyle=2, thick=3.0, /overplot
    cgPlot, [ r2, r2 ], !Y.Crange, color=cgColor('RED', !D.Table_Size), $
        linestyle=2, thick=3.0, /overplot
    device, /close
    set_plot, mydevice
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 3. Run PSFtractor 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
template = 'shuang.psfex' 
psf_psfex = name_string + '_psf.psfex'
if NOT file_test( template ) then begin 
    message, 'Can not find shuang.psfex file ! '
endif else begin 
    num_lines = file_lines( template )
    para = strarr( num_lines ) 
    openr, 10, template
    readf, 10, para 
    close, 10 
endelse
;; Change the psf size 
index = where( strpos( para, 'PSF_SIZE' ) ne -1 ) 
psize_string = strcompress( string( psf_size ), /remove_all )
para[ index ] = 'PSF_SIZE       ' + psize_string + ', ' + psize_string
;; Change the Polynom degress for each PSF group 
index = where( strpos( para, 'PSFVAR_DEGREES' ) ne -1 ) 
poly_string = strcompress( string( poly_degree ), /remove_all )
para[ index ] = 'PSFVAR_DEGREES  ' + poly_string
;; Change the FWHM range 
index = where( strpos( para, 'SAMPLE_FWHMRANGE' ) ne -1 ) 
fwhm_string = strcompress( string( min_fwhm ), /remove_all ) + ', ' + $
    strcompress( string( max_fwhm ), /remove_all )
para[ index ] = 'SAMPLE_FWHMRANGE    ' + fwhm_string
;; Change the Minimum S/N  
index = where( strpos( para, 'SAMPLE_MINSN' ) ne -1 ) 
minsn_string = string( snr ) 
para[ index ] = 'SAMPLE_MINSN ' + minsn_string 
;; Change the Number of Snapshots 
index = where( strpos( para, 'PSFVAR_NSNAP' ) ne -1 ) 
nsnap_string = string( nsnap ) 
para[ index ] = 'PSFVAR_NSNAP  ' + nsnap_string 
;; Change the maximum ellipticity for PSF 
index = where( strpos( para, 'SAMPLE_MAXELLIP' ) ne -1 ) 
maxellip_string = strcompress( string( maxellip ), /remove_all )
para[ index ] = 'SAMPLE_MAXELLIP    ' + maxellip_string
;; Change the check image type 
index = where( strpos( para, 'CHECKIMAGE_TYPE' ) ne -1 ) 
para[ index ] = 'CHECKIMAGE_TYPE  PROTOTYPES,SAMPLES,RESIDUALS,SNAPSHOTS_IMRES'
;; Change the check image name 
index = where( strpos( para, 'CHECKIMAGE_NAME' ) ne -1 ) 
prot_image = name_string + '_prot.fits'
samp_image = name_string + '_samp.fits'
resi_image = name_string + '_resi.fits'
snap_image = name_string + '_snap.fits'
para[ index ] = 'CHECKIMAGE_NAME ' + prot_image + ' ' + samp_image + ' ' + $
    resi_image + ' ' + snap_image 
;; Save the new PSFex file 
openw, 20, psf_psfex, width = 600 
for i = 0, ( num_lines - 1 ), 1 do begin 
    printf, 20, para[i] 
endfor 
close, 20 
;; Run PSFEx
spawn, 'psfex  ' + new_cat + ' -c ' + psf_psfex
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 4. Save the empirical psf image  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Read in the proto psf image  
snap_file = name_string + '_snap_psf.fits'
psf_ep_file = name_string + '_ep.fits'
if NOT file_test( snap_file ) then begin 
    message, 'Can not find the psfex output PROTO image ! '
endif else begin 
    psf_ep = mrdfits( snap_file, 0, head ) 
endelse 
;; Save the empirical PSF image 
mwrfits, psf_ep, psf_ep_file, head, /create 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; 5. Run ellipse on the empirical psf to get the profile 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set( ellipse ) then begin 
;the place you keep the data:
location_setting, setting
plot_setting, psetting
  datafile = setting.header
       loc = setting.workplace
    galfit = setting.galfit
    ximage = setting.ximage 
  isophote = setting.isophote
    ttools = setting.ttools
   ell_par = setting.ell_par 
imcopy_par = setting.imcopy_par 
 tdump_par = setting.tdump_par
;; set the initial parameters 
x0     = strcompress( string( ( psf_size + 1 ) / 2  ), /remove_all )
y0     = strcompress( string( ( psf_size + 1 ) / 2  ), /remove_all )
ellip0 = strcompress( string( 0.05 ), /remove_all )
pa0    = strcompress( string( 0.0 ), /remove_all )
sma0   = strcompress( string( 5.0 ), /remove_all )
minsma = strcompress( string( 0.0 ), /remove_all )
maxsma = strcompress( string( ( psf_size / 2.0 ) + 4 ), /remove_all )
step   = strcompress( string( 0.02 ), /remove_all ) 
mag0 = strcompress( string( magzpt ), /remove_all ) 
;; Run Ellipse and set the geometrical parameters to be free  
ellip_file = name_string + '_psf.ellipse'
output_bin = name_string + '_psf.bin'
output_asc = name_string + '_psf.prof'
output_head = 'head.dat' 
num_ellpara = file_lines( ell_par ) 
ell_para =strarr( num_ellpara )
openr, 10, ell_par
readf, 10, ell_para
close, 10
ell_para[0] = 'ellipse.input = "' + psf_ep_file + '"'
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
ell_para[44] = 'geompar.recenter = yes'
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
    spawn, ttools + ' tdump  @tdump.par'
    spawn, 'rm ' + output_head
endelse
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; read information from xml file 
xml = 'psfex.xml'
if NOT file_test( xml ) then begin 
    message, 'Sorry, can not find the xml file !' 
endif else begin 
    nlines = file_lines( xml ) 
    para = strarr( nlines ) 
    openr, 20, xml 
    readf, 20, para 
    close, 20 
    psf_meta = { catalog:'', image_ident:'', nexten:0, nstar_loaded_tot:0, $
        nstar_loaded_min:0, nstar_loaded_mean:0, nstar_loaded_max:0, $
        nstar_accept_tot:0, nstar_accept_min:0, nstar_accept_mean:0, $
        nstar_accept_max:0, fwhm_fluxrad_min:0, fwhm_fluxrad_mean:0, $
        fwhm_fluxrad_max:0, sampling_min:0.0, sampling_mean:0.0, $ 
        sampling_max:0.0, chi2_min:0.0, chi2_mean:0.0, chi2_max:0.0, $
        fwhm_min:0.0, fwhm_mean:0.0, fwhm_max:0.0, ellip_min:0.0, $
        ellip_mean:0.0, ellip_max:0.0, moffat_beta_min:0.0, $ 
        moffat_beta_mean:0.0, moffat_beta_max:0.0, resi_min:0.0, resi_mean:0.0,$
        resi_max:0.0, fwhm_pf_min:0.0, fwhm_pf_mean:0.0, fwhm_pf_max:0.0, $
        ellip_pf_min:0.0, ellip_pf_mean:0.0, ellip_pf_max:0.0, $
        moffat_beta_pf_min:0.0, moffat_beta_pf_mean:0.0, $
        moffat_beta_pf_max:0.0, resi_pf_min:0.0, resi_pf_mean:0.0, $
        resi_pf_max:0.0, asym_min:0.0, asym_mean:0.0, asym_max:0.0 }
    tdata = where( strpos( para, '<DATA><TABLEDATA>' ) ne -1 ) 
    index = tdata[0]
    line1 = para[ index + 2 ] 
    temp = strsplit( line1, '</TD><TD>', /extract ) 
    psf_meta.catalog = temp[1] 
    psf_meta.image_ident = temp[2] 
    psf_meta.nexten = long( temp[3] )
    line2 = para[ index + 3 ] 
    temp = strsplit( line2, '</TD><TD>', /extract ) 
    psf_meta.nstar_loaded_tot  = float( temp[1] ) 
    psf_meta.nstar_loaded_min  = float( temp[2] ) 
    psf_meta.nstar_loaded_mean = float( temp[3] ) 
    psf_meta.nstar_loaded_max  = float( temp[4] ) 
    line3 = para[ index + 4 ]
    temp = strsplit( line3, '</TD><TD>', /extract ) 
    psf_meta.nstar_accept_tot  = float( temp[1] ) 
    psf_meta.nstar_accept_min  = float( temp[2] ) 
    psf_meta.nstar_accept_mean = float( temp[3] ) 
    psf_meta.nstar_accept_max  = float( temp[4] ) 
    line4 = para[ index + 5 ] 
    temp = strsplit( line4, '</TD><TD>', /extract ) 
    psf_meta.fwhm_fluxrad_min  = float( temp[1] ) 
    psf_meta.fwhm_fluxrad_mean = float( temp[2] ) 
    psf_meta.fwhm_fluxrad_max  = float( temp[3] ) 
    line5 = para[ index + 6 ]
    temp = strsplit( line5, '</TD><TD>', /extract ) 
    psf_meta.sampling_min  = float( temp[1] ) 
    psf_meta.sampling_mean = float( temp[2] ) 
    psf_meta.sampling_max  = float( temp[3] ) 
    line6 = para[ index + 7 ] 
    temp = strsplit( line6, '</TD><TD>', /extract ) 
    psf_meta.chi2_min  = float( temp[1] ) 
    psf_meta.chi2_mean = float( temp[2] ) 
    psf_meta.chi2_max  = float( temp[3] ) 
    line7 = para[ index + 8 ]
    temp = strsplit( line7, '</TD><TD>', /extract ) 
    psf_meta.fwhm_min  = float( temp[1] ) 
    psf_meta.fwhm_mean = float( temp[2] ) 
    psf_meta.fwhm_max  = float( temp[3] ) 
    line8 = para[ index + 9 ] 
    temp = strsplit( line8, '</TD><TD>', /extract ) 
    psf_meta.ellip_min  = float( temp[1] ) 
    psf_meta.ellip_mean = float( temp[2] ) 
    psf_meta.ellip_max  = float( temp[3] ) 
    line9 = para[ index + 10 ] 
    temp = strsplit( line9, '</TD><TD>', /extract ) 
    psf_meta.moffat_beta_min  = float( temp[1] ) 
    psf_meta.moffat_beta_mean = float( temp[2] ) 
    psf_meta.moffat_beta_max  = float( temp[3] ) 
    line10 = para[ index + 11 ]
    temp = strsplit( line10, '</TD><TD>', /extract ) 
    psf_meta.resi_min  = float( temp[1] ) 
    psf_meta.resi_mean = float( temp[2] ) 
    psf_meta.resi_max  = float( temp[3] ) 
    line11 = para[ index + 12 ] 
    temp = strsplit( line11, '</TD><TD>', /extract ) 
    psf_meta.fwhm_pf_min  = float( temp[1] ) 
    psf_meta.fwhm_pf_mean = float( temp[2] ) 
    psf_meta.fwhm_pf_max  = float( temp[3] ) 
    line12 = para[ index + 13 ] 
    temp = strsplit( line12, '</TD><TD>', /extract ) 
    psf_meta.ellip_pf_min  = float( temp[1] ) 
    psf_meta.ellip_pf_mean = float( temp[2] ) 
    psf_meta.ellip_pf_max  = float( temp[3] ) 
    line13 = para[ index + 14 ] 
    temp = strsplit( line13, '</TD><TD>', /extract ) 
    psf_meta.moffat_beta_pf_min  = float( temp[1] ) 
    psf_meta.moffat_beta_pf_mean = float( temp[2] ) 
    psf_meta.moffat_beta_pf_max  = float( temp[3] ) 
    line14 = para[ index + 15 ]
    temp = strsplit( line14, '</TD><TD>', /extract ) 
    psf_meta.resi_pf_min  = float( temp[1] ) 
    psf_meta.resi_pf_mean = float( temp[2] ) 
    psf_meta.resi_pf_max  = float( temp[3] ) 
    line15 = para[ index + 16 ] 
    temp = strsplit( line15, '</TD><TD>', /extract ) 
    psf_meta.asym_min  = float( temp[1] ) 
    psf_meta.asym_mean = float( temp[2] ) 
    psf_meta.asym_max  = float( temp[3] ) 

    psf_sum = psf_meta
endelse 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
end
