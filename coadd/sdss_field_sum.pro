pro sdss_field_sum, logfile

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; exam the properties of SDSS images within certain region  
;;; written by Song Huang 03.28.2012 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
on_error, 2
compile_opt idl2

pix = 0.396
naxis1 = 2048 
naxis2 = 1489
;; the sdss field summary file 
sumfile = 'sdss_field_sum.fit'
if NOT file_test( sumfile ) then begin 
    message, 'Can not find the SDSS fields summary file ! '
endif else begin 
    cat = mrdfits( sumfile, 1 ) 
endelse
;; read in the input log file 
if NOT file_test( logfile ) then begin 
    message, 'Can not find the input log file !' 
endif else begin 
    readcol, logfile, run, camcol, field, ra, dec, format='I,I,I,F,F', $
        comment='#', delimiter=' ', skipline=4  
    nimage = n_elements( run ) 
    print, '#############################################' 
    print, 'INPUT LOG  : ' + logfile
    print, '##############################################'
    nlines = file_lines( logfile ) 
    temp = strarr( nlines ) 
    openr, 10, logfile 
    readf, 10, temp 
    close, 10 
    temp1 = strsplit( strcompress( temp[0], /remove_all ), '=', /extract )
    ra_cen = float( temp1[1] )
    temp1= strsplit( strcompress( temp[1], /remove_all ), '=', /extract )
    dec_cen = float( temp1[1] )
    temp1 = strsplit( strcompress( temp[2], /remove_all ), '=', /extract )
    fov = float( temp1[1] )
    fov_pix = long( fov * 60.0 * 60.0 / pix )
    temp1 = strsplit( strcompress( temp[3], /remove_all ), '=', /extract )
    filter = temp1[1]

    ;;;;;;;;
    fov = fov * 2.0 
    ;;;;;;;;

    print, 'RA  CENTER  : ' + string( ra_cen )
    print, 'DEC CENTER  : ' + string( dec_cen ) 
    print, 'FOV         : ' + string( fov ) 
    print, 'FILTER      : ' + filter
    print, '##############################################'
    print, 'Number of input fields : ' + string( nimage ) 

    logfile = strcompress( logfile, /remove_all ) 
    temp1 = strsplit( logfile, '.', /extract ) 
    nseg = n_elements( temp1 ) - 1 
    seg = temp1[0]
    for m = 1, ( nseg - 1 ), 1 do begin 
        seg = seg + '.' + temp1[m] 
    endfor
    log_string = seg
endelse
;; find the field information from the summary file 
found_arr = lonarr( nimage )
for i = 0, ( nimage - 1 ), 1 do begin 
    aa = long( run[i] )
    bb = long( camcol[i] ) 
    cc = long( field[i] ) 
    aa2 = long( cat.run ) 
    bb2 = long( cat.camcol ) 
    cc2 = long( cat.field ) 
    found = where( ( aa2 eq aa ) and ( bb2 eq bb ) and ( cc2 eq cc ) ) 
    if ( found[0] eq -1 ) then begin 
        print, 'IMAGE Number : ' + string( i + 1 ) 
        print, aa, bb, cc
        print, ' IMAGE NOT FOUND !!!!!!' 
        found_arr[i] = -1 
    endif else begin 
        ;print, ' IMAGE FOUND !  INDEX = ' + string( found[0] ) 
        found_arr[i] = found[0] 
    endelse
endfor
index = found_arr[ where( found_arr ne -1 ) ]
;index = ( where( cat.run eq long( run ) ) and $
;    where( cat.camcol eq long( camcol ) ) and $
;    where( cat.field eq long( field ) ) )
nfound = n_elements( index ) 
print, 'Number of found images : ' + string( nfound )
if ( index[0] eq -1 ) then begin 
    message, 'The input log has no useful images ! Something Wrong ! '
endif else begin 
    if ( nfound ne nimage ) then begin 
        print, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
        print, ' WARNNING : There are Missing Images !!!!'
        print, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
    endif 
endelse
;; read useful information from summary file 
run = cat[ index ].run 
camcol = cat[ index ].camcol 
field = cat[ index ].field
raMin = cat[ index ].raMin 
raMax = cat[ index ].raMax 
decMin = cat[ index ].decMin 
decMax = cat[ index ].decMax
am_ra = cat[ index ].ra 
am_dec= cat[ index ].dec
am_a = cat[ index ].a 
am_b = cat[ index ].b
am_c = cat[ index ].c
am_d = cat[ index ].d 
am_e = cat[ index ].e 
am_f = cat[ index ].f
;; psf information 
case filter of 
    'u': begin 
            psfwidth = cat[ index ].psfWidth_u
            psfSigma1 = cat[ index ].psfSigma1_u 
            psfSigma2 = cat[ index ].psfSigma2_u 
            airmass = cat[ index ].airmass_u 
            calibStatus = long( cat[ index ].calibStatus_u )
            imageStatus = long( cat[ index ].imageStatus_u )
         end 
    'g': begin 
            psfwidth = cat[ index ].psfWidth_g
            psfSigma1 = cat[ index ].psfSigma1_g 
            psfSigma2 = cat[ index ].psfSigma2_g 
            airmass = cat[ index ].airmass_g 
            calibStatus = long( cat[ index ].calibStatus_g )
            imageStatus = ( cat[ index ].imageStatus_g )
         end 
    'r': begin 
            psfwidth = cat[ index ].psfWidth_r
            psfSigma1 = cat[ index ].psfSigma1_r 
            psfSigma2 = cat[ index ].psfSigma2_r 
            airmass = cat[ index ].airmass_r 
            calibStatus = long( cat[ index ].calibStatus_r )
            imageStatus = long( cat[ index ].imageStatus_r )
         end 
endcase 

;; calculate the ra and dec at the four corners 
 ra_0 = am_a + 0.0 * am_b + 0.00 * am_c 
dec_0 = am_d + 0.0 * am_e + 0.00 * am_f 
 ra_1 = am_a + naxis1 * am_b + 0.00 * am_c 
dec_1 = am_d + naxis1 * am_e + 0.00 * am_f 
 ra_2 = am_a + naxis1 * am_b + naxis2 * am_c 
dec_2 = am_d + naxis1 * am_e + naxis2 * am_f 
 ra_3 = am_a + 0.0 * am_b + naxis2 * am_c 
dec_3 = am_d + 0.0 * am_e + naxis2 * am_f 

good_bad = make_array( nfound, VALUE=1 )
;; get information about these images 
;; psfwidth
min_psfwidth = min( psfwidth ) 
max_psfwidth = max( psfwidth )
resistant_mean, psfwidth, 3, mean_psfwidth, sm_psfwidth, nreject 
sigma_psfwidth = sm_psfwidth * sqrt( nfound - 1) 
sorted = psfwidth[ SORT( psfwidth ) ]
if ( ( n_elements( sorted ) mod 2 ) eq 0 ) then begin 
    ind = n_elements( sorted ) / 2 
    med_psfwidth = ( sorted[ ind - 1 ] + sorted[ ind ] ) / 2.0 
    lowergroup = sorted[0:(ind-1)]
    highergroup = sorted[ind:(n_elements(psfwidth)-1)] 
endif else begin 
    ind = long( n_elements( sorted ) / 2 )
    med_psfwidth = sorted[ ind ] 
    lowergroup = sorted[0:(ind-1)] 
    highergroup = sorted[(ind+1):(n_elements(psfwidth)-1)]
endelse
lq_psfwidth = median( lowergroup, /EVEN ) 
hq_psfwidth = median( highergroup, /EVEN )
irq_psfwidth = ( hq_psfwidth - lq_psfwidth )
ifen1 = ( lq_psfwidth - 1.5 * irq_psfwidth )
ifen2 = ( hq_psfwidth + 1.5 * irq_psfwidth )
ofen1 = ( lq_psfwidth - 3.0 * irq_psfwidth )
ofen2 = ( hq_psfwidth + 3.0 * irq_psfwidth )
index_largepsf = where( psfwidth ge 2.0 )
if ( index_largepsf[0] ne -1 ) then begin 
    good_bad[ index_largepsf ] = 0
    n_largepsf = n_elements( index_largepsf )
endif else begin 
    n_largepsf = 0 
endelse
print, '###################################################'
print, 'PARAMETERS :    PSF_WIDTH     '
print, 'MIN / MAX  :  ', min_psfwidth, max_psfwidth 
print, 'MEAN       :  ', mean_psfwidth
print, 'SIGMA      :  ', sigma_psfwidth 
print, 'N_REJECT   :  ', nreject 
print, 'MEDIAN     :  ', med_psfwidth 
print, 'Low/High Q :  ', lq_psfwidth, hq_psfwidth 
print, 'IRQ        :  ', irq_psfwidth
print, 'IFENCE 1/2 :  ', ifen1, ifen2
print, 'OFENCE 1/2 :  ', ofen1, ofen2
print, 'N_LARGE_PSF:  ', n_largepsf 

;; psfsigma1
min_psfsigma1 = min( psfsigma1 ) 
max_psfsigma1 = max( psfsigma1 )
resistant_mean, psfsigma1, 3, mean_psfsigma1, sm_psfsigma1, nreject 
sigma_psfsigma1 = sm_psfsigma1 * sqrt( nfound - 1) 
sorted = psfsigma1[ SORT( psfsigma1 ) ]
if ( ( n_elements( sorted ) mod 2 ) eq 0 ) then begin 
    ind = n_elements( sorted ) / 2 
    med_psfsigma1 = ( sorted[ ind - 1 ] + sorted[ ind ] ) / 2.0 
    lowergroup = sorted[0:(ind-1)]
    highergroup = sorted[ind:(n_elements(psfsigma1)-1)] 
endif else begin 
    ind = long( n_elements( sorted ) / 2 )
    med_psfsigma1 = sorted[ ind ] 
    lowergroup = sorted[0:(ind-1)] 
    highergroup = sorted[(ind+1):(n_elements(psfsigma1)-1)]
endelse
lq_psfsigma1 = median( lowergroup, /EVEN ) 
hq_psfsigma1 = median( highergroup, /EVEN )
irq_psfsigma1 = ( hq_psfsigma1 - lq_psfsigma1 )
ifen1 = ( lq_psfsigma1 - 1.5 * irq_psfsigma1 )
ifen2 = ( hq_psfsigma1 + 1.5 * irq_psfsigma1 )
ofen1 = ( lq_psfsigma1 - 3.0 * irq_psfsigma1 )
ofen2 = ( hq_psfsigma1 + 3.0 * irq_psfsigma1 )
print, '###################################################'
print, 'PARAMETERS :    PSF_SIGMA_1   '
print, 'MIN / MAX  :  ', min_psfsigma1, max_psfsigma1 
print, 'MEAN       :  ', mean_psfsigma1
print, 'SIGMA      :  ', sigma_psfsigma1 
print, 'N_REJECT   :  ', nreject 
print, 'MEDIAN     :  ', med_psfsigma1 
print, 'Low/High Q :  ', lq_psfsigma1, hq_psfsigma1 
print, 'IRQ        :  ', irq_psfsigma1
print, 'IFENCE 1/2 :  ', ifen1, ifen2
print, 'OFENCE 1/2 :  ', ofen1, ofen2

;; psfsigma2
min_psfsigma2 = min( psfsigma2 ) 
max_psfsigma2 = max( psfsigma2 )
resistant_mean, psfsigma2, 3, mean_psfsigma2, sm_psfsigma2, nreject 
sigma_psfsigma2 = sm_psfsigma2 * sqrt( nfound - 1) 
sorted = psfsigma2[ SORT( psfsigma2 ) ]
if ( ( n_elements( sorted ) mod 2 ) eq 0 ) then begin 
    ind = n_elements( sorted ) / 2 
    med_psfsigma2 = ( sorted[ ind - 1 ] + sorted[ ind ] ) / 2.0 
    lowergroup = sorted[0:(ind-1)]
    highergroup = sorted[ind:(n_elements(psfsigma2)-1)] 
endif else begin 
    ind = long( n_elements( sorted ) / 2 )
    med_psfsigma2 = sorted[ ind ] 
    lowergroup = sorted[0:(ind-1)] 
    highergroup = sorted[(ind+1):(n_elements(psfsigma2)-1)]
endelse
lq_psfsigma2 = median( lowergroup, /EVEN ) 
hq_psfsigma2 = median( highergroup, /EVEN )
irq_psfsigma2 = ( hq_psfsigma2 - lq_psfsigma2 )
ifen1 = ( lq_psfsigma2 - 1.5 * irq_psfsigma2 )
ifen2 = ( hq_psfsigma2 + 1.5 * irq_psfsigma2 )
ofen1 = ( lq_psfsigma2 - 3.0 * irq_psfsigma2 )
ofen2 = ( hq_psfsigma2 + 3.0 * irq_psfsigma2 )
print, '###################################################'
print, 'PARAMETERS :    PSF_SIGMA_2   '
print, 'MIN / MAX  :  ', min_psfsigma2, max_psfsigma2 
print, 'MEAN       :  ', mean_psfsigma2
print, 'SIGMA      :  ', sigma_psfsigma2 
print, 'N_REJECT   :  ', nreject 
print, 'MEDIAN     :  ', med_psfsigma2 
print, 'Low/High Q :  ', lq_psfsigma2, hq_psfsigma2 
print, 'IRQ        :  ', irq_psfsigma2
print, 'IFENCE 1/2 :  ', ifen1, ifen2
print, 'OFENCE 1/2 :  ', ofen1, ofen2

;; airmass
min_airmass = min( airmass ) 
max_airmass = max( airmass )
resistant_mean, airmass, 3, mean_airmass, sm_airmass, nreject 
sigma_airmass = sm_airmass * sqrt( nfound - 1) 
sorted = airmass[ SORT( airmass ) ]
if ( ( n_elements( sorted ) mod 2 ) eq 0 ) then begin 
    ind = n_elements( sorted ) / 2 
    med_airmass = ( sorted[ ind - 1 ] + sorted[ ind ] ) / 2.0 
    lowergroup = sorted[0:(ind-1)]
    highergroup = sorted[ind:(n_elements(airmass)-1)] 
endif else begin 
    ind = long( n_elements( sorted ) / 2 )
    med_airmass = sorted[ ind ] 
    lowergroup = sorted[0:(ind-1)] 
    highergroup = sorted[(ind+1):(n_elements(airmass)-1)]
endelse
lq_airmass = median( lowergroup, /EVEN ) 
hq_airmass = median( highergroup, /EVEN )
irq_airmass = ( hq_airmass - lq_airmass )
ifen1 = ( lq_airmass - 1.5 * irq_airmass )
ifen2 = ( hq_airmass + 1.5 * irq_airmass )
ofen1 = ( lq_airmass - 3.0 * irq_airmass )
ofen2 = ( hq_airmass + 3.0 * irq_airmass )
print, '###################################################'
print, 'PARAMETERS :    AIRMASS  '
print, 'MIN / MAX  :  ', min_airmass, max_airmass 
print, 'MEAN       :  ', mean_airmass
print, 'SIGMA      :  ', sigma_airmass 
print, 'N_REJECT   :  ', nreject 
print, 'MEDIAN     :  ', med_airmass 
print, 'Low/High Q :  ', lq_airmass, hq_airmass 
print, 'IRQ        :  ', irq_airmass
print, 'IFENCE 1/2 :  ', ifen1, ifen2
print, 'OFENCE 1/2 :  ', ofen1, ofen2

;; calibStatus
calibStatus = strcompress( string( calibStatus ), /remove_all )
index_photo = where( calibStatus eq '1' ) 
index_nonphotot = where( calibStatus ne '1' ) 
nphoto = n_elements( index_photo ) 
if ( index_nonphotot[0] eq -1 ) then begin 
    n_nonphoto = 0 
endif else begin 
    good_bad[ index_nonphotot ] = 0
    n_nonphoto = n_elements( index_nonphotot )
endelse
print, '##################################################' 
print, 'N_PHOTO_NIGHT : ', nphoto 
print, 'N_NON_PHOTO   : ', n_nonphoto

;; imageStatus 
imageStatus = strcompress( string( imageStatus ), /remove_all )
index_clear = where( imageStatus eq '1' ) 
index_nonclear = where( imageStatus ne '1' ) 
nclear = n_elements( index_clear ) 
if ( index_nonclear[0] eq -1 ) then begin 
    n_nonclear = 0 
endif else begin 
    good_bad[ index_nonclear ] = 0
    n_nonclear = n_elements( index_nonclear )
endelse
print, '##################################################' 
print, 'N_CLEAR_NIGHT : ', nclear 
print, 'N_NON_CLEAR   : ', n_nonclear

;; get statistics of bad images 
index_bad = where( good_bad eq 0 )
nbad = n_elements( index_bad )
print, '##################################################'
print, 'N_BAD_IMAGE   : ', nbad 
print, ' RUN     CAMCOL    FIELD    PSFWIDTH    CALIB_STATUS    IMAGE_STATUS'
if ( index_bad[0] ne -1 ) then begin 
    for i = 0, ( nbad - 1 ), 1 do begin 
        b = index_bad[i]
        print, run[b], camcol[b], field[b], psfwidth[b], calibStatus[b], $
            imageStatus[b] 
    endfor
endif
print, '##################################################'

;; get statistics of good images 
index_good = where( good_bad eq 1 )
if ( index_good[0] eq -1 ) then begin 
    print, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print, '!!!!!!!!!!!!! WARNING : NO USEFUL IMAGE FOUND !!!!!!!!!!!!!'
    print, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
endif else begin 
    ngood = n_elements( index_good )
    fgood = ( ngood / nfound )
    print, '##################################################'
    print, 'N_GOOD_IMAGE  : ', ngood
    print, 'FRAC_GOOD     : ', fgood 
    ;; save the good images list 
    good_lis = log_string + '.lis' 
    openw, 20, good_lis, width = 400 
    printf, 20, temp[0]
    printf, 20, temp[1]
    printf, 20, temp[2]
    printf, 20, temp[3]
    camcol = long( camcol )
    printf, 20, '#RUN   CAMCOL    FIELD    RA_MIN   DEC_MIN    RA_MAX    DEC_MAX '
    for j = 0, ( ngood - 1 ), 1 do begin 
        g = index_good[j]
        printf, 20, string( run[g] ) + '  ' + string( camcol[g] ) + '  ' + $ 
            string( field[g] ) + '  ' + string( raMin[g] ) + '  ' + $ 
            string( decMin[g] ) + '  ' + string( raMax[g] ) + '  ' + $
            string( decMax[g] ) 
    endfor
    close, 20 
endelse
;; make overlap plot 
overlap = log_string + '_over.ps' 
min_0 = min( [ min( ra_0 ), min( ra_1 ), min( ra_2 ), min( ra_3 ) ] ) - 0.1 
max_0 = max( [ max( ra_0 ), max( ra_1 ), max( ra_2 ), max( ra_3 ) ] ) + 0.1 
min_1 = min( [ min( dec_0 ), min( dec_1 ), min( dec_2 ), min( dec_3 ) ] ) - 0.1 
max_1 = max( [ max( dec_0 ), max( dec_1 ), max( dec_2 ), max( dec_3 ) ] ) + 0.1 
print, min_0, max_0, min_1, max_1
x0 = ( ra_cen - fov / 2.0 )
x1 = ( ra_cen + fov / 2.0 )
y0 = ( dec_cen - fov / 2.0 )
y1 = ( dec_cen + fov / 2.0 )
xcen = ra_cen 
ycen = dec_cen 
min_ra = ( x0 - 0.1 ) 
max_ra = ( x1 + 0.1 ) 
min_dec = ( y0 - 0.1 ) 
max_dec = ( y1 + 0.1 ) 
xrange = [ min_ra, max_ra ] 
yrange = [ min_dec, max_dec ]

;yx_ratio = ( ( max_dec - min_dec ) / ( max_ra - min_ra ) ) 
;if ( yx_ratio gt 1 ) then begin 
;    ysize = 24 
;    xsize = ( 24 / yx_ratio ) 
;endif else begin 
;    xsize = 24
;    ysize = ( xsize * yx_ratio ) 
;endelse
xsize = 24 
ysize = 24

xtitle = 'RA (degree)' 
ytitle = 'DEC (degree) '
if ( fov le 0.2 ) then begin 
    xtick = 0.1
    ytick = 0.05 
endif else begin
    if ( fov ge 0.55 ) then begin 
        xtick = 0.4
        ytick = 0.2 
    endif else begin 
        xtick = 0.2 
        ytick = 0.1
    endelse 
endelse

mydevice = !D.NAME
set_plot, 'PS'
device, filename=overlap, font_size=8, /encapsul, $
  /color, set_font='HELVETICA BOLD', /tt_font, xsize=xsize, ysize=ysize
;; original image
cgPlot, [ 0.0, 1.0 ], [ 0.0, 1.0 ], xstyle=1, ystyle=1, $
    xrange=xrange, yrange=yrange, xtitle=xtitle, ytitle=ytitle, $
    xtickinterval=xtick, ytickinterval=ytick, charsize=2.8, charthick=3.8, $
    xthick=3.8, ythick=3.8, /noerase, /nodata

jpeg = log_string + '.jpg' 
if file_test( jpeg ) then begin 
    print, '!!!!!!!! FOUND THE JPG FILE !!!!!!!!'
    read_jpeg, jpeg, img_jpeg 
    scale_jpeg = bytscl( img_jpeg )
    xx = ra_cen - ( 1021 * pix / 3600.0 )
    xy = ra_cen + ( 1021 * pix / 3600.0 )
    yx = dec_cen - ( 1021 * pix / 3600.0 )
    yy = dec_cen + ( 1021 * pix / 3600.0 )
    oplotimage, scale_jpeg, imgxrange=[xx,xy], imgyrange=[yx,yy]
endif

;; overplot the bad images region 
if ( index_bad[0] ne -1 ) then begin 
    for i = 0, ( nbad - 1 ), 1 do begin 
        b = index_bad[i]
        ra0 = raMin[b]
        ra1 = raMax[b]
        dec0 = decMin[b] 
        dec1 = decMax[b]
        cgPlots, [ra0,ra1,ra1,ra0,ra0], [dec0,dec0,dec1,dec1,dec0], $ 
            linestyle=1, thick=1.8, color=cgColor( 'RED', !D.Table_Size ), /data 
    endfor
endif
;; overplot the good image region
if ( index_good[0] ne -1 ) then begin 
    for i = 0, ( ngood - 1 ), 1 do begin 
        b = index_good[i]
        ra0 = raMin[b]
        ra1 = raMax[b]
        dec0 = decMin[b] 
        dec1 = decMax[b]
        cgPlots, [ra0,ra1,ra1,ra0,ra0], [dec0,dec0,dec1,dec1,dec0], $
            linestyle=0, thick=1.8, color=cgColor( 'BLUE', !D.Table_Size ), $
            /data 
    endfor
endif

cgPlots, [ xcen ], [ ycen ], psym=46, symsize=3.0, thick=4.0, $ 
    color=cgColor( 'ORG5', !D.Table_Size ), /data
cgPlots, [x0,x1,x1,x0,x0], [y0,y0,y1,y1,y0], $
    color=cgColor( 'BLACK', !D.Table_Size ), linestyle=0, thick=5.0, /data


device, /close
set_plot, mydevice

;; histogram plots 

end
