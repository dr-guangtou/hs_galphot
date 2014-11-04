pro get_sdss_image, list, filter 

wget = '/usr/bin/wget'
;wget = '/usr/local/homebrew/bin/wget'
bunzip2 = '/bin/bunzip2'
;bunzip2 = '/opt/local/bin/bunzip2'

band_array = [ 'u', 'g', 'r', 'i', 'z' ]
if ( where( band_array eq filter ) eq -1 ) then begin 
    message, 'The filter should be among u, g, r, i, z'
endif else begin 
    filter = strcompress( filter, /remove_all )
endelse 

if file_test( list ) then begin 
    readcol, list, run, camcol, field, ra, dec, format='I,I,I,F,F', $
        delimiter=' ', comment='#', skipline=4 
    n_img = n_elements( run ) 
    initial = 'http://data.sdss3.org/sas/dr8/groups/boss/photoObj/frames/301/' 
    for i = 0, ( n_img - 1 ), 1 do begin 
        para1 = run[i] 
        para2 = camcol[i] 
        para3 = field[i] 
    
        para1_s1 = strcompress( string( para1 ), /remove_all )
        para2_s1 = strcompress( string( para2 ), /remove_all )
        para3_s1 = strcompress( string( para3 ), /remove_all )
        
        if ( para1 lt 100 ) then begin 
            para1_s2 = '0000' + para1_s1
        endif else begin 
            if ( ( para1 ge 100 ) and ( para1 lt 1000 ) ) then begin 
                para1_s2 = '000' + para1_s1 
            endif else begin 
                para1_s2 = '00' + para1_s1 
            endelse 
        endelse
    
        if ( para3 lt 10 ) then begin 
            para3_s2 = '000' + para3_s1 
        endif
        if ( ( para3 ge 10 ) and ( para3 lt 100 ) ) then begin 
            para3_s2 = '00' + para3_s1 
        endif 
        if ( ( para3 ge 100 ) and ( para3 lt 1000 ) ) then begin 
            para3_s2 = '0' + para3_s1 
        endif
        print, 'IMAGE : ' + para1_s2 + '  ' + para2_s1 + '  ' + para3_s2
    
        img_loc = initial + para1_s1 + '/' + para2_s1 + '/'  
        img_file = 'frame-' + filter + '-' + para1_s2 + '-' + para2_s1 + '-' + $
            para3_s2 + '.fits.bz2' 
        
        ;; download ?
            spawn, wget + ' ' + img_loc + img_file
            if file_test( img_file ) then begin 
                spawn, bunzip2 + ' ' + img_file
            endif 
    endfor
endif

end
