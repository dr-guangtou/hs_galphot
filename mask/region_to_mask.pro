pro region_to_mask, input_image, reg_file, old_mask, savejpeg = savejpeg, overlap = overlap
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;  Add a list of square regions to the mask image:
;;  Written by Song Huang, 2011.12.08
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
on_error, 2
compile_opt idl2

if N_params() lt 2 then begin 
    print,  'Syntax - region_to_mask, input_image, reg_file, old_mask'
    print,  '       /savejpeg, /overlap' 
    return
  endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
input_image = strcompress( input_image, /remove_all )
if NOT file_test( input_image ) then begin 
    message, 'Can not find the input image !!!'
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
old_mask = strcompress( old_mask, /remove_all )
if NOT file_test( old_mask ) then begin 
    mark_old = 0 
endif else begin 
    mark_old = 1 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if NOT file_test( strcompress( reg_file, /remove_all ) ) then begin 
  message, 'Can not find the region file !!!'
endif else begin 
  lines = file_lines( reg_file )
  regions = strarr( lines ) 

  if ( lines le 4 ) then begin 
    message, 'No Usefull Region in the Reg_File !!'
  endif else begin 
    openr, 10, reg_file 
    readf, 10, regions 
    close, 10
    if ( strcompress( regions[3], /remove_all ) ne 'image' ) then begin 
      message, 'The Coordinates for Regions should be in Image Unit !!'
    endif else begin 
      regions2 = strarr( lines - 4 )
      regions2 = regions[ 4 : ( lines - 1 ) ] 
      regions = regions2
      lines = lines - 4
    endelse
  endelse
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
img_ori = mrdfits( input_image, 0, head_ori )
size_ori = size( img_ori, /dimension )
naxis1 = size_ori[0]
naxis2 = size_ori[1]
if ( mark_old eq 0 ) then begin 
    img_msk = fltarr( naxis1, naxis2 )
    header_msk = head_ori
endif else begin 
    img_msk = mrdfits( old_mask, 0, head_msk )
    header_msk = head_msk 
endelse
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
n_box = 0
n_circle = 0
n_ellipse = 0
n_polygon = 0
n_line = 0
n_vector = 0
n_annulus = 0
n_dot = 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
for i = 0, ( lines - 1 ), 1 do begin 

  region = regions[ i ]

  temp = strsplit( region, '(), #', /extract )
  print, '############# region number: ' + string( i+1 ) + ' ##############'
  print, temp
  number = n_elements( temp ) - 1
  region_type = strcompress( temp[0], /remove_all )
  case region_type of 
    'box'    :  n_box += 1
    'circle' :  n_circle += 1
    'ellipse':  n_ellipse += 1
    'polygon':  n_polygon += 1
    'line'   :  n_line += 1
    'vector' :  n_vector += 1
    'ruler'  :  n_line += 1
    'projection' : n_line += 1
    'annulus' : n_annulus += 1
    'point'  :  n_dot += 1
  else: print, '## Only Box, Circle, Ellipse, Polygon, Vector, Ruler' + $
    ', Point, Annulus is Allowed !!'
  endcase
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; For Box Region: 
  if ( region_type eq 'box' ) then begin 
    if ( number eq 4 ) then begin 
      xc = temp[ 1 ]
      yc = temp[ 2 ]
      sx = temp[ 3 ]
      sy = temp[ 4 ]
      sx2 = 0
      sy2 = 0
      angle = 0.0
    endif else begin 
      xc = temp[ 1 ]
      yc = temp[ 2 ]
      if ( number eq 5 ) then begin 
        sx = temp[ 3 ]
        sy = temp[ 4 ]
        sx2 = 0
        sy2 = 0
        angle = ( 180.0 - float( temp[ 5 ] ) ) * !DtoR
      endif else begin 
        sx = temp[ 3 ]
        sy = temp[ 4 ]
        sx2 = temp[ 5 ]
        sy2 = temp[ 6 ]
        if ( ( number / 2 ) ne 0 ) then begin 
            angle = ( 180.0 - float( temp[ 7 ] ) ) * !DtoR
        endif else begin 
            angle = 0.0 * !DtoR
        endelse
      endelse
    endelse

      x1 = ( sx / 2.0 ) 
      y1 = ( sy / 2.0 )
      x2 = ( sx / 2.0 )
      y2 = ( -1.0 * sy / 2.0 )
      x3 = ( -1.0 * sx / 2.0 )
      y3 = ( -1.0 * sy / 2.0 )
      x4 = ( -1.0 * sx / 2.0 )
      y4 = ( sy / 2.0 )
      xa = [ x1, x2, x3, x4 ]
      ya = [ y1, y2, y3, y4 ]
      xb = xa
      yb = ya
      for j = 0, 3, 1 do begin 
        xb[ j ] = xa[ j ] * cos( angle ) + ya[ j ] * sin( angle ) + xc
        yb[ j ] = ya[ j ] * cos( angle ) - xa[ j ] * sin( angle ) + yc
      endfor

      x1 = ( sx2 / 2.0 ) 
      y1 = ( sy2 / 2.0 )
      x2 = ( sx2 / 2.0 )
      y2 = ( -1.0 * sy2 / 2.0 )
      x3 = ( -1.0 * sx2 / 2.0 )
      y3 = ( -1.0 * sy2 / 2.0 )
      x4 = ( -1.0 * sx2 / 2.0 )
      y4 = ( sy2 / 2.0 )
      xa = [ x1, x2, x3, x4 ]
      ya = [ y1, y2, y3, y4 ]
      xd = xa
      yd = ya
      for j = 0, 3, 1 do begin 
        xd[ j ] = xa[ j ] * cos( angle ) + ya[ j ] * sin( angle ) + xc
        yd[ j ] = ya[ j ] * cos( angle ) - xa[ j ] * sin( angle ) + yc
      endfor
     
      if ( ( where( xb le 0 ) )[0] ne -1 ) then begin 
          xb[ where( xb le 0 ) ] = 0
      endif
      if ( ( where( yb le 0 ) )[0] ne -1 ) then begin       
          yb[ where( yb le 0 ) ] = 0
      endif
      if ( ( where( xb ge ( naxis1 - 1 ) ) )[0] ne -1 ) then begin       
          xb[ where( xb ge ( naxis1 - 1 ) ) ] = ( naxis1 - 1 )
      endif
      if ( ( where( yb ge ( naxis2 - 1 ) ) )[0] ne -1 ) then begin       
          yb[ where( yb ge ( naxis2 - 1 ) ) ] = ( naxis2 - 1 ) 
        endif

      if ( ( sx2 eq 0 ) and ( sy2 eq 0 ) ) then begin 
          index_box = polyfillv( xb, yb, naxis1, naxis2 )
          img_ori[ index_box ] = -9999.9
          img_msk[ index_box ] = 1
      endif else begin 
          index_box1= polyfillv( xb, yb, naxis1, naxis2 )
          index_box2 = polyfillv( xd, yd, naxis1, naxis2 )
          img_ori[ index_box2 ] = -9999.9
          img_msk[ index_box2 ] = 1
          img_ori[ index_box1 ] = img_cor[ index_box1 ]
          img_msk[ index_box1 ] = 0
      endelse
    endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; For Circle Region: 
  if ( region_type eq 'circle' ) then begin
          xc = temp[ 1 ]
          yc = temp[ 2 ]
          radius = temp[ 3 ]

          dist_circle, circle_mask, [ naxis1, naxis2 ], xc, yc, /double
          index_circle = where( circle_mask le radius )
          img_ori[ index_circle ] = -9999.9
          img_msk[ index_circle ] = 1
        endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; For Ellipse Region: 
  if ( region_type eq 'ellipse' ) then begin
          xc = temp[ 1 ]
          yc = temp[ 2 ]
          ra = float( temp[ 3 ] )
          rb = float( temp[ 4 ] )
          if ( number gt 4 ) then begin 
              if ( number lt 7 ) then begin 
                  angle = float( temp[ 5 ] ) + 90.0
                  rc = 1
                  rd = 1
              endif else begin 
                  angle = float( temp[ number ] ) + 90.0
                  rd = float( temp[ number - 1 ] )
                  rc = float( temp[ number - 2 ] )
              endelse
          endif else begin 
            angle = 0.0
            rc = 1
            rd = 1
          endelse
          ratio = ra / rb

          dist_ellipse, ellipse_mask, [ naxis1, naxis2 ], xc, yc, ratio, $
            angle, /double 
          if ( rc eq 1 ) then begin 
              index_ellipse = where( ellipse_mask le ra )
              img_ori[ index_ellipse ] = -9999.9
              img_msk[ index_ellipse ] = 1
          endif else begin 
              index_ellipse1 = where( ellipse_mask le ra )
              index_ellipse2 = where( ellipse_mask le rc )
              img_ori[ index_ellipse2 ] = -9999.9
              img_msk[ index_ellipse2 ] = 1
              img_ori[ index_ellipse1 ] = img_cor[ index_ellipse1 ]
              img_msk[ index_ellipse1 ] = 0
          endelse           
  endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; For Polygon Region: 
  if ( region_type eq 'polygon' ) then begin
        n_points = long( number / 2 )
        xx = fltarr( n_points )
        yy = fltarr( n_points )

        for j = 1, ( n_points ), 1 do begin 
             xx[ j-1 ] = float( temp[ 2 * j - 1 ] )
             yy[ j-1 ] = float( temp[ 2 * j ] )
           endfor

        index_poly = polyfillv( xx, yy, naxis1, naxis2 )
        img_ori[ index_poly ] = -9999.9
        img_msk[ index_poly ] = 1
      endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; For Line and Ruler Region: 
  if ( ( region_type eq 'line' ) or ( region_type eq 'ruler' ) ) then begin
        x1 = temp[ 1 ]
        y1 = temp[ 2 ]
        x2 = temp[ 3 ]
        y2 = temp[ 4 ]
        xx = [ x1, x2, ( x2 + 3 ), ( x1 - 3 ) ]
        yy = [ y1, y2, ( y2 + 3 ), ( y1 + 3 ) ]

        index_line = polyfillv( xx, yy, naxis1, naxis2 )
        img_ori[ index_line ] = -9999.9
        img_msk[ index_line ] = 1
      endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; For Vector Region: 
  if ( region_type eq 'vector' ) then begin
        x1 = temp[ 1 ]
        y1 = temp[ 2 ]
        length = float( temp[ 3 ] )
        angle = float( temp[ 4 ] ) * !DtoR
        x2 = x1 + length * cos( angle )
        y2 = y1 + length * sin( angle )
        xx = [ x1, x2, ( x2 + 4 ), ( x1 - 4 ) ]
        yy = [ y1, y2, ( y2 + 4 ), ( y1 + 4 ) ]

        index_vector = polyfillv( xx, yy, naxis1, naxis2 )
        img_ori[ index_vector ] = -9999.9
        img_msk[ index_vector ] = 1
  endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; For Projection Region: 
  if ( ( region_type eq 'projection' ) ) then begin
        x1 = float( temp[ 1 ] )
        y1 = float( temp[ 2 ] )
        x2 = float( temp[ 3 ] )
        y2 = float( temp[ 4 ] )
        thick = float( temp[ 5 ] ) > 4
        angle = atan( abs( y2 - y1 ), abs( x2 - x1 ) )
        x3 = x2 - thick * sin( angle )
        y3 = y2 + thick * cos( angle )
        x4 = x1 - thick * sin( angle )
        y4 = y1 + thick * cos( angle )

        xx = [ x1, x2, x3, x4 ]
        yy = [ y1, y2, y3, y4 ]

        index_proj = polyfillv( xx, yy, naxis1, naxis2 )
        img_ori[ index_proj ] = -9999.9
        img_msk[ index_proj ] = 1
      endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; For Point Region: 
  if ( ( region_type eq 'point' ) ) then begin
        x1 = temp[ 1 ]
        y1 = temp[ 2 ]
        if ( number lt 4 ) then begin 
            radius = 11
        endif else begin 
            radius = float( temp[ 4 ] )
        endelse

        dist_circle, point_mask, [ naxis1, naxis2 ], x1, y1, /double
        index_point = where( point_mask le radius )
        img_ori[ index_point ] = -9999.9
        img_msk[ index_point ] = 1
  endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; For Annulus Region: 
  if ( ( region_type eq 'annulus' ) ) then begin
        x1 = temp[ 1 ]
        y1 = temp[ 2 ]
        r1 = temp[ 3 ]
        r2 = temp[ 4 ]
        dist_circle, annulus_mask, [ naxis1, naxis2 ], x1, y1, /double
        index_inner = where( annulus_mask le r1 )
        index_outer = where( annulus_mask le r2 )
        img_ori[ index_outer ] = -9999.9
        img_msk[ index_outer ] = 1
        img_ori[ index_inner ] = img_cor[ index_inner ]
        img_msk[ index_inner ] = 0
  endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; For Panda Region: 
  if ( ( region_type eq 'panda' ) ) then begin
        x1 = temp[ 1 ]
        y1 = temp[ 2 ]
        r1 = temp[ 6 ]
        r2 = temp[ 7 ]
        dist_circle, panda_mask, [ naxis1, naxis2 ], x1, y1, /double
        index_inner = where( panda_mask le r1 )
        index_outer = where( panda_mask le r2 )
        img_ori[ index_outer ] = -9999.9
        img_msk[ index_outer ] = 1
        img_ori[ index_inner ] = img_cor[ index_inner ]
        img_msk[ index_inner ] = 0
  endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print, '############################################'
print, 'Number of Box     Regions: ', n_box 
print, 'Number of Circle  Regions: ', n_circle 
print, 'Number of Ellipse Regions: ', n_ellipse
print, 'Number of Polygon Regions: ', n_polygon
print, 'Number of Lines   Regions: ', n_line
print, 'Number of Point   Regions: ', n_dot
print, 'Number of Annulus Regions: ', n_annulus
print, '############################################'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
img_msk[ where( img_msk ne 0 ) ] = 1
if keyword_set( savejpeg ) then begin 
    img_ori[ where( img_msk eq 1 ) ] = - 9999.9
    jpeg_file = 'newmask.jpeg'
    write_jpeg, jpeg_file, img_ori, quality = 100
endif
if keyword_set( overlap ) then begin 
    img_ori[ where( img_msk eq 1 ) ] = - 9999.9
    over_file = 'overlap.fits'
    if file_test( over_file ) then begin 
        spawn, 'rm ' + over_file 
    endif
    mwrfits, img_ori, over_file, head_ori
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; save the new mask file: 
msk_file = 'newmask.fits'
if NOT file_test( msk_file + '_backup' ) then begin 
    spawn, 'mv  ' + msk_file + '   ' + msk_file + '_backup'
endif
spawn, 'rm  ' + msk_file 
mwrfits, img_msk, msk_file, header_msk
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end
