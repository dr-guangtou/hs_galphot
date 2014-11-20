FUNCTION hs_read_elliptab, PROF, NUM_SKIP=NUM_SKIP, $
    PROF_LINE=PROF_LINE

on_error, 2
compile_opt idl2

if N_params() lt 1 then begin 
    print,  'Syntax - sbp = read_profile( prof, SKIPLINE=SKIPLINE, prof_line=PROF_LINE'
    message, 'Check the Syntax Again !!'
    return, -1
endif

if keyword_set( NUM_SKIP ) then begin 
    NUM_SKIP = fix( NUM_SKIP ) 
endif else begin 
    NUM_SKIP = 41L 
ENDELSE
STR_SKIP = STRCOMPRESS( string( NUM_SKIP ), /remove_all )

PROFILE_FILE = STRCOMPRESS( PROF, /remove_all )

       DATA = { sma:0.0, intens:0.0, intens_err:0.0, pix_var:0.0, rms:0.0, ell:0.0, ell_err:0.0, pa:0.0, pa_err:0.0, x0:0.0,x0_err:0.0, y0:0.0, y0_err:0.0, grad:0.0, grad_err:0.0, grad_r_err:0.0, rsma:0.0, mag:0.0, mag_lerr:0.0, mag_uerr:0.0, tflux_e:0.0, tflux_c:0.0, tmag_e:0.0, tmag_c:0.0, npix_e:0L, npix_c:0L, a3:0.0, a3_err:0.0, b3:0.0, b3_err:0.0, a4:0.0, a4_err:0.0, b4:0.0, b4_err:0.0, ndata:0L, nflag:0L, niter:0L, stopp:0L, a_big:0.0, sarea:0.0 }
       
        SBP = DATA

IF NOT FILE_TEST( PROFILE_FILE ) THEN BEGIN
    
    print, ' The Profile file can not be found ! '
     
    return, -1

ENDIF ELSE BEGIN

      PROF_FILE = STRCOMPRESS( PROFILE_FILE )

      ;; By Song 12/08/2013 
      spawn, 'cp ' + PROF_FILE + ' ' + PROF_FILE + '_back'
      spawn, 'sed -i "s/INDEF/-99999.9/g"  ' + PROF_FILE + '_back'
      spawn, 'sed -i "1,' + STR_SKIP + 'd" ' + PROF_FILE + '_back'
 
      PROF_LINE = FILE_LINES( PROF_FILE + '_back' )

      PRINT, 'There are ', PROF_LINE, ' lines in the profile file'

  IF ( PROF_LINE le 0 ) THEN BEGIN
   
        PRINT, 'There is something wrong with the profile file' 
        
        return, -1        

  ENDIF ELSE BEGIN

        SBP = REPLICATE( DATA, PROF_LINE )

        OPENR, LUN, PROF_FILE + '_back', /GET_LUN
        READF, LUN, SBP
        CLOSE, LUN
        FREE_LUN, LUN

        return, sbp

  ENDELSE
 
ENDELSE

END
