;--------------------------------------------------------------------------------------;
;                                       Launch
; Mode: photon distribution mode. 'Planar'(default), 'Gauss'
; R: the radius of the laser beam in height h. Unit: m
; Seed: Seed pointer for randomU.
;--------------------------------------------------------------------------------------;
PRO Launch, Mode, R, seed, $
            X = x, Y = y
    
    CASE (Mode) OF
        'Planar': BEGIN
            x = 0 & y = 0
            WHILE (x^2+y^2 LE 1.0) DO BEGIN
                x = RandomU(*seed) & y = RandomU(*seed)
            ENDWHILE
            x = x*R & y = y*R
        END
        
        'Gauss': BEGIN
            theta = RandomU(*seed)*2.0*!PI
            r0 = RandomU(*seed, /NORMAL)*R
            x = r0*COS(theta) & y = r0*SIN(theta)
        END
        ELSE: BEGIN
            Message, 'Error in "Launch": ' + '"Mode" should be set in ' + $
                     '"Planar" and "Gauss"!'
        END
    ENDCASE

END
;--------------------------------------------------------------------------------------;