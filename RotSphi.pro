;--------------------------------------------------------------------------------------;
;                                   RotSphi
; SVIn: incident stokes vector
; phi: rotation angle. Unit: rad
;--------------------------------------------------------------------------------------;
FUNCTION RotSphi, SVIn, phi
    Return, [SVIn[0], SVIn[1]*COS(2.0*phi)+SVIn[2]*SIN(2.0*phi), $
             -SVIn[1]*SIN(2.0*phi)+SVIn[2]*COS(2.0*phi), SVIn[3]]
END
;--------------------------------------------------------------------------------------;