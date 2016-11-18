;+
; :Author:
;	yinzp
;
; :Description:
;	Calculate the scattered direction cosine of the scattered photon.
;
; :Params:
;	incDir: the direction cosine of the incident photon
;	phi: the angle between the meridian plane of the incident photon and the 
;		 scattering plane. Unit: rad
;	theta: the scatering angle. Unit: rad
;
; :Return:
;	the direction cosine of the scattered photon.
;
; :Date:
; 	2016-11-18
;-
;
;
;
FUNCTION UpdateDir, incDir, phi, theta
	IF (1.0-ABS(incDir[2]) LE 1E-12) THEN BEGIN
		Return, [SIN(theta)*COS(phi), $
				 SIN(theta)*SIN(phi), $
			 	 COS(theta)*(incDir[2] GE 0 ? 1:-1)]
	ENDIF ELSE BEGIN
		temp = Sqrt(1.0 - incDir[2]^2)
		Return, [ SIN(theta)*(incDir[0]*incDir[2]*COS(phi) - incDir[1]*SIN(phi))/temp $
			 	  + incDir[0]*COS(theta), $
		 	  	  SIN(theta)*(incDir[1]*incDir[2]*COS(phi) + incDir[0]*SIN(phi))/temp $
		 	  	  + incDir[1]*COS(theta), $
		 	  	  -SIN(theta)*COS(phi)*temp + incDir[2]*COS(theta) ]
	ENDELSE
END
