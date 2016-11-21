;+
; :Params:
;  pos: the position of the photon. Unit: m
;  FOV: the field of view. Unit: mrad
;  diaTel: the diameter of the telescope. Unit: mm
;  Dx: the distance between the centre of the telescope and the laser beam. Unit: mm
;-
FUNCTION IsInFOV, pos, FOV, diaTel, Dx
	diaH = FOV * pos[2]/1000.0 + diaTel/1000.0  ; the diameter of the FOV at H. 
													; Unit: m
	IF ((pos[0]-Dx/1000.0)^2.0+pos[1]^2.0) LT (diaH/2.0)^2 THEN Return, 1
	Return, 0
END