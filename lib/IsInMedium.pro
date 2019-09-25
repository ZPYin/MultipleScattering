;+
; :Params:
;  pos: the position of the photon. Unit: m
;  cloudBonud: the boundary of the cloud. Unit: m
;-
FUNCTION IsInMedium, pos, cloudBonud
	Return, ((pos[0]-cloudBonud[2])*(pos[0]+cloudBonud[2]) LT 0.0) AND $
		    ((pos[1]-cloudBonud[3])*(pos[1]+cloudBonud[3]) LT 0.0) AND $
		    ((pos[2]-cloudBonud[0])*(pos[2]-cloudBonud[1]) LT 0.0)
END