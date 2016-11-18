PRO ScaRotAng, incDir, scaDir, $
			   SCAANG = scaAng, ROTANG1 = rotAng1, ROTANG2 = rotAng2

	syntax = 'ScaRotAng, incDir, scaDir, ' + $
		     'SCAANG = scaAng, ROTANG1 = rotAng1, ROTANG2 = rotAng2'
	IF N_Params() NE 2 THEN $
		Message, 'Error in "ScaRotAng": ' + syntax

	phiInc = Atanxoy(incDir)   ; Azimuth. Unit: rad
	phiSca = Atanxoy(scaDir)   ; Azimuth. Unit: rad
	thetaInc = ACOS(incDir[2]) ; Zenith. Unit: rad
	thetaSca = ACOS(scaDir[2]) ; Zenith. Unit: rad

	scaAng = ACOS(Total(incDir, scaDir))   ; Scattering angle. Unit: rad

	IF ABS((SIN(scaAng))) LT 1E-12 THEN BEGIN
		rotAng1 = 0.0 & rotAng2 = 0.0
		Return
	ENDIF

	IF (1.0 - ABS(incDir[2])) LE 1E-12 THEN BEGIN
		rotAng1 = ACOS(incDir[2]*COS(phiSca - phiInc))
		rotAng2 = incDir[2]
		Return
	ENDIF

	IF (1.0 - ABS(scaDir[2])) LE 1E-12 THEN BEGIN
		rotAng1 = ACOS(scaDir[2]*COS(phiSca - phiInc))
		rotAng2 = scaDir[2]
		Return
	ENDIF

	rotAng1 = ACOS((scaDir[2] - incDir[2]*COS(scaAng))/(SIN(thetaInc)*SIN(scaAng)))
	rotAng2 = ACOS((incDir[2] - scaDir[2]*COS(scaAng))/(SIN(thetaSca)*SIN(scaAng)))

END