;+
; :Author:
;  yinzp
;
; :Description:
;   Calculate the mie scattering matrix in different water layers.
;
; :Params:
;    nLayers: the number of the layers.
;    gamma: width of the gamma distribution.
;    Reff: the effective radius. Unit: micros
;    N0: the total number density of the water droplets. Unit: m^{-3}
;    nAngles: the number of the scattering angles. 
;    m: the relative radius of the droplet as to the wavelength of the incident photon.
;    lambda: the wavelength of the incident photon. Unit: nm
;
; :Keywords:
;    FILE: If set, the output results will be saved to the file.
;    
; :Reference:
;  [1] Natasha L Miles, Johannes Verlinde, and Eugene E Clothiaux, 
;      "Cloud droplet size distributions in low-level stratiform clouds," 
;      Journal of the atmospheric sciences 57 (2), 295-311 (2000).
;
; :Examples:
;
; :History:
;  2016-11-19
;-
;
;
;
PRO WaterMieScattering, nLayers, gamma, Reff, N0, nAngles, m, lambda, $
                        FILE = file
    
    S1Rel = Fltarr(nLayers, nAngles)
    S1Img = Fltarr(nLayers, nAngles)
    S2Rel = Fltarr(nLayers, nAngles)
    S2Img = Fltarr(nLayers, nAngles)
    muSca = Fltarr(nLayers)
    muExt = Fltarr(nLayers)
    
    x = 2.0*!PI*Reff*1E3/lambda   ; Relative radius.
    angles = Sequence(0.0, !PI, NELEMENTS=nAngles)*180.0/!PI
    
    FOR iLayers = 0, nLayers-1 DO BEGIN
        Print, 'Start to caculate the '+ String(iLayers, FORMAT='(I02)') + ' Layer!'
        FOR iAngles = 0, nAngles-1 DO BEGIN
            temp = MieScattering(m, x[iLayers], angles[iAngles])
            S1Rel[iLayers, iAngles] = Real_Part(temp[0])
            S1Img[iLayers, iAngles] = Imaginary(temp[0])
            S2Rel[iLayers, iAngles] = Real_Part(temp[1])
            S2Img[iLayers, iAngles] = Imaginary(temp[1])
            muExt[iLayers] = !PI*(Reff[iLayers]*1D-6)^2*N0[iLayers]*Real_Part(temp[2])
            muSca[iLayers] = !PI*(Reff[iLayers]*1D-6)^2*N0[iLayers]*Real_Part(temp[3])
        ENDFOR
    ENDFOR
    
    IF Keyword_Set(file) THEN BEGIN
        WriteH5, S1Rel, file, /FORCE, /OVERWRITE, VARNAME = 'S1Rel'
        WriteH5, S1Img, file, /APPEND, VARNAME = 'S1Img'
        WriteH5, S2Rel, file, /APPEND, VARNAME = 'S2Rel'
        WriteH5, S2Img, file, /APPEND, VARNAME = 'S2Img'
        WriteH5, muExt, file, /APPEND, VARNAME = 'Extinction'
        WriteH5, muSca, file, /APPEND, VARNAME = 'Scattering'
        WriteH5, angles*!PI/180.0, file, /APPEND, VARNAME = 'scaAngs'        
    ENDIF
END
                                 