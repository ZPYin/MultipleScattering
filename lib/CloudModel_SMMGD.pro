;+
; :Author:
;  yinzp
;
; :Description:
;   Compute the single-mode modified gamma distribution.
;
; :Params:
;    LWC: liquid water content. Unit: kg*m^{-3}
;    Reff: effective radius. Unit: m
;    No: total number of particles. Unit: m^{-3}
;    R: the input radius: Unit: m
;
; :Reference:
;   [1] Natasha L Miles, Johannes Verlinde, and Eugene E Clothiaux, 
;   "Cloud droplet size distributions in low-level stratiform clouds," 
;   Journal of the atmospheric sciences 57 (2), 295-311 (2000).
;   [2] DP Donovan, Henk Klein Baltink, JS Henzing, SR De Roode, 
;   and AP Siebesma, "A depolarisation lidar-based method for the 
;   determination of liquid-cloud microphysical properties," 
;   Atmospheric Measurement Techniques, 8 (1), 2015 (2015).
;
; :Examples:
;    Rarr = Sequence(0.0, 8.0, 0.0001)*1E-6
;    Reff = 3.0E-6
;    LWC = 1.0E-4
;    No = 2000E6
;    P1 = Plot(Rarr, CloudModel_SMMGD(LWC, Reff, No, Rarr), XRANGE=[Rarr[0], Rarr[-1]])
;    
; :History:
;  2016-11-14
;-
;
;
;
FUNCTION CloudModel_SMMGD, LWC, Reff, No, R

;--------------------------------------------------------------------------------------;
;                                     Input check
;--------------------------------------------------------------------------------------;
    ;ON_ERROR, 2

    ; syntax checking    
    syntax = 'Res = CloudModel_SMMGD(gamma, Reff, No, R)'
    
    IF N_PARAMS() NE 4 THEN $
        MESSAGE, 'Error in ''CloudModel_SMMGD'': '+ syntax 
    IF (LWC LT 0) OR (Reff LT 0) OR (No LT 0) THEN $
        Message, 'Error in ''CloudModel_SMMGD'': '+ $
                 'Input parameters must be positive!'
;--------------------------------------------------------------------------------------;

    Rv = (LWC*3.0/(4.0*!PI*No*1E3))^(1.0/3.0)
    k = Rv^3/Reff^3
    gamma = ((1.0-4.0*k)-Sqrt(8.0*k+1.0)) / (2.0*(k-1))   ; width
    IF (gamma LT 0) OR (~Finite(gamma)) THEN $
        Message, 'Error in ''CloudModel_SMMGD'': ' + 'Gamma is negative or NaN!'
    Rm = Reff/(gamma+2.0)
    Return, No/Rm * 1.0/Factorial(gamma-1) * (R/Rm)^(gamma-1) * Exp(-R/Rm)
END