;+
; :Author:
;  yinzp
;
; :Description:
;   compute the log-gamma distribtuion.
;
; :Reference:
;   [1] EW Eloranta, presented at the Twenty-First International Laser Radar 
;   Conference Proceedings, 2002 (unpublished).
;   
; :Params:
;    Mode: choose the 'log-normal' distribution or 'gamma' distribution
;
; :Keywords:
;    BETASCA: scattering cross section. Unit: m^{-1} 
;    GAMMA: power in exponential of gamma distribution. or the width of the lognm 
;           distribution. Unit: Null
;    REFF: effective radius. Unit: microns
;    RADIUS: input radius. Unit: microns
;    ALPHA: power law parameter in the gamma distribution. Unit: Null
;
; :Examples:
;    Rarr = Sequence(0.1, 8.0, 0.01)*1E-6    
;    P1 = Plot(Rarr, CloudModel_LGD('lognm', $
;        BETASCA=1e-6, GAMMA=1, REFF=2e-6, Radius=rarr), XRANGE=[Rarr[0], Rarr[-1]])
;
; :History:
;  2016-11-16
;-
;
;
;
FUNCTION CloudModel_LGD, Mode, $
						 BETASCA=betaSca, GAMMA=gamma, REFF=Reff, RADIUS=Radius, $
						 ALPHA = alpha

	; Input check
	syntax = 'Res = CloudModel_LGD(Mode)'
	IF N_Params() NE 1 THEN $
		Message, 'Error in "CloudModel_LGD": ' + syntax
    IF (~Keyword_Set(betaSca)) OR (~Keyword_Set(gamma)) OR (~Keyword_Set(Reff)) OR $
       (~Keyword_Set(Radius)) THEN $
        Message, 'Error in "CloudModel_LGD": ' + 'Keyword "BETASCA", "GAMMA", ' + $
                 '"REFF" and "RADIUS" must be set!'
	IF (betaSca LT 0) OR (gamma LT 0) OR (Reff LT 0) THEN $
		Message, 'Error in "CloudModel_LGD": ' + 'Input parameters ' + $
		'must be non-negative!'

	CASE Mode OF
		'lognm': BEGIN
			alpha = Reff*Exp(-2.5*gamma^2.0)
			a = betaSca/(2.0*!PI*Reff^2.0)*Exp(-2.0*gamma^2.0)

			Return, a/(Sqrt(2.0*!PI)*gamma*Radius)* $
			        Exp(-0.5*(Alog(Radius / alpha)/gamma)^2.0)
		END
		'gamma': BEGIN
			IF ~Keyword_Set(alpha) THEN $
				Message, 'Error in "CloudModel_LGD": ' + 'Keyword_Set "ALPHA" must '+ $
						 'be set when "Mode" eq "gamma"!'
			p3 = (alpha+3.0)/gamma
			p4 = (alpha+4.0)/gamma
			b = ((1.0/Reff)*Exp(Lngamma(p4))/Exp(Lngamma(p3)))^gamma
			a = betaSca*gamma*b^p3/(2.0*!PI*Exp(Lngamma(p3)))

			Return, a*Radius^alpha*Exp(-Radius*R^gamma)
		END
	ENDCASE 
END