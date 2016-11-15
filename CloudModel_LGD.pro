FUNCTION CloudModel_LGD, Mode, betaSca, gamma, Reff, R

	; Input check
	syntax = 'Res = CloudModel_LGD(Mode, betaSca, gamma, Reff, R)'
	IF N_Params() LT 5 THEN $
		Message('Error in "CloudModel_LGD": ') + syntax
	IF Isa(Mode, /STRING) THEN $
		Message('Error in "CloudModel_LGD": ') + $
		Varname_Scope(Mode, LEVEL=1) + ' must be a string!'
	IF (betaSca LT 0) OR (gamma LT 0) OR (Reff LT 0) OR (R LT 0) THEN
		Message('Error in "CloudModel_LGD": ') + 'Input parameters ' + $
		'must be non-negative!'

	alpha = Reff*Exp(-2.5*gamma^2.0)
	a = betaSca/(2.0*!PI*Reff^2.0)*Exp(-2.0*gamma^2.0)

	Return, a/(Sqrt(2.0*!PI)*gamma*R)*Exp(-0.5*(Alog(R / alpha)/R)^2.0)
END