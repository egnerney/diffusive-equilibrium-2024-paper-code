function aniso_maxwellian_density, species, deltaU, Phi, Bratio

;Based on Huang and Birmingham, JGR, 97, 1511-1519, 1992
;species.T is the parallel temperature
;species.A is T_{perp,0}/T_parallel
;The perendicular temperature varies along the field line and is equal
;	to T_{perp,0) / (A + (1-A) B0 / B)
;If the reference position, 0, does not coincide with the minimum magnetic
;	field along the field line, then A must be less than B0/(B0-Bmin).
;	Otherwise the perpendicular temperature will go to infinity at
;	some point.
;;Assumes Tpar was constant along field line or species.T

	n = species.n*exp((species.m*deltaU-species.q*Phi) $
                        /species.T)/(species.A + (1 - species.A)/Bratio)
	return,n
end
