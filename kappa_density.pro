function kappa_density, species, deltaU, Phi, Bratio

;Derived independently but consistent with Meyer-Vernet, Moncuquet and Huang,
;	Icarus, 116, 202-213, 1995
;Temperature is the pressure temperature, i.e. P = n k T, not the apparent
;	temperature of the core, f''(v=0) = f''_{Maxwellain}(v=0), which is
;	[Something I'll fill in later...]
;Temperature varies along the field line, with T/T0 = (n/n0)^(-1/(kappa-1/2))

	a = sqrt(2.*species.T*(species.kappa - 1.5)/species.kappa/species.m)
	n = species.n/$
		(1. + ((2.*(deltaU-species.q*Phi/species.m)/species.kappa/a^2)>0))$
		^(species.kappa - 0.5)
;The >0 handles an unphysical case which would produce a NaN. With it, an
;iteration could temporarily wander into something unphysical but will still
;converge to the correct solution
	return,n
end
