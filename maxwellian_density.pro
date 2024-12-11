function maxwellian_density, species, deltaU, Phi, Bratio
;;Assumes Tpar was constant along field line or species.T
	n = species.n*exp((species.m*deltaU-species.q*Phi) $
                        /species.T)
	return,n
end
