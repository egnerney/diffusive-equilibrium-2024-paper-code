function fried_egg_density, species, deltaU, Phi, Bratio

  ;Based Fried Egg Dist Bagenal, Crary, Hill 2001 Unpublished. RIP Crary. I miss you friend.
;/(species.A + (1 - species.A)/Bratio))*
if (bratio lt 1) then bratio = 1d 
alpha = species.kappa/(species.kappa + 1d)
lambda_ = (Bratio - 1d)*(species.kappa + 1d)
;if (lambda_ lt 0) then print,'whooppsy'
;if (lambda_ lt 0) then stop
gfunc = GAMMA(-species.kappa)*(1 - IGAMMA(-species.kappa, lambda_));IGamma(-species.kappa,  lambda_)
g = species.kappa*(1d + lambda_/(species.kappa + 1d))*(lambda_^species.kappa)*Exp(lambda_)*gfunc
  n = species.n*g*Exp(((species.m*deltaU - species.q*Phi)/(alpha*species.kappa*species.T)) < 0 )


  ;The >(0) handles an unphysical case which would produce a NaN. With it, an
  ;iteration could temporarily wander into something unphysical but will still
  ;converge to the correct solution
  return,n
end
