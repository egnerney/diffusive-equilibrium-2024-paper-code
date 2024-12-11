function kappa_MV_version_noneq_density, species, deltaU, Phi, Bratio

  ;Meyer-Vernet, Moncuquet and Huang, Icarus, 116, 202-213, 1995
  
  
  ;a = sqrt(2.*species.T*(species.kappa - 1.5)/species.kappa/species.m)
  
  ;(2.*(deltaU-species.q*Phi/species.m)/species.kappa/a^2)
  n = species.n*((1.d + (((species.m*deltaU - species.q*Phi)/(species.kappa*species.T)) ))^(0.5d - species.kappa))
    
  ;The >0 handles an unphysical case which would produce a NaN. With it, an
  ;iteration could temporarily wander into something unphysical but will still
  ;converge to the correct solution
  return,n
end
