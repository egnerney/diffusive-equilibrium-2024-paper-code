function aniso_kappa_density, species, deltaU, Phi, Bratio

  ;If A=1 same as Meyer-Vernet, Moncuquet and Huang,Icarus, 116, 202-213, 1995
  ;Following  Huang and Birmingham get the same A relation with magnetic latitude, 
  ;;Conserve Energy and First adiabtic Invariant and use Louiville's 
  ;;theorem to relate reference distribution function to at desired location
  ;;The perendicular temperature varies along the field line and is equal
  ; to T_{perp,0) / (A + (1-A) B0 / B)
  ;Assumed kappa was same for perp and parrallel and Tpar was constant along field line 

  a = sqrt(2.*species.T*(species.kappa - 1.5)/species.kappa/species.m)
  n = species.n/$
    (1. + ((2.*(deltaU-species.q*Phi/species.m)/species.kappa/a^2)>0))$
    ^(species.kappa - 0.5)
  ;The >0 handles an unphysical case which would produce a NaN. With it, an
  ;iteration could temporarily wander into something unphysical but will still
  ;converge to the correct solution
  return,n
end
