function aniso_kappa_mv_density, species, deltaU, Phi, Bratio

  ;If A=1 same as Meyer-Vernet, Moncuquet and Huang,Icarus, 116, 202-213, 1995
  ;Following  Huang and Birmingham get the same A relation with magnetic latitude, 
  ;;Conserve Energy and First adiabtic Invariant and use Louiville's 
  ;;theorem to relate reference distribution function to at desired location
  ;;The perendicular temperature varies along the field line and is equal
  ; to T_{perp,0) / (A + (1-A) B0 / B)
  ;Assumed kappa was same for perp and parrallel and Tpar was constant along field line 

  
  n = (species.n/(species.A + (1 - species.A)/Bratio))*((1.d + (((-species.m*deltaU + species.q*Phi)/(species.kappa*species.T)) > (-1) ))^(0.5d - species.kappa))

  ;The > (-1) handles an unphysical case which would produce a NaN. With it, an
  ;iteration could temporarily wander into something unphysical but will still
  ;converge to the correct solution
  return,n
end
