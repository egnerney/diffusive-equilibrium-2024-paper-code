function diff_eq_eddie,Pos_in,Pos0_in,Bratio,species0_in,phi,$
	posc,nchoke,phic,planet=planet,fcor=fcor,if_psh=if_psh,RJ=RJ

;Returns densities of all species at Pos, based on conditions at Pos0
;Input:
;Pos_in			Position to calculate densities
;			[radius, latitude, west longitude], radius is in
;			the same units as RJ, defaulting to 1 if RJ is not set.
;			Latitude and longitude are in degrees.
;Pos0_in		Position where reference
;			densities are specified, same units and coordinates as
;			Pos_in. It is up to the user to make sure Pos_in and
;			Pos0_in are along the same field line
;			Pos0_in does not need to be at the equator (geographic,
;			magnetic or centrifugal) but for if_psh=1, it must
;			be on the equatorward side of the choke points.
;			Also care must be taken with anisotropic Maxwellians,
;			since their perpendicular temperature is greatest at
;			where the magnetic field strenght is minimum. If
;			Pos0_in is at a different location and A is too large,
;			this might result in unphysical results (e.g. infinite
;			perpendicular temperature at the equator.)
;Bratio		[nd]	Ration of magnetic field strength at Pos_in to Pos0_in
;species0_in		Array of structures defining species at Pos0_in, 
;			electrons are the last structure in the array.
;			This should be defined as something like
;species0 = replicate({species,name:'O+',m:16.0,q:1.0,T:50.0,A=1.0,
;	kappa=10.,n:1800.,type:'Maxwellian'},nspec)
;			in the calling routine. The fields are:
;	name	[string]Name of species
;	m	[AMU]	Mass of species
;	q	[e]	Charge of species
;	T	[eV]	Parallel temperature of species
;	A	[nd]	Temperature anisotropy of species
;	kappa	[nd]    Kappa index, if the distribution is a kappa function
;	n	[any]	Density of species at equator, output will have the
;				same units
;	type	[string]Velocity distribution (E.g. 'Maxwellian')
;			A routine and file <type>_density.pro must be in
;			the working directory or IDL path
;posc		position of choke point (not used if if_psh=0)
;nchoke		number of choke points (not used if if_psh=0)
;			may be 0, 1 or 2
;planet		[string]Name of the planet (Jupiter or Saturn), Jupiter by 
;			default
;fcor		[nd]	Fraction of corotation rate (for subcorotation), 1.0
;			by default
;if_psh		[nd]	Flag, set to 1 to use a phase space hole above the
;			choke point, 0 by default
;			If if_psh is set to 1 one, A routine and file
;			<type>psh_density.pro must be in the working directory
;			or IDL path
;RJ                     The radius of Jupiter in the user's units. This
;                       defines the units used in Pos_in and Pos0_in
;
;Output:
;diff_eq()	[any]	Density of all ion species at Pos
;Phi		[V]	Ambipolar potential Pos
;Phichoke	[V]	Ambipolar potential at choke point (returns 0 if
;			torusward of choke point). This is actually an input
;			if if_psh=1 and pos is planetward of a choke point
;
;Version 1.0	From a long time ago... Calculations are for isotropic
;		Maxwellians from the torus.
;Version 1.1	Working on it... First add choke point with no filling of
;		velocity space hole FJC 11 MAR 20
;Version 1.1.1	Replace root-finder with bracketing/bisection FJC 23 APR 20
;Version 1.1.2  Fixed typo in centrifugal force equations. Set limits of
;		bracketing above choke point so arguments of errofc() are
;		positive (only consider phi which puts point above the
;		choke point for all species.) Added if_psh flag. FJC 5 OCT 20
;Version 1.2	Major change to inputs. Latitude and L have been replaced by
;		Pos_in and Pos0_in has been added to spefify the location of
;		the reference point (which still has to be equatorial.) This
;		will allow generalization to non-dipolar (or tilted dipole)
;		but it does break the choke point and phase space hole
;		calculations. So don't use those options with v 1.2. And
;		the unit of distance is defined by RJ, which needs to be 
;		changed to keep the code general (it ought to be RS if planet
;		is "Saturn"... FJC 28 FEB 21
;Version 1.2.1	Changed input values at reference point to an array of 
;		structure and removed electrons. Electrons are not the last
;		element in species0 FJC 27 MAR 21
;Version 1.2.2	General cleanup a mess and reverted some changes I made in
;		a rush before a  meeting and didn't describe well... 
;		FJC 11 MAY 22

;Known quirk: If pos_in = pos0_in, you will not get diff_eq()=species0_in and
;	phi=0 as you would expect. The bisection technique will start by
;	looking at charge neutrality at phi = +-dphi, and then iterate.
;	For some reason, this doesn't convege back to phi=0 exactly, as I
;	would expect. But it does converge to a value of phi which is charge
;	neutral to within the specified tolerance.

;Set optional arguments if not specified
if (n_elements(fcor) EQ 0) then fcor=1.0
if (n_elements(planet) EQ 0) then planet='Jupiter'
if (n_elements(if_psh) EQ 0) then if_psh=0
if (n_elements(RJ) EQ 0) then RJ = 7.1492e7

;Set planet-specific variabies
if ( define_planet(planet, RP, Omega, GM) ) then stop,'Planet '+$
	planet+' is not supported'
Omega *= fcor

;Turn species definitions at reference point into a local structure with
;	SI units
species0 = species0_in
species0.m *= 1.67d-27
species0.q *= 1.6d-19
species0.T *= 1.6d-19
nspec=n_elements(species0)

;Make sure the plasma is charge neutral at Pos0
;species0[nspec-1].n=$
 ;       -total(species0[0:nspec-2].q*species0[0:nspec-2].n)$
;        /species0[nspec-1].q
;print,species0[nspec-1].n
;print,species0[4].n
species0[nspec-1].n=$
        -total(species0[0:nspec-3].q*species0[0:nspec-3].n , /double)$
        /species0[nspec-1].q       ; removing hot electrons from this, something something hot protons to achieve their neutrality...  
       

;print,species0[nspec-1].n
;stop
;print,species0[4].n
;stop
;Get the gravitational-centrifugal potential at the location of interest and
;	the reference location
U0 = -GM/RP*(RJ/Pos0_in[0])-0.5d*RP^2d*(Pos0_in[0]/RJ)^2d*$
	Cos(Pos0_in[1]*!dtor)^2d*Omega^2d
U  = -GM/RP*(RJ/Pos_in[0]) -0.5d*RP^2d*(Pos_in[0]/RJ)^2d *$
	Cos(Pos_in[1]*!dtor )^2d*Omega^2d

;This is ugly but it should work

torus_side = 1 ; torus side no choke points no phase space holes

if_psh = 0   ; no phase space hole (psh)
;torus_side:	0	Pos is planetward of the southern choke point
;		1	Pos is between the choke points, or phase space holes
;				are not being used
;		2	Pos is planetward of the north choke point

;Get potentials at the choke point

;Initialize some variables needed for the bracketing/bisection root finder
n = dblarr(nspec)
n2= dblarr(nspec)
nm= dblarr(nspec)

Phi = 0d0
Phi2= 0d0
nq = 1d0
nq2= 1d0
dPhi = 1d0 

;if (torus_side) then begin
	;Move Phi and Phi2 apart until the charge densities has opposite signs
	while (nq*nq2 GT 0d) do begin
		Phi -= dPhi
		Phi2 += dPhi
		for i=0,nspec-1 do n[i] = call_function( $
		   species0[i].type+'_density',species0[i],U0-U,Phi,Bratio)
		   
		  nq_temp = [species0[0:nspec-3].q*n[0:nspec-3],species0[nspec-1].q*n[nspec-1]] ; removing hot electrons from this
		  
		nq = total(nq_temp , /double);total(species0.q*n) 
		for i=0,nspec-1 do n2[i] = call_function( $
                   species0[i].type+'_density',species0[i],U0-U,Phi2,Bratio)
                   
                nq2_temp = [species0[0:nspec-3].q*n2[0:nspec-3],species0[nspec-1].q*n2[nspec-1]]   ; removing hot electrons from this
                nq2 = total(nq2_temp , /double);total(species0.q*n2)
	endwhile

	while (abs(nq/n[nspec-1]/species0[nspec-1].q) GT 1d-5) do begin ; change back to -5 for general use
		;Bisect until charge density is under 0.1% of electron density
		;Evaluate midpoint charge density
		Phim = 0.5*(Phi+Phi2)
		for i=0,nspec-1 do nm[i] = call_function( $
                   species0[i].type+'_density',species0[i],U0-U,Phim,Bratio)
                   
                   nqm_temp = [species0[0:nspec-3].q*nm[0:nspec-3],species0[nspec-1].q*nm[nspec-1]]; removing hot electrons from this
                nqm = total(nqm_temp , /double);total(species0.q*nm)

		;move Phi or Phi2 to Phim, so that the root is between the
		;new Phi and Phi2
		if (nq*nqm GT 0d) then begin
			Phi=Phim
			n=nm
			nq = nqm
		endif else begin
			Phi2=Phim
                        n2=nm
                        nq2 = nqm
		endelse
	endwhile



return,n

end
