FUNCTION c_index,logm,z,c0=c0,beta=beta,mstar=mstar,get_mstar=get_mstar,duffy=duffy
;+
;  NAME:
;    c_index
;
;  PURPOSE:
;    Get concentration index for an NFW profile at M, z, using
;    parameterization of Zheng et al. (2007; default) or Duffy et al. (2008)
;
;  USE:
;    c=c_index(logm,z[,c0=c0,beta=beta,mstar=mstar,get_mstar=get_mstar,duffy=duffy])
;
;  INPUT:
;    logm - Log of the halo mass(es) (in M_sun/h). 
;    z - redshift
;
;  OPTIONAL INPUT:
;    c0 - parameter controlling how concentrated distribution
;         is. Defalt: 11 (Zheng et al. 2007)
;    beta - Mass power-law index. Default: -0.13
;    mstar - Mass normalization. Default: 3.71535d+12
;  
;  OPTIONAL KEYWORDS:
;    get_mstar - If set, will get M* for your cosmology (note that
;                this is slow, so you may want to get it once and then
;                just set M*)
;    duffy - If set, uses the equations of Duffy et al. (2008) instead
;            of Zheng et al. (2007).
;
;  OUTPUT:
;    1 or 2D array of concentrations. If arrays of both z and logM are
;    supplied, 2D array has values for each mass in first dimension,
;    and values for each z in second dimension.
;
;  NOTES:
;    Uses MADs cosmology common block, initialized with load_cosmology.pro
;
;
;  HISTORY:
;    4/30/18 - Written - MAD (Dartmouth)
;-  

  COMMON cosmological_parameters

  
  IF (n_elements(duffy) EQ 0) THEN BEGIN
     ;Power-index for Mass dependence
     IF (n_elements(beta) EQ 0) THEN beta=-0.13
     ;Normalization of concentration index
     IF (n_elements(c0) EQ 0) THEN c0=11
     ;Normalization M*, set for default cosmology
     IF (n_elements(mstar) EQ 0 and n_elements(get_mstar) EQ 0) THEN mstar=3.71535d+12
  
     ;MAD Set value of Mstar, by equating
     ;sigma(M,z=0)=delta_crit(z=0)...I *think* this is right 
     ;based on Zheng et al. 2007.
     IF (n_elements(get_mstar) NE 0) THEN BEGIN
        ;MAD Set the dimensionless critical density for collapse at z=0
        delta_crit=(0.15*((12*!dpi)^(2./3))*(omega_m^(0.0055)))
        ;MAD Make some test masses, get sigma(M,z=0)
        masses=findgen(1000)/100+10
        sm=sigma_m(masses,0)
        ;MAD Find closest match to delta_crit(z)
        Mstar=10.^masses[nearest(delta_crit,sm)]
        Mstar=Mstar[0]
     ENDIF
  
     ;MAD Get concentration index, generalized for any dimensionality
     c=reform((((10.^logm)/Mstar)^beta)#(c0/(1.+z)))
  ENDIF ELSE BEGIN
     M_pivot=2d12
   
     A=7.85
     B=-0.081
     bigC=-0.71
     
     c=A*((((10.^logm)/M_pivot)^B)#((1.+z)^bigC))
     xx=where(z EQ 0, cnt)
     IF (cnt NE 0) THEN BEGIN
        A=7.96
        B=-0.091
        bigC=0
        cz0=A*((((10.^logm)/M_pivot)^B)*((1.+0.)^bigC))
        c[*,xx]=cz0
     ENDIF
     c=reform(c)
  ENDELSE

  IF (n_elements(c) EQ 1) THEN c=c[0]
  return,c
END
