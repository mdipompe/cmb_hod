FUNCTION r_vir,logm,z,delta_cr=delta_cr
;+
;  NAME:
;    r_vir
;
;  PURPOSE:
;    Calculate the virial radius of DM halos
;
;  USE:
;    rvir = r_vir(logm,z,delta_cr=delta_cr)
;
;  INPUT:
;    logm - Log of the halo mass(es) (in M_sun/h). 
;    z - redshift(s)
;
;  OPTIONAL INPUT:
;    delta_cr - multiplicative factor relative to rho_crit that mass
;               is defined for. Default: virial Delta_cr for given z.
;  
;  OPTIONAL KEYWORDS:

;  OUTPUT:
;    1 or 2D array of r_vir. If arrays of both z and logM are
;    supplied, 2D array has values for each mass in first dimension,
;    and values for each z in second dimension.
;
;  NOTES:
;    Uses MADs cosmology common block, initialized with load_cosmology.pro
;
;  HISTORY:
;    4/30/18 - Written - MAD (Dartmouth)
;-  
  COMMON cosmological_parameters

  rho_crit=rho_crit(z)
  
  ;MAD If density contrast to define halo not supplied,
  ;get it (e.g. Hill & Spergel 2014) 
  IF (n_elements(delta_cr) EQ 0) THEN BEGIN
     ;MAD Evolve omega_m, get delta_cr (Hill & Spergel 2014)
     omega_z=evolve_density(z,type='matter')
     Delta_cr=(18.*!dpi^2) + (82*(omega_z-1.)) - (39.*(omega_z-1)^2.)
  ENDIF
     
  ;MAD Get R_vir, generalized for any dimensionality
  r_vir=reform(((3.*10.^logm)#(1./(4.*!dpi*Delta_cr*rho_crit)))^(1./3))

  IF (n_elements(r_vir) EQ 1) THEN r_vir=r_vir[0]
  return, r_vir
END
