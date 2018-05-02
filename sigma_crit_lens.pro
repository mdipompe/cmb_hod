FUNCTION sigma_crit_lens,z,$
                         cmb_z=cmb_z,chi_cmb=chi_cmb
;+
;  NAME:
;    sigma_crit_lens
;
;  PURPOSE:
;    Get critical surface density for strong lensing (in *flat* cosmology!)
;
;  USE:
;    s_crit=sigma_crit_lens(z[,...])
;
;  INPUT:
;    z - redshift(s)
;
;  OPTIONAL INPUT:
;    cmb_z - redshift of CMB. Default: 1100
;    chi_cmb - co-moving distance to cmb_z. Calculated with
;              cosmocalc.pro if not provided.
;  
;  OPTIONAL KEYWORDS:

;  OUTPUT:
;    sigma_crit_lens
;
;  NOTES:
;    Uses MADs cosmology common block, initialized with
;    load_cosmology.pro.
;
;  HISTORY:
;    4/30/18 - Written - MAD (Dartmouth)
;-  
  
  COMMON cosmological_parameters

  ;MAD Put grav constant into (Mpc/h)^3/((M_sun/h) * s)
  G=(6.67408d-11)*((1.99d30/h))*(1./(3.0857d22/h))^3.

  ;MAD Put speed of light into Mpc/h/s
  c=(2.99792d8)*(1./(3.0857d22/h))

  ;MAD Get comoving distance at z (in Mpc/h)
  tmp=cosmocalc(z,d_c=chi)
  chi=chi*h
    
  ;MAD Get comoving distance to CMB (in Mpc/h)
  IF (n_elements(chi_cmb) EQ 0) THEN BEGIN
     IF (n_elements(cmb_z) EQ 0) THEN cmb_z=1100.
     tmp=cosmocalc(cmb_z,d_c=chistar)
     chistar=chistar*h
  ENDIF ELSE BEGIN
     chistar=chi_cmb
  ENDELSE

  ;MAD Get sigma_crit - note the factor of (1+z) is to get the 
  ;distance between source and lens, which is not just the difference 
  ;between the two distances.  Only works for flat cosmology!
  sigma_crit=((c^2.)/(4.*!dpi*G))*((chistar[0]*(1.+z))/(chi*(chistar[0]-chi)))

  
  return,sigma_crit

END
