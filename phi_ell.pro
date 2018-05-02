FUNCTION phi_ell,logm,z,r=r,ell=ell,c0=c0,mstar=mstar,type=type,r_limit=r_limit,n_r=n_r
;+
;  NAME:
;    phi_ell
;
;  PURPOSE:
;    Get the Fourier transform of an NFW profile
;
;  USE:
;    phi_l = phi_ell(logm,z,[...])
;
;  INPUT:
;    logm - Log of the halo mass(es) (in M_sun/h). 
;    z - redshift(s)
;
;  OPTIONAL INPUT:
;    r - radii to calculate NFW profile. Defaults to evenly
;        log-spaced, with a max at 1.5*r_vir and a minimum 10^5
;        smaller than this.
;    ell - ell values. Defaults to evenly log-spaced over 1-10,000
;    c0 - concentration parameter, passed to c_index.pro
;    mstar - mass normalization for concentration index
;    r_limit - maximum radii to get NFW profile. Default: 1.5*r_vir
;    n_r - number of radii in calculation. Default: 2000
;  
;  OPTIONAL KEYWORDS:
;    type - Dark matter ('DM'; default) or number ('num') density
;
;  OUTPUT:
;    1, 2 or 2D array of phi_ell values. If multiple M and z passed,
;    3D array has M in first dimension, z in second, ell in third
;
;  NOTES:
;    Uses MADs cosmology common block, initialized with load_cosmology.pro
;
;  HISTORY:
;    4/30/18 - Written - MAD (Dartmouth)
;-  
  COMMON cosmological_parameters

  ;MAD Set default type
  IF (n_elements(type) EQ 0) THEN type='dm'
  
  ;MAD Set some default NFW properties
  IF (n_elements(c0) EQ 0) THEN BEGIN
     IF (type EQ 'num') THEN c0=32
     IF (type EQ 'dm') THEN c0=11
  ENDIF

  ;MAD Get virial radius, concentration index, projected NFW profile
  rvir=r_vir(logm,z)
  c=c_index(logm,z,c0=c0,mstar=mstar)
  
  ;MAD Set default upper integration limit at 1.5*rvir
  IF (n_elements(r_limit) EQ 0) THEN r_limit=1.5*rvir
  IF (n_elements(logm) EQ 1) THEN r_limit=1#r_limit
  ;MAD Set default number of points in r_grid
  IF (n_elements(n_r) EQ 0) THEN n_r=2000
  
  ;MAD Set default radii and ell values
  IF (n_elements(ell) EQ 0) THEN ell=cgloggen(1000,start=1,finish=10000)
  n_ell=n_elements(ell)
  IF (n_elements(r) EQ 0) THEN BEGIN
     r=fltarr(n_elements(logm),n_elements(z),n_r)
     FOR i=0L,n_elements(logm)-1 DO BEGIN
        FOR j=0L,n_elements(z)-1 DO BEGIN
           r[i,j,*]=cgloggen(n_r,start=r_limit[i,j]/10000.,finish=r_limit[i,j])
        ENDFOR
     ENDFOR
  ENDIF
  r=rebin([r],n_elements(logm),n_elements(z),n_r,n_ell)
  ell=reform(transpose(rebin([ell],[n_ell,n_r,n_elements(z),n_elements(logm)])),$
             n_elements(logm),n_elements(z),n_r,n_ell)
  
  ;MAD Get distance to source (Mpc/h)
  dz=cosmocalc(z,d_a=da)
  da=transpose(rebin([da],n_elements(z),n_elements(logm)))
  IF (n_elements(z) EQ 1 OR n_elements(logm) EQ 1) THEN da=reform([da],n_elements(logm),n_elements(z))
  da=rebin([da],n_elements(logm),n_elements(z),n_r,n_ell)
  
  ;MAD Get NFW profile for each mass and redshift
  rho=nfw(r[*,*,*,0],logm,z,rvir=rvir,cindex=c,mstar=mstar,type=type)
  IF (n_elements(logm) EQ 1 OR n_elements(z) EQ 1) THEN rho=reform(rho,n_elements(logm),n_elements(z),n_r)
  rho=rebin(rho,n_elements(logm),n_elements(z),n_r,n_ell)

;;  MAD Some tests - do you get back the input masses? (change r_limit
;;                   above to r_vir)
;  step=diff_3d(r[*,*,*,0])
;  step=rebin(step,n_elements(logm),n_elements(z),n_r-1,n_ell)
;  mass=total(4.*!dpi*r[*,*,1:n_r-1,*]^2*rho[*,*,1:n_r-1,*]*step[*,*,*,*],3)
;  stop
  
  ;MAD set scale radius, x, angular scale
  r_s=rvir/c
  r_s=rebin(reform([r_s],n_elements(logm),n_elements(z)),$
            n_elements(logm),n_elements(z),n_r,n_ell)

  x=r/r_s
  l_s=(da*h)/r_s
  
  ;MAD Get crit surface density if not supplied
  IF (n_elements(sigma_crit) EQ 0) THEN BEGIN
     IF (type EQ 'dm') THEN BEGIN
        ;MAD Get distance to CMB, crit surface density
        sigma_crit=sigma_crit_lens(z)
        sigma_crit=rebin(transpose(rebin([sigma_crit],n_elements(z),n_elements(logm))),$
                         n_elements(logm),n_elements(z),n_r,n_ell)
     ENDIF
     IF (type EQ 'num') THEN BEGIN
        sigma_crit=fltarr(n_elements(logm),n_elements(z),n_r,n_ell)+100.
     ENDIF
  ENDIF

  step=diff_3d(x)
  integrand=(x^2.)*(sin((ell+0.5)*(x/l_s))/((ell+0.5)*(x/l_s)))*(rho/sigma_crit)
  int=total(integrand[*,*,1:(n_r-1),*]*step[*,*,*,*],3)
  int=reform(int,n_elements(logm),n_elements(z),n_ell)

;  tell=500
;  tz=2
;  tm=3
;  x_test=x[tm,tz,*,tell]
;  ell_test=ell[tm,tz,*,tell]
;  l_s_test=l_s[tm,tz,*,tell]
;  rho_test=rho[tm,tz,*,tell]
;  sigma_crit_test=sigma_crit[tm,tz,*,tell]
;  test=int_tabulated(x_test,$
;                     (x_test^2.)*(sin((ell_test+0.5)*(x_test/l_s_test))/((ell_test+0.5)*(x_test/l_s_test)))*(rho_test/sigma_crit_test),/double)
;  print,int[tm,tz,tell],test,(int[tm,tz,tell]-test)/test
;  stop
  
  ell=reform(ell[*,*,0,*],n_elements(logm),n_elements(z),n_ell)
  r_s=reform(r_s[*,*,0,*],n_elements(logm),n_elements(z),n_ell)
  l_s=reform(l_s[*,*,0,*],n_elements(logm),n_elements(z),n_ell)
  
  phi_l=(2./(ell*(ell+1)))*((4.*!dpi*r_s)/(l_s^2.))*int
;  phi_l=((4.*!dpi*r_s)/(l_s^2.))*int

  ;MAD Put r, ell back into single array so it isn't passed back as huge grid
  ell=reform(ell[0,0,*])
  r=reform(r[0,0,*,0])

  return,phi_l
END
