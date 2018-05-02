FUNCTION y_ell,logm,z,r=r,ell=ell,c0=c0,mstar=mstar,r_limit=r_limit,n_r=n_r
;+
;  NAME:
;    y_ell
;
;  PURPOSE:
;    Get the Fourier transform of an SZ pressure profile
;
;  USE:
;    y_l = y_ell(logm,z,[...])
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
;
;  OUTPUT:
;    1, 2 or 2D array of y_ell values. If multiple M and z passed,
;    3D array has M in first dimension, z in second, ell in third
;
;  NOTES:
;    Uses MADs cosmology common block, initialized with load_cosmology.pro
;
;  HISTORY:
;    4/30/18 - Written - MAD (Dartmouth)
;-  
  COMMON cosmological_parameters

  ;MAD Thompson cross section (Mpc/h)^2
  st=(6.652d-29)*(1./3.0857d22/h)^2
  ;MAD Speed of light (Mpc/h/s)
  c=2.99792458E8*(1./3.0857d22/h)
  ;MAD Electron mass in M_sun/h
  me=(9.109d-31)*(1./1.99d30/h)
  
  ;MAD Set some default NFW properties
  IF (n_elements(c0) EQ 0) THEN c0=11

  ;MAD Get virial radius for integration limits
  rvir=r_vir(logm,z)

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
  
  ;MAD Get  pressure profile for each mass
  p_e=sz_pressure(r[*,*,*,0],logm,z,c0=c0,mstar=mstar,r_200=r_s)
  p_e=rebin(p_e,n_elements(logm),n_elements(z),n_r,n_ell)
  
  ;MAD Put scale radius R_200 into array, set x
  r_s=rebin([r_s],n_elements(logm),n_elements(z),n_r,n_ell)

  x=r/r_s

  ;MAD Set angular scale
  dz=cosmocalc(z,d_a=da)
  da=transpose(rebin([da],n_elements(z),n_elements(logm)))
  IF (n_elements(z) EQ 1 OR n_elements(logm) EQ 1) THEN da=reform([da],n_elements(logm),n_elements(z))
  da=rebin([da],n_elements(logm),n_elements(z),n_r,n_ell)
  
  l_s=(da*h)/r_s
  
  step=diff_3d(x)
  integrand=(x^2.)*(sin((ell+0.5)*(x/l_s))/((ell+0.5)*(x/l_s)))*p_e
  int=total(integrand[*,*,1:(n_r-1),*]*step[*,*,*,*],3)
  int=reform(int,n_elements(logm),n_elements(z),n_ell)

;  x_test=x[0,0,*,0]
;  ell_test=ell[0,0,*,0]
;  l_s_test=l_s[0,0,*,0]
;  p_e_test=p_e[0,0,*,0]
;  test=int_tabulated(x_test,$
;                     (x_test^2.)*(sin((ell_test+0.5)*(x_test/l_s_test))/((ell_test+0.5)*(x_test/l_s_test)))*p_e_test,/double)
;  print,int[0,0,0],test,(int[0,0,0]-test)/test
;  stop
  
  r_s=reform(r_s[*,*,0,*],n_elements(logm),n_elements(z),n_ell)
  l_s=reform(l_s[*,*,0,*],n_elements(logm),n_elements(z),n_ell)
  
  y_l=(st/(me*c^2.))*((4.*!dpi*r_s)/(l_s^2.))*int

  ;MAD Put r, ell back into single array so it isn't passed back as huge grid
  ell=reform(ell[0,0,0,*])
  r=reform(r[0,0,*,0])
  
  return,y_l
END
