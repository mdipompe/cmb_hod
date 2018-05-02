FUNCTION sz_pressure,r,logm,z,cindex=cindex,rvir=rvir,c0=c0,mstar=mstar,m_200=m_200,r_200=r_200,p_200=p_200
;+
;  NAME:
;    sz_pressure
;
;  PURPOSE:
;    Gets the electron pressure profile of a DM halo for SZ effect
;    calculations.
;    See Battaglia et al. (2012)  
;
;  USE:
;    P=sz_pressure(r,logm,z[,...])
;
;  INPUT:
;    r - array of radii at which to calculate pressure 
;    logm - Log of the halo mass(es) (in M_sun/h). 
;    z - redshift(s)
;
;  OPTIONAL INPUT:
;    cindex - concentration index. If not supplied will calculate
;             using c_index.pro.
;    rvir - virial radius. If not supplied will calculate for you with r_vir.pro
;    c0 - concentration parameter, passed to c_index.pro if indices
;         not provided.
;    beta - power-law index for conentration indices, if not
;           provided.
;    mstar - mass normalization for concentration index, if not
;            provided.
;    m_200 - The virial mass converted to mass defined at
;            200*rho_crit. If not supplied, will be calculated with
;            convert_nfw.pro
;    r_200 - The virial radius converted to mass defined at
;            200*rho_crit. If not supplied, will be calculated with
;            convert_nfw.pro
;  
;  OPTIONAL KEYWORDS:
;    
;  OUTPUT:
;    1, 2, or 3D array of pressures. If arrays of both z and logM are
;    supplied, 3D array has values for each mass in first dimension,
;    and values for each z in second dimension, and r in last.
;
;    p_200 - if supplied, will pass back the pressure at r_200
;            (related to thermal and electron pressure)
;
;  NOTES:
;    Uses MADs cosmology common block, initialized with load_cosmology.pro
;
;
;  HISTORY:
;    4/30/18 - Written - MAD (Dartmouth)
;-    
  COMMON cosmological_parameters

  ;MAD Set some default NFW properties
  IF (n_elements(c0) EQ 0) THEN c0=11

  ;MAD Some matrix dimension orders are important for working with different
  ;size input combinations of M and z.
  IF (n_elements(logm) EQ 1 AND n_elements(z) GT 1) THEN dim1=n_elements(z) ELSE dim1=n_elements(logm)
  IF (n_elements(logm) EQ 1 AND n_elements(z) GT 1) THEN dim2=n_elements(logm) ELSE dim2=n_elements(z)

  ;MAD Copy r into matrix to do all M and z simultaneously
  IF (size(r,/n_dimensions) EQ 1) THEN BEGIN
     n_r=n_elements(r)
     r2=transpose(rebin([r],n_elements(r),n_elements(z),n_elements(logm)))
     IF (n_elements(logm) EQ 1 OR n_elements(z) EQ 1) THEN r2=reform([r2],n_elements(logm),n_elements(z),n_r)
  ENDIF ELSE BEGIN
     s=size(r)
     n_r=s[3]
     IF (s[1] NE n_elements(logm) OR s[2] NE n_elements(z)) THEN $
        message,'Dimensions of r not compatible with logm and/or z'
     r2=r
  ENDELSE

  ;MAD Get halo profile parameters if not supplied
  IF (n_elements(rvir) EQ 0) THEN $
     rvir=r_vir(logm,z)
  IF (n_elements(cindex) EQ 0) THEN $
     cindex=c_index(logm,z,c0=c0,mstar=mstar)
  
  ;MAD Get critical density at z
  rho_crit=rho_crit(z)
  
  rho_crit=reform(transpose(rebin([rho_crit],n_elements(z),n_elements(logm))))
  rho_crit=reform(rebin([rho_crit],dim1,dim2,n_r))
  
  ;MAD Get R_200, M_200
  IF (n_elements(m_200) EQ 0 OR n_elements(r_200) EQ 0) THEN $
     m_200=convert_nfw(logm,z,delta_cr=200.,r_deltacr=r_200)
  m_200=rebin([m_200],n_elements(logm),n_elements(z),n_r)
  r_s=r_200
  r_s=rebin([r_s],n_elements(logm),n_elements(z),n_r)
  x=r2/r_s

  ;MAD Set pressure profile parameters (from Battaglia et al. 2012
  ;equations 10 & 11 and table 1)
  zarr=transpose(rebin([z],n_elements(z),n_elements(logm)))
  IF (n_elements(logm) EQ 1 OR n_elements(z) EQ 1) THEN zarr=reform([zarr],n_elements(logm),n_elements(z))
  zarr=rebin([zarr],n_elements(logm),n_elements(z),n_r)
  
  gamma=-0.3
  P00=18.1
  alpha_mP=0.154
  alpha_zP=-0.758
  P0=P00*((m_200/10^14.)^alpha_mP)*((1.+zarr)^alpha_zP)

  x0=0.497
  alpha_mx=-0.00865
  alpha_zx=0.731
  xc=x0*((m_200/10^14.)^alpha_mx)*((1.+zarr)^alpha_zx)

  b0=4.35
  alpha_mb=0.0393
  alpha_zb=0.415
  beta=b0*((m_200/10^14.)^alpha_mb)*((1.+zarr)^alpha_zb)

  ;MAD Thermal pressure
  P_th_norm=P0*((x/xc)^gamma)*(1+(x/xc))^(beta*(-1.))

  ;MAD Get P200, multiply in to get pressure in physical units
  ;(see e.g. Battaglia et al. 2014 Eq 11)
  G=(6.67408d-11)*((1.99d30/h))*(1./(3.0857d22/h))^3. ;Grav constant in (Mpc/h)/((M_sun/h) * s^2)
  f_b=omega_b/omega_m
  p_200=(G*m_200*200.*rho_crit*f_b)/(2*r_s)
  P_th=P_th_norm*p_200
  
  ;Electron pressure
  P_e=P_th/1.932

  ;MAD Remove dimensions from M_200, P_200 for
  ;passing back
  m_200=m_200[*,*,0]
  p_200=p_200[*,*,0]
  
  return,P_e
END
