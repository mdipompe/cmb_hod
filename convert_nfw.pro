FUNCTION convert_nfw,logm,z,delta_cr=delta_cr,r_deltacr=r_deltacr,c_deltacr=c_deltacr,wrt_m=wrt_m
;+
;  NAME:
;    convert_nfw
;
;  PURPOSE:
;    Convert parameters (M, r, c) of an NFW profile from virial values
;    to ones relative to a different overdensity limit. Defaults to be
;    relative to critical density, unless wrt_m is set, then relative
;    to mean density.
;
;  USE:
;    new_nfw=convert_nfw(logm,z[,delta_cr=delta_cr,r_deltacr=r_deltacr,c_deltacr=c_deltacr,wrt_m=wrt_m])
;
;  INPUT:
;    logm - Log of the halo mass(es) (in M_sun/h). 
;    z - redshift
;
;  OPTIONAL INPUT:
;    delta_cr - the new relative overdensity. Default: 200.
;  
;  OPTIONAL KEYWORDS:
;    wrt_m - set to do calculations with delta_cr relative to mean,
;            rather than critical, density
;
;  OUTPUT:
;    Function outputs new masses that correspond to input virial
;    masses for the new definition.
;    r_deltacr - halo radii with new definiton
;    c_deltacr - concentration indices with new definition
;
;  NOTES:
;    Uses MADs cosmology common block, initialized with load_cosmology.pro
;
;  HISTORY:
;    4/30/18 - Written - MAD (Dartmouth)
;-  
  COMMON cosmological_parameters

  ;MAD Set default factor converting to
  IF (n_elements(delta_cr) EQ 0) THEN delta_cr=200.

  ;MAD Get concentration index, virial radius
  c=c_index(logm,z)
  rvir=r_vir(logm,z)

  ;MAD Get virial Delta
  omega_z=evolve_density(z,type='matter')
  Delta_cr_v=(18.*!dpi^2) + (82*(omega_z-1.)) - (39.*(omega_z-1)^2.)
  

  ;MAD If converting just to new delta but still relative to the 
  ;MAD critical density, use Appendix eqs of Hu and Kravtsov 99
  ;MAD (faster and more accurate)
  ;MAD https://arxiv.org/pdf/astro-ph/0203169.pdf
  IF (n_elements(wrt_m) EQ 0) THEN BEGIN

     delta_cr_v=transpose(rebin([delta_cr_v],n_elements(z),n_elements(logm)))
     
     f=(delta_cr/delta_cr_v)*(1./c^3.)*(alog(1+c)-(1./(1+(1./c))))
  
     a1=0.5116
     a2=-0.4283
     a3=-3.13d-3
     a4=-3.52d-5
     p=a2+a3*alog(f)+a4*(alog(f)^2)
     xf=((a1*f^(2.*p)+(3./4.)^2.)^(-0.5))+(2*f)

     rd_rv=(1./(c*xf))
     r_deltacr=reform(rd_rv*rvir)

     tmp=rebin([logm],n_elements(logm),n_elements(z))
     m_deltacr=reform((10.^tmp)*(delta_cr/delta_cr_v)*(rd_rv)^3.)
     
     c_deltacr=reform(c*rd_rv)
  ENDIF ELSE BEGIN
     ;MAD if converting to relative to mean density, use eq 4 of
     ;MAD Hill & Pajer 13, with the integral solved analytically.
     ;MAD (still have to solve the equation numerically, which can be 
     ;MAD less accurate than the equations above for some M and z, because
     ;MAD of gridding in r).
     ;MAD Set NFW scale radius

     IF (n_elements(logm) EQ 1 AND n_elements(z) GT 1) THEN dim1=n_elements(z) ELSE dim1=n_elements(logm)
     IF (n_elements(logm) EQ 1 AND n_elements(z) GT 1) THEN dim2=n_elements(logm) ELSE dim2=n_elements(z)

     r=cgloggen(100000,start=1d-5,finish=10)
     n_r=n_elements(r)
     r=reform(transpose(rebin([r],n_r,n_elements(z),n_elements(logm))))
     
     r_s=rvir/c
     r_s=reform(rebin([r_s],dim1,dim2,n_r))
     c=reform(rebin([c],dim1,dim2,n_r))
     
     ;MAD Get mean matter density
     rho_c=rho_crit(z)
     rho_c=transpose(rebin([rho_c],n_elements(z),n_elements(logm)))
     rho_c=reform(rebin([rho_c],n_elements(logm),n_elements(z),n_r))
     rho_out=rho_crit(z,/physical)*omega_z
     rho_out=transpose(rebin([rho_out],n_elements(z),n_elements(logm)))
     target=reform(rho_out*delta_cr)
     rho_out=reform(rebin([rho_out],n_elements(logm),n_elements(z),n_r))

     delta_cr_v=transpose(rebin([delta_cr_v],n_elements(z),n_elements(logm)))
     delta_cr_v=reform(rebin([delta_cr_v],n_elements(logm),n_elements(z),n_r))
     delta_c=(Delta_cr_v/3.)*((c^3)/((alog(1+c)-(c/(1+c)))))

     term1=(r_s/(r_s+r))+alog(r_s+r)
     term2=1.+alog(r_s)

     lhs=delta_c*((r_s/r)^3.)*(rho_c/rho_out)*(term1-term2)
     rhs=(1/3.)*delta_cr

     IF (n_elements(logm) EQ 1 OR n_elements(z) EQ 1) THEN BEGIN
        n=max([n_elements(logm),n_elements(z)])
        r_deltacr=fltarr(n)
        FOR i=0L,n-1 DO BEGIN
           counter,i,n
           yy=closest(reform(lhs[i,*]),rhs)
           r_deltacr[i]=reform(r[i,yy])
        ENDFOR
     ENDIF ELSE BEGIN
        r_deltacr=fltarr(n_elements(logm),n_elements(z))
        FOR i=0L,n_elements(logm)-1 DO BEGIN
           counter,i,n_elements(logm)
           FOR j=0L,n_elements(z)-1 DO BEGIN
              yy=closest(reform(lhs[i,j,*]),rhs)
              r_deltacr[i,j]=reform(r[i,j,yy])
           ENDFOR
        ENDFOR
     ENDELSE
        
     m_deltacr=(4./3.)*!dpi*(r_deltacr^3.)*target
  ENDELSE

  return,m_deltacr
END
