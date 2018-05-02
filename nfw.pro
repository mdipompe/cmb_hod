FUNCTION nfw,r,logm,z,cindex=cindex,rvir=rvir,$
             c0=c0,beta=beta,mstar=mstar,$
             M_min=M_min,sig_logM=sig_logM,M1=M1,alpha=alpha,$
             project=project,type=type,delta_cr=delta_cr
;+
;  NAME:
;    nfw
;
;  PURPOSE:
;    Generate an NFW (97) density profile
;
;  USE:
;    rho=nfw(r,logm,z[,...])
;
;  INPUT:
;    r - array of radii at which to calculate density 
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
;    M_min - min mass for HOD parameterization.
;    sig_logM - for HOD parameterization
;    M1 - mass normalization for HOD parameterization
;    alpha - power-law index for HOD parameterization
;    delta_cr - the multiplicative factor for defining
;               halos. Defaults to virial parameterization of
;               Hill & Spergel (2014) (amongst many others)
;  
;  OPTIONAL KEYWORDS:
;    project - set to return 2D projected density profile
;    type - Dark matter ('DM'; default) or number ('num') density
;
;  OUTPUT:
;    1, 2, or 3D array of densities. If arrays of both z and logM are
;    supplied, 3D array has values for each mass in first dimension,
;    and values for each z in second dimension, and r in last.
;
;  NOTES:
;    Uses MADs cosmology common block, initialized with load_cosmology.pro
;
;
;  HISTORY:
;    4/30/18 - Written - MAD (Dartmouth)
;-  
  COMMON cosmological_parameters
  
  ;MAD Set default type (dark matter)
  IF (n_elements(type) EQ 0) THEN type='dm'

  ;MAD Some matrix dimension orders are important for working with different
  ;size input combinations of M and z.
  IF (n_elements(logm) EQ 1 AND n_elements(z) GT 1) THEN dim1=n_elements(z) ELSE dim1=n_elements(logm)
  IF (n_elements(logm) EQ 1 AND n_elements(z) GT 1) THEN dim2=n_elements(logm) ELSE dim2=n_elements(z)

  ;MAD Copy r into matrix to do all M and z simultaneously
  IF (size(r,/n_dimensions) EQ 1) THEN BEGIN
     n_r=n_elements(r)
     r2=reform(transpose(rebin([r],n_elements(r),n_elements(z),n_elements(logm))))
  ENDIF ELSE BEGIN
     s=size(r)
     n_r=s[3]
     IF (s[1] NE n_elements(logm) OR s[2] NE n_elements(z)) THEN $
        message,'Dimensions of r not compatible with logm and/or z'
     r2=reform(r)
  ENDELSE
  
  ;MAD Get virial radius if not supplied
  IF (n_elements(rvir) EQ 0) THEN $
     rvir=r_vir(logm,z,delta_cr=delta_cr)
  
  ;MAD Get concentration index if not supplied
  IF (n_elements(cindex) EQ 0) THEN $
     cindex=c_index(logm,z,c0=c0,beta=beta,mstar=mstar)

  ;MAD Set scale radius, put into matrix
  r_s=rvir/cindex
  r_s=reform(rebin([r_s],dim1,dim2,n_r))

  ;MAD If doing DM, get critical density and 
  ;delta_c, normalization for densifty profile
  IF (type EQ 'dm') THEN BEGIN
     rho_crit=rho_crit(z)
     ;MAD Need to expand this into a matrix for multiplication with 
     ;other factors that depend on M and z
     rho_crit=reform(transpose(rebin([rho_crit],n_elements(z),n_elements(logm))))
     ;MAD add another dimension for r values
     rho_crit=reform(rebin([rho_crit],dim1,dim2,n_r))

     IF (n_elements(delta_cr) EQ 0) THEN BEGIN
        ;MAD From Hill & Spergel 14, get nonlinear overdensity for collapse
        omega_z=evolve_density(z,type='matter')
        Delta_cr=(18.*!dpi^2) + (82*(omega_z-1.)) - (39.*(omega_z-1)^2.)
     ENDIF
     
     ;Also needs to be expanded into a matrix for multiplication
     delta_cr=reform(transpose(rebin([delta_cr],n_elements(z),n_elements(logm))))

     ;MAD Equation 2 from NFW modified for evolving critical overdensity
     ;for collapse, rather than just 200.
     delta_c=(delta_cr/3.)*(cindex^3.)/(alog(1.+cindex)-(cindex/(1.+cindex)))
     delta_c=reform(rebin([delta_c],dim1,dim2,n_r))

     ;MAD get 3D rho for DM, first putting everything into the 
     ;appropriately sized data cubes
     IF (n_elements(project) EQ 0) THEN BEGIN
        rho=(rho_crit*delta_c)/((r2/r_s)*((1.+(r2/r_s))^2.))
     ENDIF ELSE BEGIN
        ;MAD Calculate and return projected rho(r) (Bartelmann 1996, eq 7,8)
        fx=fltarr(n_elements(logm),n_elements(z),n_r)
        xx=where(r2/r_s GT 1., cnt)
        IF (cnt NE 0) THEN $
           fx[xx]=1.-((2./sqrt(((r2[xx]/r_s[xx])^2)-1))*atan(sqrt(((r2[xx]/r_s[xx])-1)/((r2[xx]/r_s[xx])+1))))
        xx=where(r2/r_s LT 1., cnt)
        IF (cnt NE 0) THEN $
           fx[xx]=1.-((2./sqrt(1.-(r2[xx]/r_s[xx])^2.))*atanh(sqrt((1.-r2[xx]/r_s[xx])/(1.+r2[xx]/r_s[xx]))))
        xx=where(r/r_s EQ 1., cnt)
        IF (cnt NE 0) THEN $
           fx[xx]=0.

        rho=((2.*rho_crit*delta_c*r_s)/(((r2/r_s)^2.)-1))*fx
     ENDELSE
  ENDIF

  ;MAD If doing number density (e.g. galaxies), get <N(M)> and
  ;set normalization density
  IF (type EQ 'num') THEN BEGIN
     N_m=hod(logm,m_min=m_min,sig_logm=sig_logm,m1=m1,alpha=alpha)
     N_m=reform(rebin([N_m],n_elements(logm),n_elements(z)))
     
     ;MAD Get normalization for NFW profile (set such that integral
     ;of profile gives you the mean number)
     rho_s=(N_m/(4.*!dpi*r_s^3.))*(((1./(cindex+1.))+alog(cindex+1.)-1.)^(-1.))
     rho_s=reform(rebin([rho_s],dim1,dim2,n_r))
   
     
     ;MAD Get 3D rho for number density
     IF (n_elements(project) EQ 0) THEN BEGIN
        rho=rho_s/((r2/r_s)*((1.+(r2/r_s))^2.))
     ENDIF ELSE BEGIN
        ;MAD Calculate and return projected rho(r) (Bartelmann 1996, eq 7,8)
        fx=fltarr(n_elements(logm),n_elements(z),n_r)
        xx=where(r2/r_s GT 1., cnt)
        IF (cnt NE 0) THEN $
           fx[xx]=1.-((2./sqrt(((r2[xx]/r_s[xx])^2)-1))*atan(sqrt(((r2[xx]/r_s[xx])-1)/((r2[xx]/r_s[xx])+1))))
        xx=where(r2/r_s LT 1., cnt)
        IF (cnt NE 0) THEN $
           fx[xx]=1.-((2./sqrt(1.-(r2[xx]/r_s[xx])^2.))*atanh(sqrt((1.-r2[xx]/r_s[xx])/(1.+r2[xx]/r_s[xx]))))
        xx=where(r2/r_s EQ 1., cnt)
        IF (cnt NE 0) THEN $
           fx[xx]=0.

        rho=((2.*rho_s*r_s)/(((r2/r_s)^2.)-1))*fx     
     ENDELSE
  ENDIF

  return, reform(rho)

END
