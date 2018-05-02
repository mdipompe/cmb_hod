PRO go

  r=findgen(1000)/100+0.01
  logm=13.
  z=1.

  h=0.702
  omega_m=0.275
  omega_l=0.725
  
  ;MAD Set default type (dark matter)
  IF (n_elements(type) EQ 0) THEN type='dm'
  
  ;MAD Set default cosmology params
     IF (n_elements(h) EQ 0) THEN h=0.702
     IF (n_elements(omega_m) EQ 0) THEN omega_m=0.275
     IF (n_elements(omega_l) EQ 0) THEN omega_l=0.725

  ;MAD Set some defaults for concentration index
  ;Power index for Mass dependence
  IF (n_elements(beta) EQ 0) THEN beta=-0.13
  ;Normalization of concentration index
  IF (n_elements(c0) EQ 0) THEN c0=11

  ;MAD Get virial radius if not supplied
;  IF (n_elements(r_vir) EQ 0) THEN $
  r_vir=r_vir(logm,z,h=h,omega_m=omega_m,omega_l=omega_l)

  ;MAD Get concentration index
  IF (n_elements(c) EQ 0) THEN BEGIN
     IF (n_elements(mstar) EQ 0) THEN $
        c=c_index(logm,z,h=h,omega_m=omega_m,c0=c0,beta=beta) ELSE $
           c=c_index(logm,z,h=h,omega_m=omega_m,c0=c0,beta=beta,mstar=mstar)
  ENDIF
  
  ;MAD Set scale radius
  r_s=R_vir/c

  ;MAD If doing DM, get critical density and 
  ;delta_c, normalization for densifty profile
  IF (type EQ 'dm') THEN BEGIN
     rho_crit=rho_crit(z,h=h,omega_m=omega_m,omega_l=omega_l,/physical)
     omega_z=(omega_m*(1.+z)^3.)/((omega_m*(1.+z)^3.) + omega_l)
     Delta_cr=(18.*!dpi^2) + (82*(omega_z-1.)) - (39.*(omega_z-1)^2.)
     delta_c=(delta_cr/3.)*(c^3.)/(alog(1.+c)-(c/(1.+c)))
  ENDIF

  
  ;MAD Get rho projected to surface density for DM
  ;MAD Calculate and return projected rho(r) (Bartelmann 1996, eq 7,8)
  fx=fltarr(n_elements(r))
  xx=where(r/r_s GT 1., cnt)
  IF (cnt NE 0) THEN $
     fx[xx]=1.-((2./sqrt(((r[xx]/r_s)^2)-1))*atan(sqrt(((r[xx]/r_s)-1.)/((r[xx]/r_s)+1))))
  xx=where(r/r_s LT 1., cnt)
  IF (cnt NE 0) THEN $
     fx[xx]=1.-((2./sqrt(1.-(r[xx]/r_s)^2.))*atanh(sqrt((1.-r[xx]/r_s)/(1.+r[xx]/r_s))))
  xx=where(r/r_s EQ 1., cnt)
  IF (cnt NE 0) THEN $
     fx[xx]=0.
  
  rho1=((2.*rho_crit*delta_c*r_s)/(((r/r_s)^2.)-1))*fx     

;  tmp=((10.^logm)/(4.*!dpi*(r_vir^3.)))*c^3*(1./(alog(1.+c)-(c/(1.+c))))  
;  rho1=((2.*tmp*r_s)/(((r/r_s)^2.)-1))*fx 
  print,((10.^logm)/((4./3)*!dpi*r_vir^3))/200.

  norm1=(2.*rho_crit*delta_c*r_s)
;  norm1=(2.*tmp*r_s)



  ;Mad instead use Lokas & Mamon 2001
  gc=1./(alog(1+c)-(c/(1.+c)))
  Cinv=fltarr(n_elements(r))
  xx=where(r GT r_s, cnt)
  IF (cnt NE 0) THEN $
     Cinv[xx]=acos(1./(c*r[xx]/r_vir))
  xx=where(r LT r_s, cnt)
  IF (cnt NE 0) THEN $
     Cinv[xx]=acosh(1./(c*r[xx]/r_vir))
  
  term1=((c^2.)*gc/(2*!dpi))
  term2=(10.^logm)/(r_vir^2)
  term3=(1.-(abs((c^2)*((r/r_vir)^2)-1.)^(-0.5))*Cinv)/(((c^2)*(r/r_vir)^2)-1.)

  norm2=term1*term2
  
  rho2=term1*term2*term3
  xx=where(r EQ r_s, cnt)
  IF (cnt NE 0) THEN $
     rho2[xx]=term1*term2*(1./3.)
  

  plot,r,rho1,linestyle=1,/xlog,/ylog,xtit='R',ytit=textoidl('\Sigma(R)'),charsize=2
  oplot,r,rho2,linestyle=2
  oplot,[r_vir,r_vir],[1d1,1d20],linestyle=1
  legend,['L&M01','B96'],linestyle=[2,1],box=0,charsize=1.5,/top,/right
  
  
  stop
END
