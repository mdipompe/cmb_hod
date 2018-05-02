;+
;  NAME:
;    cosmocalc
;  PURPOSE:
;    Get cosmological distances, etc.
;
;  USE:
;    res=cosmocalc(z)
;
;  INPUT:
;    z - redshift
;
;  OPTIONAL INPUT:
;    
;  KEYWORDS:
;        
;  OUTPUT:
;    structure with 8 tags: d_h - the hubble distance                      
;                           d_m - the transverse comiving distance (Mpc)   
;                           d_L - the luminosity distance (Mpc)            
;                           d_a - the angular size distnace (Mpc)          
;                           d_c - the comoving distance (Mpc)              
;                           v_c - the comoving volume (Mpc^3)              
;                           t_l - the lookback time (Gyrs)                 
;                           t_h - the hubble time (Gyrs)                   
;  HISTORY:
;    2013 - Written - MAD (UWyo)
;    2015 - Cleaned and documented - MAD (UWyo)
;    1-25-16 - Fixed z=0 Bug - MAD (Dartmouth)
;    12-1-17 - Updated to use common block set in load_cosmology.pro
;-
FUNCTION cosmocalc,z,d_h=d_h,d_l=d_l,d_a=d_a,d_c=d_c,v_c=v_c,t_l=t_l

COMMON cosmological_parameters

;MAD Set speed of light, in km/s
c=2.99792458E5

;MAD Get curvature, set H0
omega_k=1.-omega_m-omega_l
H0=h*100.

;MAD Convert hubble constant to units of s^-1, calculate hubble distance
;(Mpc) and hubble time (yrs)
H0_conv=H0*(1./(3.0859E19))
d_h=(c/H0_conv)*(1./3.0859E19)
t_h=(1./H0_conv)*(1./3600.)*(1./24.)*(1./365.)

;MAD IF z=0, then all distances are 0.
IF (max(z) EQ 0.) THEN BEGIN
   d_m=0.
   d_L=0.
   d_a=0.
   d_c=0.
   V_c=0.
   t_l=0.
ENDIF ELSE BEGIN    ;MAD If z NE 0, integrate.
   ;MAD Integrate
   IF (max(z) LT 10.) THEN zvals=dindgen(max(z)*100000.)/100000. ELSE zvals=dindgen(max(z)*10000.)/10000.
   E=1./SQRT(omega_m*((1+zvals)^(3.0))+Omega_k*((1+zvals)^(2.0))+omega_l)
   E2=1./((1.+zvals)*SQRT(omega_m*((1+zvals)^(3.0))+Omega_k*((1+zvals)^(2.0))+Omega_l))

   y=fltarr(n_elements(z))
   y2=y
   FOR i=0L,n_elements(y)-1 DO BEGIN
      y[i]=int_tabulated(zvals[where(zvals LE z[i])],E[where(zvals LE z[i])],/double)
      y2[i]=int_tabulated(zvals[where(zvals LE z[i])],E2[where(zvals LE z[i])],/double)
   ENDFOR
   
   ;MAD Find comoving distance
   d_c=d_h*y

   ;MAD Set d_m
   IF Omega_k GT 0 THEN d_m=d_h*(1.0/SQRT(Omega_k))*SINH(SQRT(Omega_k)*(d_c/d_h))
   IF Omega_k EQ 0 THEN d_m=d_c
   IF Omega_k LT 0 THEN d_m=d_h*(1.0/SQRT((-1)*Omega_k))*SIN(SQRT((-1)*Omega_k)*(d_c/d_h))

   ;MAD Calculate angular size distance
   d_a=d_m/(1.+z)

   ;MAD Calculate luminosity distance
   d_L=d_m*(1.+z)

   ;MAD Calculate comoving volume within object
   IF Omega_k GT 0 THEN v_c=((4.*!dpi*d_h^3.)/(2.*omega_k))*(((d_m/d_h)*SQRT(1.+omega_k*(d_m^2./d_h^2)))-((1./SQRT(abs(omega_k)))*ASINH(SQRT(abs(omega_k))*(d_m/d_h))))
   IF Omega_k EQ 0 THEN v_c=((4.*!dpi)/3.)*d_m^3.
   IF Omega_k LT 0 THEN v_c=((4.*!dpi*d_h^3.)/(2.*omega_k))*(((d_m/d_h)*SQRT(1.+omega_k*(d_m^2./d_h^2)))-((1./SQRT(abs(omega_k)))*ASIN(SQRT(abs(omega_k))*(d_m/d_h))))

   ;MAD Calculate lookback time
   t_l=t_h*y2*(1e-9)
ENDELSE

;MAD Build return structure
dists={d_h:0.D, d_m:0.D, d_L:0.D, d_a:0.D, d_c:0.D, v_c:0.D, t_l:0.D, t_h:0.D}
dists=replicate(dists,n_elements(z))
dists.d_h=d_h
dists.d_m=d_m
dists.d_L=d_L
dists.d_a=d_a
dists.d_c=d_c
dists.v_c=V_c
dists.t_l=t_l
dists.t_h=t_h*(1e-9)

return,dists
END


