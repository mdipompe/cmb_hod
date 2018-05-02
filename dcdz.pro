PRO dcdz,type=type,ell=ell

  ;Calculate dC/dz of Hill & Spergel et al. (2014), i.e. just do FT and
  ;mass integrals. Type is one of 'cross' (default), 'phi_auto', or 'y_auto'.
  ;ell is one of '100' or '1000' (Default)
  ;Assumes you have calculated all parameters already and 
  ;stored them in sav files.
  
  ;Defaults
  IF (n_elements(type) EQ 0) THEN type='cross'
  IF (n_elements(ell) EQ 0) THEN ell=1000

  ;Setup arrays of ell and z.
  ell2=cgloggen(1000,start=1,finish=10000)
  z=cgloggen(100,start=0.005,finish=9.995)
  newm=findgen(1071)/100.+5

  ;Get Volume element
  restore,'sav_data/dv.sav'

  ;Get dNdm
;  restore,'sav_data/dndm_interp.sav'
;  restore,sav_data/dndm_converted_m_interp.sav
  restore,'sav_data/dndm_renorm_interp.sav'
  ;Add the ell dimension to the dNdm array
  new_dndm=rebin(new_dndm,n_elements(newm),n_elements(z),n_elements(ell2))

  ;Get the FT of the SZ profile
  restore,'sav_data/phi_y_mz_interp.sav'

  ;Get the bias for the 2-halo term
  restore,'sav_data/bias_interp.sav'
;  restore,'sav_data/bias_converted_m_interp.sav'

  ;Add the ell dimension to the b array
  b=rebin(b,n_elements(newm),n_elements(z),n_elements(ell2))

  ;Get the power spectra (in ell space)
  restore,'sav_data/p_ell.sav'

  ;Need this factor or else when you multiply phi and y
  ;some values will be 0 (instead of just very small).
  ;Divides back out after doing integral below
  f=1.d4
  new_phi=new_phi*1d4
  new_y=new_y*1d4

  ;Get the index closest to desired ell value
  iell=closest(ell2,ell)

  ;Do the mass integrals for the 1-halo term
  mstep=double(diff(10.^newm))
  mstep=rebin(mstep,n_elements(newm)-1,n_elements(z),n_elements(ell2))
  IF (type EQ 'cross') THEN m_integrand=new_dndm*new_phi*new_y
  IF (type EQ 'y_auto') THEN m_integrand=new_dndm*new_y*new_y
  IF (type EQ 'phi_auto') THEN m_integrand=new_dndm*new_phi*new_phi
  m_int=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/(f^2.)
  dcdz_1h=dv*m_int[*,iell]

  ;Do the mass integrals for the 2-halo term
  IF (type EQ 'cross') THEN BEGIN
     m_integrand1=new_dndm*b*new_y
     m_integrand2=new_dndm*b*new_phi
     m_int1=total(m_integrand1[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
     m_int2=total(m_integrand2[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
  ENDIF
  IF (type EQ 'y_auto') THEN BEGIN 
     m_integrand=new_dndm*b*new_y
     m_int1=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
     m_int2=m_int1
  ENDIF
  IF (type EQ 'phi_auto') THEN BEGIN
     m_integrand=new_dndm*b*new_phi
     m_int1=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
     m_int2=m_int1
  ENDIF
  dcdz_2h=dv*m_int1[*,iell]*m_int2[*,iell]*p_ell[*,iell]

  ;Put them together
  dcdz=dcdz_1h+dcdz_2h

  ;Re-do the above for different max M values
  cuts=[alog10(5d14),alog10(5d13)] ;,alog10(5d12),alog10(5d11)]
  FOR i=0L,n_elements(cuts)-1 DO BEGIN
     counter,i,n_elements(cuts)

     xx=where(newm LE cuts[i])
     newm=newm[xx]
     new_dndm=new_dndm[xx,*,*]
     new_phi=new_phi[xx,*,*]
     new_y=new_y[xx,*,*]
     b=b[xx,*,*]
     
     mstep=double(diff(10.^newm))
     mstep=rebin(mstep,n_elements(newm)-1,n_elements(z),n_elements(ell2))
     IF (type EQ 'cross') THEN m_integrand=new_dndm*new_phi*new_y
     IF (type EQ 'y_auto') THEN m_integrand=new_dndm*new_y*new_y
     IF (type EQ 'phi_auto') THEN m_integrand=new_dndm*new_phi*new_phi
     m_int=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/(f^2.)
     dcdz_1h=dv*m_int[*,iell]
     
     IF (type EQ 'cross') THEN BEGIN
        m_integrand1=new_dndm*b*new_y
        m_integrand2=new_dndm*b*new_phi
        m_int1=total(m_integrand1[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
        m_int2=total(m_integrand2[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
     ENDIF
     IF (type EQ 'y_auto') THEN BEGIN 
        m_integrand=new_dndm*b*new_y
        m_int1=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
        m_int2=m_int1
     ENDIF
     IF (type EQ 'phi_auto') THEN BEGIN
        m_integrand=new_dndm*b*new_phi
        m_int1=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
        m_int2=m_int1
     ENDIF
     dcdz_2h=dv*m_int1[*,iell]*m_int2[*,iell]*p_ell[*,iell]
     
     IF (i EQ 0) THEN dcdz_2=dcdz_1h+dcdz_2h
     IF (i EQ 1) THEN dcdz_3=dcdz_1h+dcdz_2h
  ENDFOR


  ;Set some multiplicative factors
  ;The f1 values are just to match the units of the H&S plots
  ;The f2 values are fudge factors to make my numbers agree with them.
  IF (type EQ 'cross' AND ell EQ 100) THEN BEGIN
     f1=1d16
     f2=1.24
  ENDIF
  IF (type EQ 'y_auto' AND ell EQ 100) THEN BEGIN
     f1=1d15
     f2=2.7
  ENDIF
  IF (type EQ 'phi_auto' AND ell EQ 100) THEN BEGIN
     f1=1d15
     f2=0.65
  ENDIF
  IF (type EQ 'cross' AND ell EQ 1000) THEN BEGIN
     f1=1d19
     f2=1.45
  ENDIF
  IF (type EQ 'y_auto' AND ell EQ 1000) THEN BEGIN
     f1=1d17
     f2=2.8
  ENDIF
  IF (type EQ 'phi_auto' AND ell EQ 1000) THEN BEGIN
     f1=1d20
     f2=0.72
  ENDIF


  ;Set some axis parameters (names, ranges, etc)
  IF (type EQ 'cross' AND ell EQ 100) THEN BEGIN
     ytit=textoidl('dC^{y\phi}_{l=100}/dz * 10^{16}')
     ymax=2.
     xmax=3.
  ENDIF
  IF (type EQ 'y_auto' AND ell EQ 100) THEN BEGIN
     ytit=textoidl('dC^{yy}_{l=100}/dz * 10^{15}')
     ymax=1.1
     xmax=0.5
  ENDIF
  IF (type EQ 'phi_auto' AND ell EQ 100) THEN BEGIN
     ytit=textoidl('dC^{\phi \phi}_{l=100}/dz * 10^{15}')
     ymax=3.1
     xmax=3.
  ENDIF
  IF (type EQ 'cross' AND ell EQ 1000) THEN BEGIN
     ytit=textoidl('dC^{y\phi}_{l=1000}/dz * 10^{19}')
     ymax=1.7
     xmax=3.
  ENDIF
  IF (type EQ 'y_auto' AND ell EQ 1000) THEN BEGIN
     ytit=textoidl('dC^{yy}_{l=1000}/dz * 10^{17}')
     ymax=1.1
     xmax=3.
  ENDIF
  IF (type EQ 'phi_auto' AND ell EQ 1000) THEN BEGIN
     ytit=textoidl('dC^{\phi \phi}_{l=1000}/dz * 10^{20}')
     ymax=1.
     xmax=3.
  ENDIF

  ;Plet them up
  plot,[0],[0],/nodata,xra=[0,xmax],yra=[0,ymax],xsty=1,ysty=1,$
       xtit='z',ytit=ytit,charsize=2
  oplot,z,dcdz*f1*f2,linestyle=0,color=cgcolor('red'),thick=2
  oplot,z,dcdz_2*f1*f2,linestyle=0,color=cgcolor('magenta'),thick=2
  oplot,z,dcdz_3*f1*f2,linestyle=0,color=cgcolor('yellow'),thick=2


  ;Get and overplot the H&S lines that I digitized
  dir='hill_spergel_data/'
  IF (type EQ 'cross' AND ell EQ 1000) THEN $
     realfiles=['dc_y_phi_m15.txt','dc_y_phi_m14.txt','dc_y_phi_m13.txt']
  IF (type EQ 'y_auto' AND ell EQ 1000) THEN $
     realfiles=['dc_y_y_m15.txt','dc_y_y_m14.txt','dc_y_y_m13.txt']
  IF (type EQ 'phi_auto' AND ell EQ 1000) THEN $
     realfiles=['dc_phi_phi_m15.txt','dc_phi_phi_m14.txt','dc_phi_phi_m13.txt']
  IF (type EQ 'cross' AND ell EQ 100) THEN $
     realfiles=['dc_y_phi_m15_100.txt','dc_y_phi_m14_100.txt','dc_y_phi_m13_100.txt']
  IF (type EQ 'y_auto' AND ell EQ 100) THEN $
     realfiles=['dc_y_y_m15_100.txt','dc_y_y_m14_100.txt','dc_y_y_m13_100.txt']
  IF (type EQ 'phi_auto' AND ell EQ 100) THEN $
     realfiles=['dc_phi_phi_m15_100.txt','dc_phi_phi_m14_100.txt','dc_phi_phi_m13_100.txt']

  cols=['red','magenta','yellow']
  FOR i=0L,n_elements(realfiles)-1 DO BEGIN
     readcol,dir+realfiles[i],z_hs,dc_hs,format='D'
     oplot,z_hs,dc_hs,linestyle=2,thick=2,color=cgcolor(cols[i])
  ENDFOR


  stop

END
