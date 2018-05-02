PRO dcdz,type=type,ell=ell

  IF (n_elements(type) EQ 0) THEN type='cross'
  IF (n_elements(ell) EQ 0) THEN ell=1000
  
  ell2=cgloggen(1000,start=1,finish=10000)
  z=cgloggen(100,start=0.005,finish=9.995)
  newm=findgen(1071)/100.+5
  
  restore,'../dv.sav'
  restore,'../dndm_interp.sav'
  new_dndm=rebin(new_dndm,n_elements(newm),n_elements(z),n_elements(ell2))
  restore,'../phi_y_mz_interp.sav'
  undefine,new_y
  restore,'phi_gal_mz_interp.sav'
  restore,'../bias_interp.sav'
  b=rebin(b,n_elements(newm),n_elements(z),n_elements(ell2))
  restore,'../p_ell.sav'

  f=1.d4
  new_phi=new_phi*f
  new_phi_gal=new_phi_gal*f

  iell=closest(ell2,ell)
  
  mstep=double(diff(10.^newm))
  mstep=rebin(mstep,n_elements(newm)-1,n_elements(z),n_elements(ell2))
  IF (type EQ 'cross') THEN m_integrand=new_dndm*new_phi*new_phi_gal
  IF (type EQ 'phi_gal_auto') THEN m_integrand=new_dndm*new_phi_gal*new_phi_gal
  IF (type EQ 'phi_auto') THEN m_integrand=new_dndm*new_phi*new_phi
  m_int=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/(f^2.)
  dcdz_1h=dv*m_int[*,iell]

  IF (type EQ 'cross') THEN BEGIN
     m_integrand1=new_dndm*b*new_phi_gal
     m_integrand2=new_dndm*b*new_phi
     m_int1=total(m_integrand1[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
     m_int2=total(m_integrand2[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
  ENDIF
  IF (type EQ 'phi_gal_auto') THEN BEGIN 
     m_integrand=new_dndm*b*new_phi_gal
     m_int1=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
     m_int2=m_int1
  ENDIF
  IF (type EQ 'phi_auto') THEN BEGIN
     m_integrand=new_dndm*b*new_phi
     m_int1=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
     m_int2=m_int1
  ENDIF
  dcdz_2h=dv*m_int1[*,iell]*m_int2[*,iell]*p_ell[*,iell]

  dcdz=dcdz_1h+dcdz_2h

  cuts=[alog10(5d14),alog10(5d13)] ;,alog10(5d12),alog10(5d11)]
  FOR i=0L,n_elements(cuts)-1 DO BEGIN
     counter,i,n_elements(cuts)

     xx=where(newm LE cuts[i])
     newm=newm[xx]
     new_dndm=new_dndm[xx,*,*]
     new_phi=new_phi[xx,*,*]
     new_phi_gal=new_phi_gal[xx,*,*]
     b=b[xx,*,*]
     
     mstep=double(diff(10.^newm))
     mstep=rebin(mstep,n_elements(newm)-1,n_elements(z),n_elements(ell2))
     IF (type EQ 'cross') THEN m_integrand=new_dndm*new_phi*new_phi_gal
     IF (type EQ 'phi_gal_auto') THEN m_integrand=new_dndm*new_phi_gal*new_phi_gal
     IF (type EQ 'phi_auto') THEN m_integrand=new_dndm*new_phi*new_phi
     m_int=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/(f^2.)
     dcdz_1h=dv*m_int[*,iell]
     
     IF (type EQ 'cross') THEN BEGIN
        m_integrand1=new_dndm*b*new_phi_gal
        m_integrand2=new_dndm*b*new_phi
        m_int1=total(m_integrand1[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
        m_int2=total(m_integrand2[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/f
     ENDIF
     IF (type EQ 'phi_gal_auto') THEN BEGIN 
        m_integrand=new_dndm*b*new_phi_gal
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

  IF (type EQ 'cross' AND ell EQ 100) THEN BEGIN
     f1=1d16
     f2=1.
  ENDIF
  IF (type EQ 'phi_gal_auto' AND ell EQ 100) THEN BEGIN
     f1=1d15
     f2=1.
  ENDIF
  IF (type EQ 'phi_auto' AND ell EQ 100) THEN BEGIN
     f1=1d15
     f2=1.
  ENDIF
  IF (type EQ 'cross' AND ell EQ 1000) THEN BEGIN
     f1=1d21
     f2=1.
  ENDIF
  IF (type EQ 'phi_gal_auto' AND ell EQ 1000) THEN BEGIN
     f1=1d17
     f2=1.
  ENDIF
  IF (type EQ 'phi_auto' AND ell EQ 1000) THEN BEGIN
     f1=1d20
     f2=1.
  ENDIF


  
  IF (type EQ 'cross' AND ell EQ 100) THEN BEGIN
     ytit=textoidl('dC^{\phi_{g} \phi}_{l=100}/dz * 10^{16}')
     ymax=2.
     xmax=3.
  ENDIF
  IF (type EQ 'phi_gal_auto' AND ell EQ 100) THEN BEGIN
     ytit=textoidl('dC^{\phi_{g} \phi_g}_{l=100}/dz * 10^{15}')
     ymax=1.1
     xmax=0.5
  ENDIF
  IF (type EQ 'phi_auto' AND ell EQ 100) THEN BEGIN
     ytit=textoidl('dC^{\phi \phi}_{l=100}/dz * 10^{15}')
     ymax=3.1
     xmax=3.
  ENDIF
  IF (type EQ 'cross' AND ell EQ 1000) THEN BEGIN
     ytit=textoidl('dC^{\phi_{g} \phi}_{l=1000}/dz * 10^{21}')
     ymax=2.5
     xmax=3.
  ENDIF
  IF (type EQ 'phi_gal_auto' AND ell EQ 1000) THEN BEGIN
     ytit=textoidl('dC^{\phi_{g} \phi_{g}}_{l=1000}/dz * 10^{21}')
     ymax=2.7
     xmax=3.
  ENDIF
  IF (type EQ 'phi_auto' AND ell EQ 1000) THEN BEGIN
     ytit=textoidl('dC^{\phi \phi}_{l=1000}/dz * 10^{20}')
     ymax=1.
     xmax=3.
  ENDIF
    
  plot,[0],[0],/nodata,xra=[0,xmax],yra=[0,ymax],xsty=1,ysty=1,$
       xtit='z',ytit=ytit,charsize=2
  oplot,z,dcdz*f1*f2,linestyle=0,color=cgcolor('red'),thick=2
  oplot,z,dcdz_2*f1*f2,linestyle=0,color=cgcolor('magenta'),thick=2
  oplot,z,dcdz_3*f1*f2,linestyle=0,color=cgcolor('yellow'),thick=2

  stop
END
