PRO fig_5,type=type,ell=ell

  ;Reproduce fig 5 of Hill & Spergel (2014). Can specify
  ;types of 'cross' (default), 'phi_auto', and 'y_auto'.
  ;Can take 2 ell values, either 100 or 1000 (default)
  
  IF (n_elements(type) EQ 0) THEN type='cross'
  IF (n_elements(ell) EQ 0) THEN ell=1000
  
  ell2=cgloggen(1000,start=1,finish=10000)
  z=cgloggen(100,start=0.005,finish=9.995)
  newm=findgen(1071)/100.+5
  newz=cgloggen(10000,start=0.005,finish=9.995)

  ;Get dV, add ell dimension
  restore,'sav_data/dv_interp.sav'
  new_dv=rebin(new_dv,n_elements(newz),n_elements(ell2))

  ;Get dNdM, add ell dimension
;  restore,'sav_data/dndm_interp.sav'
  restore,'sav_data/dndm_renorm_interp.sav'
;  restore,'sav_data/dndm_hmfcalc.sav'  
  new_dndm=rebin(new_dndm,n_elements(newm),n_elements(z),n_elements(ell2))

  ;Get phi_ell and y_ell
  restore,'sav_data/phi_y_mz_interp.sav'

  ;Get bias, add ell dimension
  restore,'sav_data/bias_interp.sav'
;  restore,'sav_data/bias_converted_m_interp.sav'
  b=rebin(b,n_elements(newm),n_elements(z),n_elements(ell2))

  ;Get P(ell,z)
  restore,'sav_data/p_ell.sav'

  ;Need this factor for when you multiply terms together or
  ;you'll get 0s. Divided back out later
  f=1.d4
  new_phi=new_phi*1d4
  new_y=new_y*1d4

  ;Get the ell value closest to your input
  iell=closest(ell2,ell)

  ;Get baseline model (all M)
  xx=where(newm LE alog10(5d15))
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

  ;Now that mass integral is done, interpolate over z
  new_m_int=findgen(n_elements(newz),n_elements(ell2))
  FOR i=0L,n_elements(ell2)-1 DO BEGIN
     counter,i,n_elements(ell2)
     new_m_int[*,i]=spline(z,m_int[*,i],newz)     
  ENDFOR

  ;Integrate z to get 1-halo term
  zstep=diff(newz)
  zstep=rebin(zstep,n_elements(newz)-1,n_elements(ell2))
  z_integrand=new_dv*new_m_int
  c_l_1h=total(z_integrand[1:n_elements(newz)-1,*]*zstep[*,*],1,/double)

  ;Now do integrals for 2-halo term
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
  
  new_m_int1=findgen(n_elements(newz),n_elements(ell2))
  new_m_int2=findgen(n_elements(newz),n_elements(ell2))
  new_p_ell=findgen(n_elements(newz),n_elements(ell2))
  FOR i=0L,n_elements(ell2)-1 DO BEGIN
     counter,i,n_elements(ell2)
     new_m_int1[*,i]=spline(z,m_int1[*,i],newz)
     new_m_int2[*,i]=spline(z,m_int2[*,i],newz)
     new_p_ell[*,i]=spline(z,p_ell[*,i],newz)
  ENDFOR
  z_integrand=new_dv*new_m_int1[*,*]*new_m_int2[*,*]*new_p_ell[*,*]
  c_l_2h=total(z_integrand[1:n_elements(newz)-1,*]*zstep[*,*],1,/double)

  ;Put 1 and 2 halo terms together
  cl=c_l_1h[iell]+c_l_2h[iell]



  ;Now re-do everything for different mass cuts
  cuts=[alog10(1d15),alog10(1d14),alog10(1d13),alog10(1d12)] ;,alog10(5d12),alog10(5d11)]
  FOR j=0L,n_elements(cuts)-1 DO BEGIN
     counter,j,n_elements(cuts)

     xx=where(newm LE cuts[j])
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

     new_m_int=findgen(n_elements(newz),n_elements(ell2))
     FOR i=0L,n_elements(ell2)-1 DO BEGIN
;        counter,i,n_elements(ell2)
        new_m_int[*,i]=spline(z,m_int[*,i],newz)     
     ENDFOR
     
     zstep=diff(newz)
     zstep=rebin(zstep,n_elements(newz)-1,n_elements(ell2))
     z_integrand=new_dv*new_m_int
     c_l_1h=total(z_integrand[1:n_elements(newz)-1,*]*zstep[*,*],1,/double)

     

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

     new_m_int1=findgen(n_elements(newz),n_elements(ell2))
     new_m_int2=findgen(n_elements(newz),n_elements(ell2))
     new_p_ell=findgen(n_elements(newz),n_elements(ell2))
     FOR i=0L,n_elements(ell2)-1 DO BEGIN
;        counter,i,n_elements(ell2)
        new_m_int1[*,i]=spline(z,m_int1[*,i],newz)
        new_m_int2[*,i]=spline(z,m_int2[*,i],newz)
        new_p_ell[*,i]=spline(z,p_ell[*,i],newz)
     ENDFOR
     z_integrand=new_dv*new_m_int1[*,*]*new_m_int2[*,*]*new_p_ell[*,*]
     c_l_2h=total(z_integrand[1:n_elements(newz)-1,*]*zstep[*,*],1,/double)
     
     IF (j EQ 0) THEN cl_1=c_l_1h[iell]+c_l_2h[iell]
     IF (j EQ 1) THEN cl_2=c_l_1h[iell]+c_l_2h[iell]
     IF (j EQ 2) THEN cl_3=c_l_1h[iell]+c_l_2h[iell]
     IF (j EQ 3) THEN cl_4=c_l_1h[iell]+c_l_2h[iell]
  ENDFOR
  cls=[cl_1,cl_2,cl_3,cl_4]


  dir='hill_spergel_data/'
  IF (type EQ 'cross' and ell EQ 100) THEN realfile='cl_Mmax_cross_100.txt'
  IF (type EQ 'y_auto' and ell EQ 100) THEN realfile='cl_Mmax_yy_100.txt'
  IF (type EQ 'phi_auto' and ell EQ 100) THEN realfile='cl_Mmax_phiphi_100.txt'
  IF (type EQ 'cross' and ell EQ 1000) THEN realfile='cl_Mmax_cross_1000.txt'
  IF (type EQ 'y_auto' and ell EQ 1000) THEN realfile='cl_Mmax_yy_1000.txt'
  IF (type EQ 'phi_auto' and ell EQ 1000) THEN realfile='cl_Mmax_phiphi_1000.txt'
  


  IF (ell EQ 100) THEN ytit=textoidl('C_{l=100}(<M_{max}) / C_{l=100}(total)')
  IF (ell EQ 1000) THEN ytit=textoidl('C_{l=1000}(<M_{max}) / C_{l=1000}(total)')
  plot,[0],[0],/nodata,charsize=2,xra=[1d11,5d15],yra=[0,1.],xsty=1,ysty=1,/xlog,$
       xtit=textoidl('M_{max}'),ytit=ytit
  oplot,10.^cuts,cls/cl,psym=6,symsize=2

  readcol,dir+realfile,x,y,format='D'
  oplot,x,y,linestyle=2,thick=2

  quadterp,x,y,10.^cuts,vals
  oplot,10.^cuts,vals,psym=1,symsize=2
  
  FOR i=0L,n_elements(cuts)-1 DO BEGIN
     print,'Mine/H&S at M < '+strtrim(cuts[i],2)+': '+strtrim((cls[i]/cl)/vals[i],2)
  ENDFOR
  
  
  stop
  
END
