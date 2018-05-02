PRO two_halo

  ;Get the 2-halo C(ell), after calculating all the individual
  ;peices and making .sav files.
  ;Hill & Spergel (2014) Equation 11.
  
  COMMON cosmological_parameters

  ell=cgloggen(1000,start=1,finish=10000)
  z=cgloggen(100,start=0.005,finish=9.995)
  newz=cgloggen(10000,start=0.005,finish=9.995)
  logm=findgen(108)/10.+5
  newm=findgen(1071)/100.+5
  pspec_ell=file_search('power_spec/ell/camb_matterpower*')

  ;Get dNdM, add ell dimension
  restore,'sav_data/dndm_renorm_interp.sav'
  new_dndm=rebin(new_dndm,n_elements(newm),n_elements(z),n_elements(ell))

  ;Get phi_ell, y_ell
  restore,'sav_data/phi_y_mz_interp.sav'

  ;Get dV, add ell dimension
  restore,'sav_data/dv.sav'
  dv=rebin(dv,n_elements(z),n_elements(ell))

  restore,'sav_data/dv_interp.sav'
  new_dv=rebin(new_dv,n_elements(newz),n_elements(ell))

  ;Read in power spectra, put into array
  p_ell=dblarr(n_elements(z),n_elements(ell))
  FOR i=0L,n_elements(pspec_ell)-1 DO BEGIN
     readcol,pspec_ell[i],l,p,format='D',/silent
     p_ell[i,*]=p
  ENDFOR
  save,p_ell,filename='sav_data/p_ell.sav'

  ;Get bias factor
  restore,'sav_data/bias_interp.sav'
  b=rebin(b,n_elements(newm),n_elements(z),n_elements(ell))
  
 
  ;Do M integral
  mstep=diff(10.^newm)
  mstep=rebin(mstep,n_elements(newm)-1,n_elements(z),n_elements(ell))
  m_integrand_phi=new_dndm*b*new_phi
  m_integrand_y=new_dndm*b*new_y
  m_int_phi=total(m_integrand_phi[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)
  m_int_y=total(m_integrand_y[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)

  ;Interpolate z dimension to finer grid
  new_m_int_phi=findgen(n_elements(newz),n_elements(ell))
  new_m_int_y=findgen(n_elements(newz),n_elements(ell))
  new_p_ell=findgen(n_elements(newz),n_elements(ell))
  FOR i=0L,n_elements(ell)-1 DO BEGIN
     counter,i,n_elements(ell)
     new_m_int_phi[*,i]=spline(z,m_int_phi[*,i],newz)
     new_m_int_y[*,i]=spline(z,m_int_y[*,i],newz)
     new_p_ell[*,i]=spline(z,p_ell[*,i],newz)
  ENDFOR

  ;Do z integral
  zstep=diff(newz)
  zstep=rebin(zstep,n_elements(newz)-1,n_elements(ell))
  z_integrand=new_dv*new_m_int_phi*new_m_int_y*new_p_ell
  c_l_2h=total(z_integrand[1:n_elements(newz)-1,*]*zstep[*,*],1,/double)

  zstep2=diff(z)
  zstep2=rebin(zstep2,n_elements(z)-1,n_elements(ell))
  z_integrand2=dv*m_int_phi*m_int_y*p_ell
  c_l_2h_2=total(z_integrand2[1:n_elements(z)-1,*]*zstep2[*,*],1,/double)

  ;Output
  save,c_l_2h,filename='sav_data/c_l_2h.sav'
  
  stop
END
