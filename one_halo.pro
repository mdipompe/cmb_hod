PRO one_halo

  ;Calculate the one_halo term, from Hill & Spergel (2014)
  ;equation 8. Need to have created all the necessary
  ;.sav files for all the parameters first

  ;Also has some tests for checking integrals
  
  COMMON cosmological_parameters

  ell=cgloggen(1000,start=1,finish=10000)
  z=cgloggen(100,start=0.005,finish=9.995)
  logm=findgen(108)/10.+5
  newm=findgen(1071)/100.+5
  newz=cgloggen(10000,start=0.005,finish=9.995)

  ;Get dV, add ell dimension
  restore,'sav_data/dv.sav'
  dv=rebin(dv,n_elements(z),n_elements(ell))

  ;Get interpolated dv, add ell dimension
  restore,'sav_data/dv_interp.sav'
  new_dv=rebin(new_dv,n_elements(newz),n_elements(ell))

  ;Get dndm, add ell dimension
  restore,'sav_data/dndm_renorm.sav'
  dndm=rebin(dndm,n_elements(logm),n_elements(z),n_elements(ell))

  ;Get interpolated dndm, add ell dimension
  restore,'sav_data/dndm_renorm_interp.sav'
  new_dndm=rebin(new_dndm,n_elements(newm),n_elements(z),n_elements(ell))

  ;Get phi_ell, y_ell
  restore,'sav_data/phi_y_mz.sav'
  restore,'sav_data/phi_y_mz_interp.sav'

  ;Need a multiplicative factor so you don't get 0s when combining
  ;very small terms. Divided back out after integrating
  f=1d4  
  phi_mz=phi_mz*f
  y_mz=y_mz*f
  new_phi=new_phi*f
  new_y=new_y*f

  ;Do the M integral using fine grid
  mstep=diff(10.^newm)
  mstep=rebin(mstep,n_elements(newm)-1,n_elements(z),n_elements(ell))
  m_integrand=new_dndm*new_phi*new_y
  m_int=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/(f^2.)

  ;Do the M integral using course grid
  mstep2=diff(10.^logm)
  mstep2=rebin(mstep2,n_elements(logm)-1,n_elements(z),n_elements(ell))
  m_integrand2=dndm*phi_mz*y_mz
  m_int2=total(m_integrand2[1:n_elements(logm)-1,*,*]*mstep2[*,*,*],1,/double)/(f^2)

  ;Do the M integral using int_tabulated
  ell_test=100
  z_test=0
  test=int_tabulated(10.^newm,m_integrand[*,z_test,ell_test],/double)/(f^2.)
  test2=int_tabulated(10.^logm,m_integrand2[*,z_test,ell_test],/double)/(f^2.)

  ;Compare integration methods
  print,'Mass integral:'
  print,'Interpolated int_tab: ',test
  print,'Interpolated simple: ',m_int[z_test,ell_test]
  print,'Original int_tab: ',test2
  print,'Original simple: ',m_int2[z_test,ell_test]

  stop
  
  ;Interpolate M integral to finer grid in z
  new_m_int=findgen(n_elements(newz),n_elements(ell))
  FOR i=0L,n_elements(ell)-1 DO BEGIN
     counter,i,n_elements(ell)
     new_m_int[*,i]=spline(z,m_int[*,i],newz)     
  ENDFOR

  ;Do z integral with fine grid
  zstep=diff(newz)
  zstep=rebin(zstep,n_elements(newz)-1,n_elements(ell))
  z_integrand=new_dv*new_m_int
  c_l_1h=total(z_integrand[1:n_elements(newz)-1,*]*zstep[*,*],1,/double)

  
  save,c_l_1h,filename='sav_data/c_l_1h.sav'

  ;Do z integral with course grid
  zstep2=diff(z)
  zstep2=rebin(zstep2,n_elements(z)-1,n_elements(ell))
  z_integrand2=dv*m_int
  c_l2=total(z_integrand2[1:n_elements(z)-1,*]*zstep2[*,*],1,/double)

  ;Do z integrals with int_tabulated
  test=int_tabulated(newz,z_integrand[*,ell_test],/double)
  test2=int_tabulated(z,z_integrand2[*,ell_test],/double)

  ;Compare methods
  print,'z integral:'
  print,'Interpolated int_tab: ',test
  print,'Interpolated simple: ',c_l_1h[ell_test]
  print,'Original int_tab: ',test2
  print,'Original simple: ',c_l2[ell_test]

  stop
END
