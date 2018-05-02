PRO int_logm_z

  ;This is just comparing the accuracy of integrals,
  ;here using dNdlogM.
  
  COMMON cosmological_parameters

  ell=cgloggen(1000,start=1,finish=10000)
  z=cgloggen(100,start=0.005,finish=9.995)
  logm=findgen(108)/10.+5
  newm=findgen(1071)/100.+5

  restore,'sav_data/dndlogm_interp.sav'
  restore,'sav_data/dndlogm.sav'
  restore,'sav_data/phi_y_mz.sav'
  restore,'sav_data/phi_y_mz_interp.sav'

  new_dndm=rebin(new_dndlogm,n_elements(newm),n_elements(z),n_elements(ell))

  ;Need this factor to avoid getting 0s when you multiply 
  ;very small numbers together. Factored back out later.
  f=1d4  
  phi_mz=phi_mz*f
  y_mz=y_mz*f
  new_phi=new_phi*f
  new_y=new_y*f

  ;Do the M integral by calculating area using fine bins
  mstep=diff(newm)
  mstep=rebin(mstep,n_elements(newm)-1,n_elements(z),n_elements(ell))
  m_integrand=new_dndm*new_phi*new_y
  m_int=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)

  ;Do the M integral by calculating area using broader bins
  mstep2=diff(logm)
  mstep2=rebin(mstep2,n_elements(logm)-1,n_elements(z),n_elements(ell))
  m_integrand2=rebin(dndlogm,n_elements(logm),n_elements(z),n_elements(ell))*phi_mz*y_mz
  m_int2=total(m_integrand2[1:n_elements(logm)-1,*,*]*mstep2[*,*,*],1,/double)

  ;Do the M integral with int_tabulated for fine, broad bins
  ell_test=100
  z_test=0
  test=int_tabulated(newm,m_integrand[*,z_test,ell_test],/double)
  test2=int_tabulated(logm,m_integrand2[*,z_test,ell_test],/double)
  ;Compare the different methods
  print,'Mass integral:'
  print,'Interpolated int_tab: ',test
  print,'Interpolated simple: ',m_int[z_test,ell_test]
  print,'Original int_tab: ',test2
  print,'Original simple: ',m_int2[z_test,ell_test]

  stop


  ;Make finer z grid, interpolate
  newz=cgloggen(10000,start=0.005,finish=9.995)
  new_m_int=findgen(n_elements(newz),n_elements(ell))
  FOR i=0L,n_elements(ell)-1 DO BEGIN
     counter,i,n_elements(ell)
     new_m_int[*,i]=spline(z,m_int[*,i],newz)     
  ENDFOR

  ;Restore volume element
  restore,'sav_data/dv_interp.sav'
  new_dv=rebin(new_dv,n_elements(newz),n_elements(ell))
  restore,'sav_data/dv.sav'
  dv=rebin(new_dv,n_elements(newz),n_elements(ell))
  
  ;Integrate z using finer grid
  zstep=diff(newz)
  zstep=rebin(zstep,n_elements(newz)-1,n_elements(ell))
  z_integrand=new_dv*new_m_int
  c_l=total(z_integrand[1:n_elements(newz)-1,*]*zstep[*,*],1,/double)
  
  ;Integrate z using course grid
  zstep2=diff(z)
  zstep2=rebin(zstep2,n_elements(z)-1,n_elements(ell))
  z_integrand2=dv*m_int
  c_l2=total(z_integrand2[1:n_elements(z)-1,*]*zstep2[*,*],1,/double)

  ;Integrate z using int_tabulated
  test=int_tabulated(newz,z_integrand[*,ell_test],/double)
  test2=int_tabulated(z,z_integrand2[*,ell_test],/double)

  ;Comare methods
  print,'z integral:'
  print,'Interpolated int_tab: ',test
  print,'Interpolated simple: ',c_l[ell_test]
  print,'Original int_tab: ',test2
  print,'Original simple: ',c_l2[ell_test]

  stop
END
