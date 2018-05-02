PRO get_hmf

  ;Get the halo mass function for your cosmology, z, and M values
  
  COMMON cosmological_parameters

  z=cgloggen(100,start=0.005,finish=9.995)
  logm=findgen(108)/10.+5
  newm=findgen(1071)/100.+5
  pspec=file_search('power_spec/camb_matterpower*')

  ;This is a hack to properly re-normalize the Tinker HMF
  readcol,'dndm_norm.txt',nz,norm

  ;Since HMF from Tinker 2010 defined relative to rho_mean,
  ;need to figure out what Delta gives same mass as one I'm using
  ;everywhere else. This is just the ratio between my Delta and omega
  ;matter.
  omega_z=evolve_density(z,type='matter')   
  Delta_cr=(18.*!dpi^2) + (82*(omega_z-1.)) - (39.*(omega_z-1)^2.)
  delta=Delta_cr/omega_z

;  ;Since HMF from Tinker 2010 defined relative to rho_mean,
;  ;convert to equivalent mass
;  restore,'sav_data/Mvir_Mmean200.sav'
  
  print,'Working on dndm, normal grid:'  
  dndm=findgen(n_elements(logm),n_elements(z))
  FOR i=0L,n_elements(z)-1 DO BEGIN
     counter,i,n_elements(z)
;     dndm[*,i]=halo_mass_function(logm[*,i], z[i],power_spec=pspec[i],/silent)
     dndm[*,i]=halo_mass_function(logm, z[i],power_spec=pspec[i],delta=delta[i],/silent)/norm[i]
;     dndm[*,i]=halo_mass_function(logm_mean[*,i], z[i],power_spec=pspec[i],/silent)
  ENDFOR
;  save,dndm,filename='sav_data/dndm.sav'
  save,dndm,filename='sav_data/dndm_renorm.sav'
;  save,dndm,filename='sav_data/dndm_converted_m.sav


  print,'Working on dndm, fine grid:'
  new_dndm=fltarr(n_elements(newm),n_elements(z))
  FOR i=0L,n_elements(z)-1 DO BEGIN
     counter,i,n_elements(z)
     new_dndm[*,i]=spline(logm,dndm[*,i],newm)
  ENDFOR
;  save,new_dndm,filename='sav_data/dndm_interp.sav'
  save,new_dndm,filename='sav_data/dndm_renorm_interp.sav'
;  save,new_dndm,filename='sav_data/dndm_converted_m_interp.sav'
  

  stop

  ;Doing some testing in log space, to see if that makes a difference in
  ;integral accuracy
  print,'Working on dndlogm, normal grid:'
  dndlogm=findgen(n_elements(logm),n_elements(z))
  FOR i=0L,n_elements(z)-1 DO BEGIN
     counter,i,n_elements(logm)
     dndlogm[*,i]=halo_mass_function(logm, z[i],power_spec=pspec[i],/logm,delta=delta[i],/silent)
  ENDFOR
  save,dndlogm,filename='sav_data/dndlogm.sav'

  print,'Working on dndlogm, fine grid:'
  new_dndlogm=fltarr(n_elements(newm),n_elements(z))
  FOR i=0L,n_elements(z)-1 DO BEGIN
     counter,i,n_elements(z)
     new_dndlogm[*,i]=spline(logm,dndlogm[*,i],newm)
  ENDFOR
  save,new_dndlogm,filename='sav_data/dndlogm_interp.sav'
  
  stop
END
