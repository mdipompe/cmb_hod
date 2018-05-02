PRO get_b

  ;Get the bias for the 2-halo calculations
  
  COMMON cosmological_parameters

  z=cgloggen(100,start=0.005,finish=9.995)
  logm=findgen(108)/10.+5
  newm=findgen(1071)/100.+5
  pspec=file_search('power_spec/camb_matterpower*')

  ;Since HMF from Tinker 2010 defined relative to rho_mean,
  ;need to figure out what Delta gives same mass as one I'm using
  ;everywhere else. This is just the ratio between my Delta and omega
  ;matter.
  omega_z=evolve_density(z,type='matter')   
  Delta_cr=(18.*!dpi^2) + (82*(omega_z-1.)) - (39.*(omega_z-1)^2.)
  delta=Delta_cr/omega_z

;  restore,'sav_data/Mvir_Mmean200.sav'
  
  b=fltarr(n_elements(logm),n_elements(z))
  FOR i=0L,n_elements(z)-1 DO BEGIN
     counter,i,n_elements(z)
     b[*,i]=mhalo2bias(logm,z[i],power_spec=pspec[i],delta=delta[i],/silent)
;     b[*,i]=mhalo2bias(logm_mean[*,i],z[i],power_spec=pspec[i],/silent)
  ENDFOR
  save,b,filename='sav_data/bias.sav'
;  save,b,filename='sav_data/bias_converted_m.sav'


  ;Interpolate in mass for more accuate integrals
  b=fltarr(n_elements(newm),n_elements(z))
  FOR i=0L,n_elements(z)-1 DO BEGIN
     counter,i,n_elements(z)
;     quadterp,logm,logm_mean[*,i],newm,newm_mean
     b[*,i]=mhalo2bias(newm,z[i],power_spec=pspec[i],delta=delta[i],/silent)
;     b[*,i]=mhalo2bias(newm_mean,z[i],power_spec=pspec[i],/silent)
  ENDFOR
  save,b,filename='sav_data/bias_interp.sav'
;  save,b,filename='sav_data/bias_converted_m_interp.sav'

  stop
END
