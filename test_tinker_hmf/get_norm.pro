PRO go

  COMMON cosmological_parameters
  
  z=cgloggen(100,start=0.005,finish=9.995)
  logm=findgen(108)/10.+5
  newm=findgen(1071)/100.+5

  omega_z=evolve_density(z,type='matter')   
  Delta_cr=(18.*!dpi^2) + (82*(omega_z-1.)) - (39.*(omega_z-1)^2.)
  delta=Delta_cr/omega_z

  pspec=file_search('../power_spec/camb_matterpower*')

  norm=fltarr(n_elements(z))
  FOR i=0L,n_elements(z)-1 DO BEGIN
     counter,i,n_elements(z)
     out=mhalo2bias(newm,z[i],power_spec=pspec[i],delta=delta[i],/silent)
     nu=out[*,0]
     b=out[*,1]
     out=halo_mass_function(newm,z[i],power_spec=pspec[i],delta=delta[i],/silent)
     fnu=out[*,1]
     norm[i]=int_tabulated(alog10(nu),b*fnu*alog(10)*nu,/double)
  ENDFOR


  plot,z,1./norm,line=0

  openw,lun,'dndm_norm.txt',/get_lun
  FOR i=0L,n_elements(z)-1 DO printf,lun,strtrim(z[i],2) + '  ' + strtrim(norm[i],2),format='(A)'
  free_lun,lun
  
  stop
END
