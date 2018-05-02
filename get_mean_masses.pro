PRO get_mean_masses

  ;Convert virial masses to masses based on 200*rho_mean for HMF calculation.

  z=cgloggen(100,start=0.005,finish=9.995)
  logm=findgen(108)/10.+5
  
  logm_mean=fltarr(n_elements(logm),n_elements(z))
  FOR i=0L,n_elements(z)-1 DO BEGIN
     counter,i,n_elements(z)
     logm_mean[*,i]=alog10(convert_nfw(logm,z[i],/wrt_m))
  ENDFOR 
  save,logm_mean,filename='sav_data/Mvir_Mmean200.sav'
  restore,'sav_data/Mvir_Mmean200.sav'
  

  stop
END
