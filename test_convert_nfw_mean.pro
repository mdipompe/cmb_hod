PRO go

  COMMON cosmological_parameters

  z=3.
;  logm=findgen(50)/10+11.
  logm=findgen(12)+5
  rvir=r_vir(logm,z)
  c=c_index(logm,z)
  
  r2=fltarr(n_elements(logm))
  logm2=fltarr(n_elements(logm))
  c2=fltarr(n_elements(logm))

  FOR i=0L,n_elements(logm)-1 DO BEGIN
     counter,i,n_elements(logm)
     test=convert_nfw_mean(logm[i],z,delta_cr=200.)
     r2[i]=test[0]
     logm2[i]=test[1]
  ENDFOR


  stop
END
