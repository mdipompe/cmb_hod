PRO go

  COMMON cosmological_parameters
  
  z=cgloggen(100,start=0.005,finish=9.995)
  logm=findgen(108)/10.+5
  newm=findgen(1071)/100.+5

  rho_c=rho_crit(z,/physical)
  om_z=evolve_density(z,type='matter')
  rho_m=rho_c*om_z
  
  restore,'bias_interp.sav'
  restore,'dndm_interp.sav'
  
  int=fltarr(n_elements(z))
  int2=fltarr(n_elements(z))
  FOR i=0L,n_elements(z)-1 DO BEGIN
     integrand=new_dndm[*,i]*b[*,i]*((10.^newm)/rho_m[i])
     int[i]=int_tabulated(10.^newm,integrand,/double)

     step=diff(10.^newm)
     int2[i]=total(step*integrand)
  ENDFOR
  
  plot,z,int,linestyle=2
  stop
  plot,z,int2,linestyle=2


  stop

  restore,'dndlogm_interp.sav'
  
  int=fltarr(n_elements(z))
  int2=fltarr(n_elements(z))
  FOR i=0L,n_elements(z)-1 DO BEGIN
     integrand=new_dndlogm[*,i]*b[*,i]*((newm)/rho_m[i])
     int[i]=int_tabulated(newm,integrand,/double)

     step=diff(newm)
     int2[i]=total(step*integrand)
  ENDFOR

  plot,z,int,linestyle=2
  stop
  plot,z,int2,linestyle=2

  
  stop
END
