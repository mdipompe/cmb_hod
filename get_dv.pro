PRO get_dv

  ;Get the differential comoving volume element
  
  COMMON cosmological_parameters

  z=cgloggen(100,start=0.005,finish=9.995)
  newz=cgloggen(10000,start=0.005,finish=9.995)

  ;Need to set h=1 for my cosmology calculator,
  ;so you get per h units
  orig_h=h
  h=1.
  d=cosmocalc(z)
  E_z=sqrt((omega_m*(1.+z)^3.) + omega_l)
  dv=(d.d_h*h)*((1.+z)^2.)*((d.d_a*h)^2.)/E_z
  ;Then set it back
  h=orig_h
  save,dv,filename='sav_data/dv.sav'

  ;And now do it for the finer grid of z
  orig_h=h
  h=1.
  d=cosmocalc(newz)
  E_z=sqrt((omega_m*(1.+newz)^3.) + omega_l)
  new_dv=(d.d_h*h)*((1.+newz)^2.)*((d.d_a*h)^2.)/E_z
  h=orig_h
  save,new_dv,filename='sav_data/dv_interp.sav'

  
  stop
END
