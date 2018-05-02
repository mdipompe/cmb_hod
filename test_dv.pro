PRO go

  COMMON cosmological_parameters

  h=1
  omega_m=0.2
  omega_l=0.8
  
  z=cgloggen(100,start=0.005,finish=9.995)
  d=cosmocalc(z)
  v=d.v_c
  dv1=diff(d.v_c)/diff(z)
  
  E_z=sqrt((omega_m*(1.+z)^3.) + omega_l)

  dv2=(d.d_h*h)*((1.+z)^2.)*((d.d_a*h)^2.)/E_z

  plot,z,dv1,linestyle=0
  oplot,z,dv2*4*!dpi,linestyle=2
  
;  plot,z,dv2/(d.d_h^3.),linestyle=2
  
  stop
END
