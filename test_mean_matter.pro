PRO go

  z=findgen(10)/10.

  COMMON cosmological_parameters
  
  omega_m_z=evolve_density(z,type='matter')

  rc=rho_crit(z)

  ;The extra factor of (1+z)^3 is because the critical density 
  ;is in physical units (or something like that...)
  rho_z=rc*omega_m_z/((1.+z)^3.)
  
  plot,z,rho_z,linestyle=2,xtit='z',ytit=textoidl('\rho_m(z) [M_sun/h / (Mpc/h)^3]'),charsize=2
  stop
  plot,z,rc/max(rc),linestyle=2,xtit='z',charsize=2
  oplot,z,omega_m_z/max(omega_m_z),linestyle=1
  legend,[textoidl('\Omega_m(z)/max(\Omega_m(z))'),textoidl('\rho_{crit}(z)/max(\rho_{crit}(z))')],linestyle=[1,2],box=0,$
         /bottom,/right,charsize=1.5

  stop
END
