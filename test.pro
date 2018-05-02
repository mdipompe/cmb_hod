PRO go
  
  logm=[11.]
  z=[1.]
  rvir=r_vir(logm,z)
  
  newm=convert_nfw(logm,z,r_deltacr=r_deltacr,/rho_m)

  
  stop
END
