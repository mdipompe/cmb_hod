PRO go

  COMMON cosmological_parameters
  
  logm=[11,12,13.]
  z=[0,1.,2]

  r_all=r_vir(logm,z)


  r00=r_vir(logm[0],z[0])
  r01=r_vir(logm[0],z[1])
  r02=r_vir(logm[0],z[2])

  r10=r_vir(logm[1],z[0])
  r11=r_vir(logm[1],z[1])
  r12=r_vir(logm[1],z[2])
  
  r20=r_vir(logm[2],z[0])
  r21=r_vir(logm[2],z[1])
  r22=r_vir(logm[2],z[2])
  

  print,r_all[0,0],r00
  print,r_all[0,1],r01
  print,r_all[0,2],r02
  print,r_all[1,0],r10
  print,r_all[1,1],r11
  print,r_all[1,2],r12
  print,r_all[2,0],r20
  print,r_all[2,1],r21
  print,r_all[2,2],r22
  
  print,' '

  r_m=r_vir(logm,z[0])
  print,r_m[0],r00
  print,r_m[1],r10
  print,r_m[2],r20

  print,' '
  
  r_z=r_vir(logm[0],z)
  print,r_z[0],r00
  print,r_z[1],r01
  print,r_z[2],r02
  
  
  
  
  
  
  stop
END
