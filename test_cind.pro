PRO go

  COMMON cosmological_parameters
  
  logm=[11,12,13.]
  z=[0,1,2.]

  c_all=c_index(logm,z)


  c00=c_index(logm[0],z[0],mstar=mstar)
  c01=c_index(logm[0],z[1],mstar=mstar)
  c02=c_index(logm[0],z[2],mstar=mstar)

  c10=c_index(logm[1],z[0],mstar=mstar)
  c20=c_index(logm[2],z[0],mstar=mstar)

  c11=c_index(logm[1],z[1],mstar=mstar)
  c22=c_index(logm[2],z[2],mstar=mstar)


  print,c_all[0,0],c00
  print,c_all[0,1],c01
  print,c_all[0,2],c02
  print,c_all[1,0],c10
  print,c_all[2,0],c20
  print,c_all[1,1],c11
  print,c_all[2,2],c22

  print,' '

  c_m=c_index(logm,z[0])
  print,c_m[0],c00
  print,c_m[1],c10
  print,c_m[2],c20

  print,' '

  c_z=c_index(logm[0],z)
  print,c_z[0],c00
  print,c_z[1],c01
  print,c_z[2],c02

  stop
END
