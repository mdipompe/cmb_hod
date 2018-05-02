PRO go

  COMMON cosmological_parameters

  logm=[5.,12,16.]
  z=[0.5,2,9.]

  test1=convert_nfw(logm,z,r_deltacr=newr)
  test2=convert_nfw(logm,z,r_deltacr=newr2,/wrt_m)

  print,test1/test2
  print,newr/newr2


  stop
  
  z=1.
  logm=findgen(50)/10+11.
  rvir=r_vir(logm,z)
  c=c_index(logm,z)
  
  r2=fltarr(n_elements(logm))
  logm2=fltarr(n_elements(logm))
  c2=fltarr(n_elements(logm))

  r3=r2
  logm3=logm2
  c3=c2

  r=findgen(100000)/10000.+0.0001
  
  FOR i=0L,n_elements(logm)-1 DO BEGIN
     counter,i,n_elements(logm)
     logm2[i]=alog10(convert_nfw(logm[i],z,delta_cr=200.,r_deltacr=r_tmp,c_deltacr=c_tmp))
     r2[i]=r_tmp
     c2[i]=c_tmp
     test=convert_nfw_slow(r,logm[i],z,delta_cr=200.)
     r3[i]=test[0]
     logm3[i]=test[1]
     c3[i]=test[2]
  ENDFOR


  stop
END
