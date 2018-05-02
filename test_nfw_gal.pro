PRO go

  COMMON cosmological_parameters

  r=findgen(100000)/10000.+0.0001

  ;Try for multiple M and z at once
  z=[1,3.]
  logm=[11,15.]
  
  rvir=r_vir(logm,z)

  rho=nfw(r,logm,z,rvir=rvir,cindex=cindex,type='num',c0=32)

  xx=where(r LE rvir[1,1])
  n_test=int_tabulated(r[xx],4.*!dpi*r[xx]^2.*rho[1,1,*],/double)
;  xx=where(r LE rvir)
;  n_test=int_tabulated(r[xx],4.*!dpi*r[xx]^2.*rho,/double)

  n_hod=hod(logm)

  print,n_test,n_hod[1]
  

  stop
  
  plot,alog10(r),alog10(rho[1,0]),linestyle=2,xra=[0.001,4],$
       xtit='Log(R [Mpc])',ytit=textoidl('Log(n [1/Mpc^3])'),$
       charsize=2
  oplot,alog10(r),alog10(rho[0,0]),linestyle=1
  
  legend,['LogM=15','LogM=11'],linestyle=[2,1],box=0,$
         charsize=1.5,charthick=1.5,/top,/right



  
  stop
END
