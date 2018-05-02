PRO go

  COMMON cosmological_parameters
  
  r=findgen(100000)/10000.+0.0001

  ;Try for multiple M and z at once
;  z=[1,2,3,4,5,6.]
;  z=[0,0.5,1,1.5,2,2.5,3,3.5]
  z=[1]
  logm=[11.]

  r_vir=r_vir(logm,z)
  cindex=c_index(logm,z)
  
  rho=nfw(r,logm,z,rvir=r_vir,mstar=mstar,cindex=cindex)
  stop
  xx=where(r LE r_vir[0,0])
  mtest=int_tabulated(r[xx],4*!dpi*r[xx]^2.*rho[0,0,xx],/double)
  print,10.^(logm[0]),mtest
  
  xx=where(r LE r_vir[1,0])
  mtest=int_tabulated(r[xx],4*!dpi*r[xx]^2.*rho[1,0,xx],/double)
  print,10.^(logm[1]),mtest

  xx=where(r LE r_vir[0,1])
  mtest=int_tabulated(r[xx],4*!dpi*r[xx]^2.*rho[0,1,xx],/double)
  print,10.^(logm[0]),mtest
  
  xx=where(r LE r_vir[1,1])
  mtest=int_tabulated(r[xx],4*!dpi*r[xx]^2.*rho[1,1,xx],/double)
  print,10.^(logm[1]),mtest


  plot,r,rho[0,0,*],linestyle=0,/xlog,/ylog,$
       xtit='r [Mpc/h]',ytit=textoidl('\rho [(M_{sun}/h)/(Mpc/h)^3]'),$
       charsize=2
  oplot,r,rho[0,0,*],linestyle=0,color=cgcolor('red'),thick=2
  oplot,r,rho[1,0,*],linestyle=2,color=cgcolor('red'),thick=2

  oplot,r,rho[0,1,*],linestyle=0,color=cgcolor('magenta'),thick=2
  oplot,r,rho[1,1,*],linestyle=2,color=cgcolor('magenta'),thick=2
  
  legend,['Log(M)=11','Log(M)=15'],$
         linestyle=[0,2],box=0,/bottom,/left,$
         thick=2,charsize=2
  circsym,/fill
  legend,['z=1','z=6'],$
         psym=[8,8],color=[cgcolor('red'),cgcolor('magenta')],box=0,$
         thick=2,charsize=2,$
         position=[0.001,8d9],/data

  
  
  stop
END
