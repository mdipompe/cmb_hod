PRO go

  COMMON cosmological_parameters

  r=findgen(100000)/10.+1
  z=1.
  mstar=3.71535d+12
  r_mpc=(r/1000.)
  

  ;Try for multiple M and z at once
  z=[1.]
  logm=[11,13,15.]
  r=r_mpc
  
  r_vir=r_vir(logm,z)

  rho=nfw(r,logm,z,r_vir=r_vir,mstar=mstar,/project,type='num',c0=32)


  

  r=findgen(100000)/10.+1
  z=1.
  mstar=3.71535d+12
  
  r_mpc=(r/1000.)

  r_vir15=r_vir(15.,z)
  r_vir11=r_vir(11.,z)
  
  rho_m15=nfw(r_mpc,15.,z,r_vir=r_vir15,mstar=mstar,c0=32,type='num',/project)*h*(1./(1000.*h))^3.
  rho_m11=nfw(r_mpc,11.,z,r_vir=r_vir11,mstar=mstar,c0=32,type='num',/project)*h*(1./(1000.*h))^3.

  rho_m15_mpc=nfw(r_mpc,15.,z,r_vir=r_vir15,mstar=mstar,c0=32,type='num',/project)
  rho_m11_mpc=nfw(r_mpc,11.,z,r_vir=r_vir11,mstar=mstar,c0=32,type='num',/project)
  stop
  use15=where(r_mpc LE r_vir15*2./!dpi)
  N15=int_tabulated(r_mpc[use15],2*!dpi*r_mpc[use15]*rho_m15_mpc,/double)
  use11=where(r_mpc LE r_vir11*2./!dpi)
  N11=int_tabulated(r_mpc[use11],2*!dpi*r_mpc[use11]*rho_m11_mpc,/double)
  print,N15,N11

  
  
  ;plot,r,rho_m15,linestyle=2,/xlog,/ylog,xra=[1d-1,10000],$
  ;     xtit='Log(R [Mpc/h])',ytit=textoidl('\rho [M_{sun}/Mpc^3 h^2]'),$
  ;     charsize=2
  ;oplot,r,rho_m12,linestyle=1
  plot,alog10(r),alog10(rho_m15),linestyle=2,xra=[0.1,4],yra=[-14,-5],$
       xtit='Log(R [kpc])',ytit=textoidl('Log(\Sigma [1/kpc^2])'),$
       charsize=2
  oplot,alog10(r),alog10(rho_m11),linestyle=1
  
  legend,['LogM=15','LogM=11'],linestyle=[2,1],box=0,$
         charsize=1.5,charthick=1.5,/top,/right



  
  stop
END
