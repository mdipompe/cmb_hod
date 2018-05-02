PRO go

  COMMON cosmological_parameters

  omega_m=0.25
  omega_l=0.75
  sigma8=1.3
  h=0.75
  
  r=findgen(1000000)/10.+1
  r_mpc=(r/1000.)

  z=0.
  
  m=[4.1d13*0.01,4.1d13*10]
  logm=alog10(m)

  rho_crit=rho_crit(z)
  
;  cindex=c_index(logm,z,mstar=mstar,/get_mstar)
;  stop
  
  r_vir=r_vir(logm,z,delta_cr=200)

  ;These are based off of Figure 6 in NFW
  cindex=[11.,5.]
  
  rho=nfw(r_mpc,logm,z,rvir=r_vir,delta_cr=200,cindex=cindex)
  rho_noh=rho*h*(1./(1000.*h))^3.
  
  plot,alog10(r*h),alog10(rho_noh[0,*]/1d10),linestyle=0,$
       xra=[0.1,4.1],yra=[-9.,-2],xsty=1,ysty=1,$
       xtit='Log(R [kpc])',ytit=textoidl('Log(\rho/10^{10} [M_{sun}/kpc^3])'),$
       charsize=2
  oplot,alog10(r*h),alog10(rho_noh[1,*]/1d10),linestyle=2


  xx=where(r_mpc LE r_vir[0])
  m_test=int_tabulated(r_mpc[xx],4*!dpi*r_mpc[xx]^2.*rho[0,xx],/double)
  print,m_test,m[0]

  xx=where(r_mpc LE r_vir[1])
  m_test=int_tabulated(r_mpc[xx],4*!dpi*r_mpc[xx]^2.*rho[1,xx],/double)
  print,m_test,m[1]

  
  
  
  stop
END
