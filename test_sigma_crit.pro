PRO go

  COMMON cosmological_parameters
  h=0.704
  sigma8=0.81
  omega_m=0.272
  omega_b=0.0455
  omega_l=1-omega_m
  slope=0.967
  
  zcmb=findgen(16)/10.+0.02

  z=[1,2,3,5,8]/10.

  s=findgen(n_elements(zcmb),n_elements(z))

  FOR i=0L,n_elements(zcmb)-1 DO BEGIN
     s[i,*]=sigma_crit_lens(z,cmb_z=zcmb[i])/((1d6)^2)/(1.+z)
  ENDFOR

  xx=where(s LT 0)
  s[xx]=20000
  
  plot,[0],[0],/nodata,xra=[0.1,1.5],yra=[1000,10000],xsty=1,ysty=1,$
       xtit=textoidl('z_{s}'),ytit=textoidl('\Sigma_{crit} [h M_{sun} pc^{-2}]'),$
       charsize=2

  ls=[0,1,2,3,4]
  FOR i=0L,n_elements(z)-1 DO BEGIN
     oplot,zcmb,s[*,i],linestyle=ls[i]
  ENDFOR



  stop


  zcmb=[1.,0.8,1.2]
  z=findgen(30)/100+0.01

  s=findgen(n_elements(zcmb),n_elements(z))
  FOR i=0L,n_elements(zcmb)-1 DO BEGIN
     s[i,*]=1./(sigma_crit_lens(z,cmb_z=zcmb[i])/(1d15))
  ENDFOR

  plot,[0],[0],/nodata,xra=[0.01,0.3],yra=[0,0.3],xsty=1,ysty=1,$
       xtit=textoidl('z_{L}'),ytit=textoidl('\Sigma_{crit}^{-1} [Mpc/h 10^{15}M_{sun}]'),$
       charsize=2

  ls=[0,2,3]
  FOR i=0L,n_elements(zcmb)-1 DO BEGIN
     oplot,z,s[i,*],linestyle=ls[i]
  ENDFOR

  
  
  stop
END
