PRO go

  COMMON cosmological_parameters
  
  ;Matching battaglia et al. cosmology
  omega_b=0.043
  omega_m=0.25
  omega_l=0.75
  sigma8=0.8
  h=0.7
  
  r=findgen(100000)/10000.+0.0001
  
  logm=[11,12,13,14,15.]
  z=[0.,1,2.]

  p=sz_pressure(r,logm,z,m_200=m_200,r_200=r_200,p_200=p_200)
  r_s=r_200


  loadct,33
  
  PS_start,filename='sz_test.png'
  
  plot,[0],[0],/nodata,xra=[1d-4,10.],yra=[1d-35,1d-15],xsty=1,ysty=1,/xlog,/ylog,$
       xtit=textoidl('r [Mpc/h]'),ytit=textoidl('P_e [M_{sun}/h / Mpc/h s^2]'),$
       charsize=2,color=cgcolor('black')

  
  FOR i=0L,n_elements(logm)-1 DO BEGIN
     FOR j=0L,n_elements(z)-1 DO Begin
        oplot,r,p[i,j,*],linestyle=j,color=i*50,thick=2
     ENDFOR
  ENDFOR
  
  legend,['z=0','z=1','z=2'],linestyle=[0,1,2],$
         box=0,/bottom,/left,thick=2,charsize=2
  circsym,/fill
  legend,['logM=11','logm=15'],psym=[8,8],color=[0,i*50],$
         box=0,/top,/right

  PS_end,/png

  j=0
  ptot=fltarr(n_elements(logm))
  FOR i=0L,n_elements(logm)-1 DO BEGIN
;     xx=where(r LE r_200[i,j])
;     ptot[i]=int_tabulated(r[xx],p[i,j,xx],/double)

;     xx=closest(r,r_200[i,j])
;     ptot[i]=p[i,j,xx]

     ptot[i]=p[i,j,1000]
     
  ENDFOR

  plot,logm,alog10(ptot),psym=1
  fit=linfit(logm,alog10(ptot))
  x=findgen(100)
  y=fit[0]+fit[1]*x
  oplot,x,y,linestyle=2
  


  
  stop
  
END
