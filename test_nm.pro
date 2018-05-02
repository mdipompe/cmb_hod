PRO go

  logm=findgen(60)/10.+11
  
  M_cut=10.^5    
  ;MAD Set other default HOD parameters (Richardson et al. 2013)
  IF (n_elements(M_min) EQ 0) THEN M_min=10.^13.65
  logM_min=alog10(M_min)
  IF (n_elements(sig_logM) EQ 0) THEN sig_logM=0.78
  IF (n_elements(M1) EQ 0) THEN M1=10.^14.32
  IF (n_elements(alpha) EQ 0) THEN alpha=2.59

  ;MAD Calculate mean N(M) from HOD
  N_m_cent=0.5*(1+erf((logm-logM_min)/sig_logM))
  N_m_sat=(((10.^logm)/M1)^alpha)*exp((-1.)*M_cut/(10.^logm))
  N_m=N_m_cent+N_m_sat

  plot,10.^logm,n_m_cent,linestyle=2,/ylog,/xlog,$
       xra=[10.^12,10.^15],xsty=1,$
       yra=[10.^(-3),10.],ysty=1,$
       xtit='M [M_sun/h]',ytit='<N(M)>',charsize=2
  oplot,10.^logm,n_m_sat,linestyle=1
  oplot,10.^logm,n_m,linestyle=3

  legend,['Central','Satellite','Total'],linestyle=[2,1,0],$
         box=0,charsize=1.5,charthick=1.5,/top,/left
  
  print,n_m[where(logm EQ 15.)]
  print,n_m[where(logm EQ 11.)]
  
  stop
END
