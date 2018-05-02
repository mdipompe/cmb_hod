PRO go

  COMMON cosmological_parameters
  
  ;Matching battaglia et al. cosmology
  omega_b=0.043
  omega_m=0.25
  omega_l=0.75
  sigma8=0.8
  h=0.7
  
  r=findgen(100000)/10000.+0.0001
  
  logm=[13.85,15.5]
  z=[0.]

  ;The factor of 1.932 puts this back as the thermal pressue,
  ;which is what is plotted in Battaglia et al., while the code
  ;returns the electron pressure P_e.
  p=sz_pressure(r,logm,z,m_200=m_200,r_200=r_200,p_200=p_200)*1.932
  r_s=r_200

  r_vir=r_vir(logm,z)
  print,r_vir/r_s

  
  plot,[0],[0],/nodata,xra=[0.03,2.5],yra=[0.0008,0.22],xsty=1,ysty=1,/xlog,/ylog,$
       xtit=textoidl('r/R_{200}'),ytit=textoidl('P/P_{200} * (r/R_{200})^3'),$
       charsize=2

  readcol,'sz_data_lowM.txt',r_low,p_low,format='D',delim=','
  readcol,'sz_data_highM.txt',r_high,p_high,format='D',delim=','

  oplot,r_low,p_low,linestyle=2,color=cgcolor('red'),thick=3
  oplot,r_high,p_high,linestyle=2,color=cgcolor('magenta'),thick=3

  oplot,r/r_s[0],(p[0,0,*]/p_200[0])*((r/r_s[0])^3.),linestyle=0,color=cgcolor('red'),thick=2
  oplot,r/r_s[1],(p[1,0,*]/p_200[1])*((r/r_s[1])^3.),linestyle=0,color=cgcolor('magenta'),thick=2

  circsym,/fill
  legend,[textoidl('M_{200}=5.6 \times 10^{13}'),textoidl('M_{200}=2.4 \times 10^{15}')],$
         psym=[8,8],color=[cgcolor('red'),cgcolor('magenta')],box=0,/bottom,/right,$
         thick=2,charsize=2
  legend,['Battaglia et al.','Mine'],$
         linestyle=[2,0],box=0,$
         thick=2,charsize=2,$
         position=[0.3,0.005],/data
  


  stop

  undefine,m_200,r_200,p_200,r_s,r_vir



  logm=[14.25]
  z=[0.,1]

  ;The factor of 1.932 puts this back as the thermal pressue,
  ;which is what is plotted in Battaglia et al., while the code
  ;returns the electron pressure P_e.
  p=sz_pressure(r,logm,z,m_200=m_200,r_200=r_200,p_200=p_200)*1.932
  r_s=r_200

  r_vir=r_vir(logm,z)
  print,r_vir/r_s
  
  plot,[0],[0],/nodata,xra=[0.03,2.5],yra=[0.0008,0.22],xsty=1,ysty=1,/xlog,/ylog,$
       xtit=textoidl('r/R_{200}'),ytit=textoidl('P/P_{200} * (r/R_{200})^3'),$
       charsize=2

  readcol,'sz_data_lowz.txt',r_low,p_low,format='D',delim=','
  readcol,'sz_data_highz.txt',r_high,p_high,format='D',delim=','

  oplot,r_low,p_low,linestyle=2,color=cgcolor('red'),thick=3
  oplot,r_high,p_high,linestyle=2,color=cgcolor('magenta'),thick=3

  oplot,r/r_s[0,0],(p[0,0,*]/p_200[0,0])*((r/r_s[0,0])^3.),linestyle=0,color=cgcolor('red'),thick=2
  oplot,r/r_s[0,1],(p[0,1,*]/p_200[0,1])*((r/r_s[0,1])^3.),linestyle=0,color=cgcolor('magenta'),thick=2

  circsym,/fill
  legend,[textoidl('z=0'),textoidl('z=1')],$
         psym=[8,8],color=[cgcolor('red'),cgcolor('magenta')],box=0,/bottom,/right,$
         thick=2,charsize=2
  legend,['Battaglia et al.','Mine'],$
         linestyle=[2,0],box=0,$
         thick=2,charsize=2,$
         position=[0.3,0.005],/data
  

  
END
