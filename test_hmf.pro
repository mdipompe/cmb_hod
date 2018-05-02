PRO go

  COMMON cosmological_parameters
  
  z=cgloggen(100,start=0.005,finish=9.995)
;  logm=findgen(100)/10.+5
;  logm=findgen(991)/100.+5
;  logm=findgen(9901)/1000.+5
  logm=findgen(99001)/10000.+5.
  
  pspec=file_search('power_spec/camb_matterpower*')

  hmf1=halo_mass_function(logm,z[70],power_spec=pspec[70],model='tinker08')
  hmf2=halo_mass_function(logm,z[70],power_spec=pspec[70],model='tinker08',/logm)
;  plot,10^logm,hmf1,linestyle=0,/xlog,/ylog,charsize=2,thick=2
;  oplot,10^logm,hmf2/((10.^logm)*alog(10)),linestyle=2,color=cgcolor('red'),thick=2

  
  test1=int_tabulated(10.^logm,hmf1,/double)
  test2=int_tabulated(logm,hmf2,/double)
  print,'Int tabulated linear: ',test1
  print,'Int tabulated log: ',test2
  step1=diff(10.^logm)
  test3=total(hmf1*step1,/double)
  step2=diff(logm)
  test4=total(hmf2*step2,/double)
  print,'Total linear: ',test3
  print,'Total log: ',test4
  
  stop
    
  logm2=findgen(991)/100.+5
  hmf2=spline(logm,hmf1,logm2)

  plot,logm,hmf1,/ylog
  oplot,logm2,hmf2,linestyle=2,color=cgcolor('red')
  
  print,int_tabulated(10.^logm2,hmf2,/double)

  step=diff(10.^logm)
  print,total(hmf1*step,/double)
  step2=diff(10.^logm2)
  print,total(hmf2*step2,/double)
  
  stop
END
