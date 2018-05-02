PRO go

  z=1.
  logm=11.

  r_real=dindgen(30000000)/1000000.+0.0001
  step_real=diff(r_real)
  rho_real=nfw(r_real,logm,z)
  int_real=int_tabulated(r_real,rho_real,/double)
;  int_est=int_simp(r_real,rho_real)
  int_est=total(rho_real[1:n_elements(r_real)-1]*step_real)

  print,'Full int tabulated, int summed: ',int_real,int_est
  print,int_simp(r_real,rho_real)
  
  rlin=findgen(3000)/100.+0.0001
  rlog=cgloggen(3000,start=0.0001,finish=30,/double)
  
  step_lin=diff(rlin)
  step_log=diff(rlog)

  nfw_lin=nfw(rlin,logm,z)
  nfw_log=nfw(rlog,logm,z)

  plot,r_real,rho_real,linestyle=0,/xlog,/ylog
  oplot,rlin,nfw_lin,linestyle=2,color=cgcolor('red')
  oplot,rlog,nfw_log,linestyle=3,color=cgcolor('blue')

  int_lin=int_tabulated(rlin,nfw_lin)
  int_log=int_tabulated(rlog,nfw_log)

  int_lin_tot=total(nfw_lin[1:n_elements(rlin)-1]*step_lin)
  int_log_tot=total(nfw_log[1:n_elements(rlog)-1]*step_log)
  
  print,'Low res int tabulated linear, log: ', int_lin,int_log
  print,'Low res int summed linear, log: ',int_lin_tot,int_log_tot


  stop
  newr=dindgen(30000000)/1000000.+0.0001
  newstep=diff(newr)

  rho2=spline(r1,nfw1,newr)
  oplot,newr,rho2,linestyle=1,color=cgcolor('yellow')
  
  stop
END
