PRO go

  COMMON cosmological_parameters

  z=[10.]
  logm=[6,10,13.,16]
;  logm=[11,12,13,14,15.]
  
  phi_y=y_ell(logm,z,ell=ell)


  
  plot,ell,phi_y[0,0,*],linestyle=0,/xlog,yra=[2.1632d-26,2.16327d-26],ysty=1,xtit='l',ytit='y_l',charsize=2
  plot,ell,phi_y[2,0,*],linestyle=0,/xlog,ysty=1,xtit='l',ytit='y_l',charsize=2
  

;  stop
  
  plot,ell,phi_y[0,0,*]/max(phi_y[0,0,*]),linestyle=0,/xlog,ysty=1,xtit='l',ytit='y_l',charsize=2,/ylog
  oplot,ell,phi_y[1,0,*]/max(phi_y[1,0,*]),linestyle=1
    oplot,ell,phi_y[2,0,*]/max(phi_y[2,0,*]),linestyle=2
  
  
  stop
END
