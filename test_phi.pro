PRO go

  COMMON cosmological_parameters

  z=[1.,2,3,10]
  logm=[6.,11.,13,16]

  phi_m_dm=phi_ell(logm,z,ell=ell,type='dm')

  zind=0
  
  plot,ell,phi_m_dm[0,zind,*]*ell*(ell+1.),linestyle=0,/xlog,ysty=1,xtit='l',ytit='phi_l',charsize=2
;  plot,ell,phi_m_dm[0,0,*],linestyle=0,/xlog,/ylog,ysty=1,xtit='l',ytit='phi_l',charsize=2

  stop
  
  test=phi_m_dm[0,zind,*]*ell*(ell+1)/2.
  test2=phi_m_dm[1,zind,*]*ell*(ell+1)/2.
  test3=phi_m_dm[2,zind,*]*ell*(ell+1)/2.
  test4=phi_m_dm[3,zind,*]*ell*(ell+1)/2.
;  plot,ell,test/max(test),linestyle=0,/xlog,yra=[7.5d-14,8.5d-14],ysty=1,xtit='l',ytit='l(l+1)*phi_l',charsize=2

  
  plot,ell,test/max(test),linestyle=0,/xlog,ysty=1,xtit='l',ytit='l(l+1)*phi_l',charsize=2
  oplot,ell,test2/max(test2),linestyle=2
  oplot,ell,test3/max(test3),linestyle=3
  oplot,ell,test4/max(test4),linestyle=4
  
  
  
  stop
END
