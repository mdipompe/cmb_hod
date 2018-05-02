PRO go

  logm=[11,13,15.]
  z=[0,1,3,5.]

  c_zheng=c_index(logm,z)
  c_duffy=c_index(logm,z,/duffy)

  plot,logm,c_zheng[*,0],psym=1,xra=[10,16]
  oplot,logm,c_duffy[*,0],psym=6

  stop

  plot,z,c_zheng[1,*],psym=1,xra=[-0.5,6]
  oplot,z,c_duffy[1,*],psym=6

  r=findgen(100000)/10000.+0.0001
  rho_zheng=nfw(r,logm,z,rvir=rvir_z,cindex=c_zheng)
  rho_duffy=nfw(r,logm,z,rvir=rvir_d,cindex=c_duffy)

  plot,r,rho_zheng[1,2,*],linestyle=0,/xlog,/ylog
  oplot,r,rho_duffy[1,2,*],linestyle=2

  
  stop
END

  
