PRO one_halo

  COMMON cosmological_parameters

  ell=cgloggen(1000,start=1,finish=10000)
  z=cgloggen(100,start=0.005,finish=9.995)
  logm=findgen(108)/10.+5
  newm=findgen(1071)/100.+5
  newz=cgloggen(10000,start=0.005,finish=9.995)
  
  restore,'../sav_data/dv_interp.sav'
  new_dv=rebin(new_dv,n_elements(newz),n_elements(ell))
  
  restore,'../sav_data/dndm_renorm_interp.sav'
  new_dndm=rebin(new_dndm,n_elements(newm),n_elements(z),n_elements(ell))

  restore,'../sav_data/phi_y_mz_interp.sav'
  undefine,new_y

  restore,'../sav_data/phi_gal_mz_interp.sav'

  ;Get dNdz
  readcol,'all_dndz.txt',wise_z
  fit_dndz,wise_z,newz,dndz,binsize=0.15
  dndz=rebin(dndz,n_elements(newz),n_elements(ell))

  ;Need a factor so that you don't get a bunch of 0s for very small 
  ;numbers when multiplying for integral. Factored back out later.
  f=1d4
  new_phi_gal=new_phi_gal*f

  ;Do the mass integral
  mstep=diff(10.^newm)
  mstep=rebin(mstep,n_elements(newm)-1,n_elements(z),n_elements(ell))
  m_integrand=new_dndm*new_phi*new_phi_gal
  m_int=total(m_integrand[1:n_elements(newm)-1,*,*]*mstep[*,*,*],1,/double)/(f^2.)

  new_m_int=findgen(n_elements(newz),n_elements(ell))
  FOR i=0L,n_elements(ell)-1 DO BEGIN
     counter,i,n_elements(ell)
     new_m_int[*,i]=spline(z,m_int[*,i],newz)     
  ENDFOR

  ;Do the z integral
  zstep=diff(newz)
  zstep=rebin(zstep,n_elements(newz)-1,n_elements(ell))
  z_integrand=dndz*new_dv*new_m_int
  c_l_1h=total(z_integrand[1:n_elements(newz)-1,*]*zstep[*,*],1,/double)

  save,c_l_1h,filename='c_l_1h.sav'
  
  xtit='l'
  ytit=textoidl('l^2(l+1)C_l^{y\phi}/2\pi')
  plot,ell,((ell^2)*(ell+1))*c_l_1h/(2.*!dpi),/xlog,$
       xtit=xtit,ytit=ytit,charsize=2,xra=[10,10000],xsty=1,thick=2,ysty=1

  
  stop
END
