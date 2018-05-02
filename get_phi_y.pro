PRO get_phi_y

  z=cgloggen(100,start=0.005,finish=9.995)
  logm=findgen(108)/10.+5

  phi_mz=fltarr(n_elements(logm),n_elements(z),1000)
  y_mz=fltarr(n_elements(logm),n_elements(z),1000)
  FOR i=0L,n_elements(z)-1 DO BEGIN
     counter,i,n_elements(z)
     phi_mz[*,i,*]=phi_ell(logm,z[i],ell=ell_phi,type='dm')
     y_mz[*,i,*]=y_ell(logm,z[i],ell=ell_y)
  ENDFOR
  save,phi_mz,y_mz,filename='sav_data/phi_y_mz.sav'
  restore,'sav_data/phi_y_mz.sav'

  newm=findgen(1071)/100.+5
  new_phi=fltarr(n_elements(newm),n_elements(z),n_elements(ell_phi))
  new_y=fltarr(n_elements(newm),n_elements(z),n_elements(ell_y))
  FOR i=0L,n_elements(z)-1 DO BEGIN
     counter,i,n_elements(z)
     FOR j=0L,n_elements(ell_phi)-1 DO BEGIN
        new_phi[*,i,j]=spline(logm,phi_mz[*,i,j],newm)
        new_y[*,i,j]=spline(logm,y_mz[*,i,j],newm)
     ENDFOR
  ENDFOR
  save,new_phi,new_y,filename='sav_data/phi_y_mz_interp.sav'

  
  stop
END
