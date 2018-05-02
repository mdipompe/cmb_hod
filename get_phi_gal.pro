PRO get_phi_gal

  ;Get the phi_ell terms for the galaxy/agn distribution
  ;Right now this only differs from DM by setting type='num' 
  ;in all to phi_ell.
  
  z=cgloggen(100,start=0.005,finish=9.995)
  logm=findgen(108)/10.+5

  phi_gal_mz=fltarr(n_elements(logm),n_elements(z),1000)
  FOR i=0L,n_elements(z)-1 DO BEGIN
     counter,i,n_elements(z)
     phi_gal_mz[*,i,*]=phi_ell(logm,z[i],ell=ell_phi,type='num')
  ENDFOR
  save,phi_gal_mz,filename='sav_data/phi_gal_mz.sav'
;  restore,'sav_data/phi_gal_mz.sav'

  newm=findgen(1071)/100.+5
  new_phi_gal=fltarr(n_elements(newm),n_elements(z),n_elements(ell_phi))
  FOR i=0L,n_elements(z)-1 DO BEGIN
     counter,i,n_elements(z)
     FOR j=0L,n_elements(ell_phi)-1 DO BEGIN
        new_phi_gal[*,i,j]=spline(logm,phi_gal_mz[*,i,j],newm)
     ENDFOR
  ENDFOR
  save,new_phi_gal,filename='sav_data/phi_gal_mz_interp.sav'

  
  stop
END
