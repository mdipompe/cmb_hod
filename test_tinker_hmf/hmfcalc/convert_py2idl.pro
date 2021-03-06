PRO go

  z=cgloggen(100,start=0.005,finish=9.995)
  newm=findgen(1071)/100.+5
  
  files=file_search('hmf_dndm_z*')

  new_dndm=fltarr(n_elements(newm),n_elements(z))
  FOR i=0L,n_elements(files)-1 DO BEGIN
     counter,i,n_elements(files)
     readcol,files[i],m,sig,fsig,dndm,format='D'
     new_dndm[*,i]=dndm
  ENDFOR
  
  save,new_dndm,filename='../../sav_data/dndm_hmfcalc.sav'
  
  stop
END
