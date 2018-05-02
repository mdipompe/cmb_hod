PRO go

  COMMON cosmological_parameters
  
  zarray=cgloggen(100,start=0.005,finish=9.995)
  paramfile='default_params.ini'

  root='camb'
  revz=reverse(zarray)
  numloop=floor((n_elements(revz)/150))
  IF (numloop NE 0) THEN BEGIN
     numrem=(n_elements(revz) MOD 150)
     FOR i=0L,numloop-1 DO BEGIN
        indx=indgen(150)+(i*150)
        matter_power_spec,paramfile,zarray[indx],h0=h0,omega_b=omega_b,$
                          omega_dm=omega_m-omega_b,omega_l=omega_l,maxk=50
     ENDFOR
     IF (numrem NE 0) THEN BEGIN
        matter_power_spec,paramfile,zarray[(n_elements(zarray)-numrem):n_elements(zarray)-1],$
                          h0=h,omega_b=omega_b,omega_dm=omega_m-omega_b,omega_l=omega_l,maxk=50
     ENDIF
  ENDIF ELSE BEGIN
     matter_power_spec,paramfile,zarray,h0=h,omega_b=omega_b,omega_dm=omega_m-omega_b,omega_l=omega_l,maxk=50
  ENDELSE
  
   
  stop
END
