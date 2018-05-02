FUNCTION diff_3d,x

  ;A somewhat hacky way to get the difference between
  ;adjacent values in a particular dimension of a multi-dimensional array. 
  ;Not very general, but set up specifically for the use cases where
  ;it is needed in the rest of these codes.
  
  s=size(x)
  IF (s[0] EQ 4) THEN BEGIN
     diff=fltarr(s[1],s[2],s[3]-1,s[4])
     FOR i=0L,n_elements(x[*,0,0,0])-1 DO BEGIN
        FOR j=0L,n_elements(x[0,*,0,0])-1 DO BEGIN
           diff[i,j,*,*]=x[i,j,1:s[3]-1,*]-x[i,j,0:s[3]-2,*]
        ENDFOR
     ENDFOR
  ENDIF

  IF (s[0] EQ 3) THEN BEGIN
     diff=fltarr(s[1],s[2],s[3]-1)
     FOR i=0L,n_elements(x[*,0,0])-1 DO BEGIN
        FOR j=0L,n_elements(x[0,*,0])-1 DO BEGIN
           diff[i,j,*]=x[i,j,1:s[3]-1]-x[i,j,0:s[3]-2]
        ENDFOR
     ENDFOR
  ENDIF


  IF (s[0] EQ 2) THEN BEGIN
     diff=fltarr(s[1],s[2]-1)
     FOR i=0L,n_elements(x[*,0])-1 DO BEGIN
        diff[i,*]=x[i,1:s[2]-1]-x[i,0:s[2]-2]
     ENDFOR
  ENDIF
  
  return,diff
END
