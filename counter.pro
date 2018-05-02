;+
;  NAME:
;    counter
;  PURPOSE:
;    prints the current element a loop is working on, overwriting the
;    last one so that output doesn't build on the screen.  On
;    the last element will leave the counter displayed.
;
;  USE:
;    counter,i,n_elements(data)
;
;  INPUT:
;    index - current index of loop
;    total - number of elements looping over
;    n - number of points in the mock sample
;
;  OPTIONAL INPUT:
;    str - a string that will display at the end of each count (so you
;          could add a list of e.g. object names, and input one of
;          these list elements to counter each time so you know which
;          object, not just number, it is working on).
;
;  KEYWORDS:        
;
;  OUTPUT:
;    prints to screen
;
;  NOTES:
;
;  HISTORY:
;    10-31-14 - Written - MAD (UWyo)
;    11-06-14 - Added str option - MAD (UWyo)
;    2015 - Cleaned and documented - MAD (UWyo)
;-
PRO counter,index,total,str=str

IF ~keyword_set(str) THEN BEGIN
   IF (index+1 NE total) THEN BEGIN
      print, format = '("Working on ",i8," of ",i8,a1,$)', $
             index+1,total,string(13B)
   ENDIF ELSE BEGIN
      print, format = '("Working on ",i8," of ",i8,"...counter finished")', $
             index+1,total
   ENDELSE
ENDIF ELSE BEGIN
   IF (index+1 NE total) THEN BEGIN
      print, format = '("Working on ",i8," of ",i8,A10,a1,$)', $
             index+1,total,str,string(13B)
   ENDIF ELSE BEGIN
      print, format = '("Working on ",i8," of ",i8,A10,"...counter finished")', $
             index+1,total,str
   ENDELSE
ENDELSE

END
