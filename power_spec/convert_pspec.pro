PRO convert_pspec

  ;Convert power spectra from P(k) to P(ell).
  ;Note also the power spectra from CAMB don't cover the full ell range
  ;at all z. Only on the high ell side, where things are a nice
  ;power law, so just fit and extrapolate...
  
  COMMON cosmological_parameters
  
  ell=cgloggen(1000,start=1,finish=10000)
  z=cgloggen(100,start=0.005,finish=9.995)
  d=cosmocalc(z,d_c=dc)
  dc=dc*h
  pspec=file_search('camb_matterpower*')

  FOR i=0L,n_elements(z)-1 DO BEGIN
     readcol,pspec[i],k,p,format='D'
     ell2=(k*dc[i])-0.5
     
     IF (min(ell) LT min(ell2)) THEN stop
     
     IF (max(ell2) LT max(ell)) THEN BEGIN
        plot,ell2,p,linestyle=0,/xlog,/ylog,xra=[1,10000],yra=[1d-6,1d4],xsty=1
        fit=linfit(alog10(ell2[n_elements(ell2)-100:n_elements(ell2)-1]),$
                   alog10(p[n_elements(ell2)-100:n_elements(ell2)-1]))
        newp=fit[0]+fit[1]*alog10(ell)
        xx=where(ell GT max(ell2))
        ell2=[ell2,ell[xx]]
        p=[p,10.^newp[xx]]
        quadterp,ell2,p,ell,p_out
        oplot,ell,p_out,linestyle=2,color=cgcolor('red')
;        stop
     ENDIF ELSE BEGIN
        quadterp,ell2,p,ell,p_out
     ENDELSE
     
     outfile='ell/'+strtrim(pspec[i],2)
     openw,lun,outfile,/get_lun
     FOR j=0L,n_elements(ell)-1 DO BEGIN
        printf,lun,strtrim(ell[j],2) + '   ' + strtrim(p_out[j],2),format='(A)'
     ENDFOR
     free_lun,lun
  ENDFOR
  stop
END
