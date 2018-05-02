PRO plot_final

  ;Plot final 1 and 2-halo C_ell, compared to H&S 2014 digitized versions

  ell=cgloggen(1000,start=1,finish=10000)
  
  restore,'sav_data/c_l_1h.sav'
  restore,'sav_data/c_l_2h.sav'

  ;Multiplicative factors. The 1d11 part is just to match the axis 
  ;scaling of H&S, but the mantissa should be 1 if things are right... 
  f1=2.2d11
  f2=2.9d11

  f1=1.55d11
  f2=1.2d11
  
  xtit='l'
  ytit=textoidl('l^2(l+1)C_l^{y\phi}/2\pi')
  plot,ell,((ell^2)*(ell+1))*c_l_1h/(2.*!dpi)*f1,/xlog,$
       xtit=xtit,ytit=ytit,charsize=2,xra=[10,10000],xsty=1,thick=2,yra=[0,4.],ysty=1
  oplot,ell,((ell^2)*(ell+1))*c_l_1h/(2.*!dpi)*f1,linestyle=0,thick=2,color=cgcolor('steel blue')
  oplot,ell,((ell^2)*(ell+1))*c_l_2h/(2.*!dpi)*f2,linestyle=0,thick=2,color=cgcolor('red')

  dir='hill_spergel_data/'
  readcol,dir+'one_halo_final.txt',hs_1h_ell,hs_1h_cl,format='D'
  oplot,hs_1h_ell,hs_1h_cl,linestyle=2,thick=2,color=cgcolor('steel blue')

  readcol,dir+'two_halo_final.txt',hs_2h_ell,hs_2h_cl,format='D'
  oplot,hs_2h_ell,hs_2h_cl,linestyle=2,thick=2,color=cgcolor('red')

  circsym,/fill
  legend,['1h*'+strtrim(f1/1d11,2),'2h*'+strtrim(f2/1d11,2)],psym=[8,8],$
         color=[cgcolor('steel blue'),cgcolor('red')],/top,/left,box=0,thick=2,charsize=2
  

  
  stop
END
