PRO compare


  ell=cgloggen(1000,start=1,finish=10000)
  restore,'c_l_2h.sav'
  restore,'c_l_1h.sav'
  readcol,'dipompeo_2017_xcorr_model.txt',$
          ell_old,cl_old


  plot,ell,c_l_2h,linestyle=0,thick=2,$
       xtit='ell',ytit=textoidl('C_l'),charsize=2,$
       xra=[1,5000],/xlog,yra=[1d-30,1d-1],/ylog

  oplot,ell,c_l_1h,linestyle=2,thick=2

  oplot,ell,c_l_1h+c_l_2h,linestyle=2,thick=2,color=cgcolor('yellow')
  
  oplot,ell_old,cl_old,linestyle=2,color=cgcolor('red'),thick=2



  
  stop
END
