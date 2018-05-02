PRO go

  COMMON cosmological_parameters
  
  z=cgloggen(100,start=0.005,finish=9.995)
  logm=findgen(108)/10.+5
  pspec=file_search('../../power_spec/camb_matterpower*')

  readcol,'hmf_dndm_z0.005_200.txt',hmf_m1,hmf_sigma_1,hmf_f_1,hmf_dndm_1,format='D'
  readcol,'hmf_dndm_z5.00851_200.txt',hmf_m2,hmf_sigma_2,hmf_f_2,hmf_dndm_2,format='D'
  readcol,'hmf_dndm_z9.995_200.txt',hmf_m3,hmf_sigma_3,hmf_f_3,hmf_dndm_3,format='D'
;  readcol,'hmf_dndm_z0.005_200_t08.txt',hmf_m1,blah,blah,blah,hmf_f_1,format='D'
;  readcol,'hmf_dndm_z5.0_200_t08.txt',hmf_m2,blah,blah,blah,hmf_f_2,format='D'
;  readcol,'hmf_dndm_z9.9_200_t08.txt',hmf_m3,blah,blah,blah,hmf_f_3,format='D'
;  readcol,'hmf_dndm_z0.005_200_smt.txt',hmf_m1,blah,blah,blah,hmf_f_1,format='D'
;  readcol,'hmf_dndm_z5.0_200_smt.txt',hmf_m2,blah,blah,blah,hmf_f_2,format='D'
;  readcol,'hmf_dndm_z9.9_200_smt.txt',hmf_m3,blah,blah,blah,hmf_f_3,format='D'
  
  f1=halo_mass_function(alog10(hmf_m1), z[0],power_spec=pspec[0],/silent)
  f2=halo_mass_function(alog10(hmf_m2), z[90],power_spec=pspec[90],/silent)
  f3=halo_mass_function(alog10(hmf_m3), z[99],power_spec=pspec[99],/silent)

  s1=sigma_m(alog10(hmf_m1), z[0],power_spec=pspec[0],/silent)
  s2=sigma_m(alog10(hmf_m2), z[90],power_spec=pspec[90],/silent)
  s3=sigma_m(alog10(hmf_m3), z[99],power_spec=pspec[99],/silent)

;  f1=halo_mass_function(alog10(hmf_m1), z[0],power_spec=pspec[0],/silent,model='tinker08')
;  f2=halo_mass_function(alog10(hmf_m2), z[90],power_spec=pspec[90],/silent,model='tinker08')
;  f3=halo_mass_function(alog10(hmf_m3), z[99],power_spec=pspec[99],/silent,model='tinker08')
;  f1=halo_mass_function(alog10(hmf_m1), z[0],power_spec=pspec[0],/silent,model='smt')
;  f2=halo_mass_function(alog10(hmf_m2), z[90],power_spec=pspec[90],/silent,model='smt')
;  f3=halo_mass_function(alog10(hmf_m3), z[99],power_spec=pspec[99],/silent,model='smt')



  plot,hmf_m1,s1,linestyle=0,thick=2,/xlog,xra=[1d5,1d16],xsty=1,yra=[1d-2,10],ysty=1,charsize=2
  oplot,hmf_m2,s2,linestyle=2,thick=2
  oplot,hmf_m3,s3,linestyle=3,thick=2

  oplot,hmf_m1,hmf_sigma_1,linestyle=0,thick=2,color=cgcolor('red')
  oplot,hmf_m2,hmf_sigma_2,linestyle=2,thick=2,color=cgcolor('red')
  oplot,hmf_m2,hmf_sigma_3,linestyle=3,thick=2,color=cgcolor('red')

  stop

  plot,hmf_m1,hmf_sigma_1/s1,linestyle=0,thick=2,/xlog,xra=[1d5,1d16],yra=[0.95,1.05],xsty=1,ysty=1
  oplot,hmf_m2,hmf_sigma_2/s2,linestyle=2,thick=2
  oplot,hmf_m3,hmf_sigma_3/s3,linestyle=3,thick=2


  stop


  
  plot,hmf_m1,f1,linestyle=0,thick=2,/xlog,/ylog,xra=[1d5,1d16],xsty=1,yra=[1d-4,0.4],ysty=1,charsize=2
  oplot,hmf_m2,f2,linestyle=2,thick=2
  oplot,hmf_m3,f3,linestyle=3,thick=2

  oplot,hmf_m1,hmf_f_1,linestyle=0,thick=2,color=cgcolor('red')
  oplot,hmf_m2,hmf_f_2,linestyle=2,thick=2,color=cgcolor('red')
  oplot,hmf_m2,hmf_f_3,linestyle=3,thick=2,color=cgcolor('red')

  stop

  plot,hmf_m1,hmf_f_1/f1,linestyle=0,thick=2,/xlog,xra=[1d5,1d16],yra=[0.6,1.1],xsty=1,ysty=1
  oplot,hmf_m2,hmf_f_2/f2,linestyle=2,thick=2
  oplot,hmf_m3,hmf_f_3/f3,linestyle=3,thick=2


  stop
END
