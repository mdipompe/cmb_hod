PRO go

  COMMON cosmological_parameters


  pspec='../power_spec/camb_matterpower_0.00500000.dat'
;  pspec='../camb_matterpower_5.00851.dat'
  readcol,'hmfcalc/hmf_dndm_z0.005.txt',hmf_m,hmf_sig,hmf_f,hmf_dndm,format='D'
  readcol,'hmfcalc/hmf_pspec_z0.005.txt',hmf_k,hmf_pk,format='D'
  readcol,pspec,k,pk,format='D'

;  pspec='hmf_pspec_z0.005.txt'

  readcol,'dndm_norm.txt',nz,norm,format='D'

  
  z=0.005
;  z=5.00851
;  z=0.0
  
;  newm=findgen(108)/10.+5
  newm=findgen(1071)/100.+5
;  newm=findgen(2000)/100.+1.
  
  out=mhalo2bias(newm,z,power_spec=pspec,model='tinker10',delta=200.)
  nu=out[*,0]
  b=out[*,1]

  out2=halo_mass_function(newm,z,power_spec=pspec,model='tinker10',delta=200.)
  nu2=out2[*,0]
  fnu=out2[*,1]
  dndm=out2[*,2]

;  xx=closest(nz,z)
;  norm=norm[xx]
;  fnu=fnu/norm
  
  
  sig=1.686/nu2

  plot,alog10(nu),b,linestyle=0,xra=[-0.5,0.7],thick=2,xtit='log(nu)',ytit='b',charsize=2
  readcol,'t10_b.txt',t_nu,t_b,format='D'
  oplot,t_nu,t_b,linestyle=2,color=cgcolor('red'),thick=2


  print,int_tabulated(nu2,b*fnu,/double)
;  print,int_tabulated(nu2,b*(hmf_f/nu2),/double)

  print,int_tabulated(alog10(nu2),b*fnu*alog(10)*nu2,/double)
;  print,int_tabulated(alog10(nu2),b*(hmf_f/nu2)*alog(10)*nu2,/double)

  step=diff(alog10(nu2))
  print,total(b*fnu*alog(10)*nu2*step)
  
  stop
END
