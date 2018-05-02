PRO load_cosmology, set_omega_l=set_omega_l, $
                    set_omega_m=set_omega_m, $
                    set_omega_b=set_omega_b, $
                    set_omega_cdm=set_omega_cdm, $
                    set_omega_r=set_omega_r, $
                    set_h=set_h, $
                    set_spec_ind=set_spec_ind, $
                    set_sigma8=set_sigma8
                    
  COMMON cosmological_parameters, $
     omega_l,omega_m,omega_b, $
     omega_cdm, omega_r, $
     h, spec_ind, sigma8

  
  IF (n_elements(set_omega_l) EQ 0) THEN omega_l=0.725 ELSE $
     omega_l=set_omega_l
  IF (n_elements(set_omega_m) EQ 0) THEN omega_m=0.275 ELSE $
     omega_m=set_omega_m
  IF (n_elements(set_omega_b) EQ 0) THEN omega_b=0.046 ELSE $
     omega_b=set_omega_b
  IF (n_elements(set_omega_l) EQ 0) THEN omega_cdm=omega_m-omega_b ELSE $
     omega_cdm=set_omega_cdm
  IF (n_elements(set_omega_r) EQ 0) THEN omega_r=8.4d-5 ELSE $
     omega_r=set_omega_r
  IF (n_elements(set_h) EQ 0) THEN h=0.702 ELSE $
     h=set_h
  IF (n_elements(set_spec_ind) EQ 0) THEN spec_ind=0.96 ELSE $
     spec_ind=set_spec_ind
  IF (n_elements(set_sigma8) EQ 0) THEN sigma8=0.82 ELSE $
     sigma8=set_sigma8

  return
  
END
