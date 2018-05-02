FUNCTION hod,logm,N_cent=N_cent,N_sat=N_sat,$
             m_cut=m_cut,m_min=m_min,sig_logM=sig_logM,m1=m1,alpha=alpha
;+
;  NAME:
;    hod
;
;  PURPOSE:
;    Get N(m) (N_central + N_satellite) for the halo model. Default
;    parameters from Richardson et al. (2013).
;
;  USE:
;    Nm=hod(logm[,N_cent=N_cent,N_sat=N_sat,$
;                 m_cut=m_cut,m_min=m_min,sig_logM=sig_logM,m1=m1,alpha=alpha])
;
;  INPUT:
;    logm - Log of the halo mass(es) (in M_sun/h). 
;
;  OPTIONAL INPUT:
;    m_cut - mass for exponential cutoff. Default: 10^5
;    m_min - minimum log(mass). Default: 13.65
;    sig_logM - Default: 0.78
;    m1 - power law mass normalization. Default: 10^14.32
;    alpha - power law index. Default: 2.59
;  
;  OPTIONAL KEYWORDS:
;    
;  OUTPUT:
;    Total N(M)
;    N_cent - number of centrals
;    N_sat - number of satellites
;
;  NOTES:
;
;  HISTORY:
;    4/30/18 - Written - MAD (Dartmouth)
;-
  
  ;MAD As in Richardson et al. (2013) and others, set M_cut to
  ;value << 10^12, such that exponential term is essentially 1.
  ;Only left in here for completeness and in case ever want to change.
  IF (n_elements(m_cut) EQ 0) THEN m_cut=10.^5    
  ;MAD Set other default HOD parameters (Richardson et al. 2013)
  IF (n_elements(m_min) EQ 0) THEN m_min=13.65
  IF (n_elements(sig_logM) EQ 0) THEN sig_logM=0.78
  IF (n_elements(M1) EQ 0) THEN M1=10.^14.32
  IF (n_elements(alpha) EQ 0) THEN alpha=2.59

  N_cent=0.5*(1+erf((logm-M_min)/sig_logM))
  N_sat=(((10.^logm)/M1)^alpha)*exp((-1.)*M_cut/(10.^logm))
  N_m=N_cent+N_sat

  return,N_m
END
