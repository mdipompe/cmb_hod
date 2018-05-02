from hmf import MassFunction
from hmf import cosmo
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from astropy.cosmology import FlatLambdaCDM
import pdb
import math

Mmax=15.72
Mmin=5.
Om0=0.275
H0=70.2
hmf_model='Tinker10'
delta_wrt='mean'
delta_h=1.686
sigma_8=0.82
n=0.96
Tcmb0=2.7255

Lambda=0.725

#zs=[0.005, 5.00851, 9.99500]
zs=np.logspace(np.log10(0.005),np.log10(9.995),num=100,endpoint=True)

E_z=np.sqrt((Om0*(1.+zs)**3.) + Lambda)
Omz=Om0*((1.+zs)**3.)/(E_z**2.)

Delta_cr=(18.*np.pi**2) + (82*(Omz-1.)) - (39.*(Omz-1)**2.)
deltas=Delta_cr/Omz

cosmo_model = FlatLambdaCDM(H0 = H0,
                            Om0=Om0,
                            Tcmb0 = Tcmb0,
                            Ob0 = 0.046)


hmf = MassFunction(Mmax=15.72,
                   Mmin=5.,
                   cosmo_model=cosmo_model,
                   hmf_model='Tinker10',
                   delta_wrt='mean',
                   delta_h=200.,
                   sigma_8=0.82,
                   n=0.96)

i=0
for i in range(len(zs)-1):
    hmf.update(z=zs[i],delta_h=deltas[i])

    mass_func = hmf.dndm

#    outfile1='hmf_dndm_z'+str(zs[i])+'_200.txt'
    outfile1='hmf_dndm_z'+str(zs[i])+'.txt'
    dndmfile=open(outfile1,'w')
#    outfile2='hmf_pspec_z'+str(zs[i])+'_200.txt'
    outfile2='hmf_pspec_z'+str(zs[i])+'.txt'
    pspecfile=open(outfile2,'w')
    
    for i in range(0,len(hmf.m)-1):
        dndmfile.write(str(hmf.m[i]) + '  ' +
                       str(hmf.sigma[i]) + '  ' +
                       str(hmf.fsigma[i]) + '  ' +
                       str(hmf.dndm[i])+'\n')
        

    for i in range(0,len(hmf.k)-1):
        pspecfile.write(str(hmf.k[i]) + '  ' +
                        str(hmf.power[i])+'\n')

    dndmfile.close()
    pspecfile.close()


#fig=plt.figure(figsize=(8,5))
#ax=fig.add_subplot(1,1,1)
#plt_dndm, = ax.plot(hmf.m,hmf.dndm,'-',color='grey')
#ax.set_yscale('log')
#ax.set_xscale('log')
#plt.xlabel(r'$M_{h}$',fontsize=18)
#plt.xlabel(r'$dN/dM$')
#ax.set_xlim(1e5,1e16)
#ax.grid(True)
#plt.show()

#file="../../power_spec/camb_matterpower_0.00500000.dat"
#dtype=np.dtype([('k', 'f'), ('pk','f')])
#pspec = np.loadtxt(file, dtype=dtype)

#fig=plt.figure(figsize=(10,7))
#ax=fig.add_subplot(1,1,1)
#plt_power, = ax.plot(hmf.k,hmf.power,'-',color='grey')
#plt_mypower, = ax.plot(pspec['k'],pspec['pk'],'--',color='red')                          
#ax.set_yscale('log')
#ax.set_xscale('log')
#plt.xlabel(r'$k$',fontsize=18)
#plt.xlabel(r'$P(k)$')
#ax.grid(True)
#plt.show()

