# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:00:46 2020

@author: Azzurra
"""

#Density matrix calculations for photocurrent in a quantum cascade detector of 6 levels. Temperature, Electric field and spectral energy dependance.
#with some modifications in functions, it can be adapted to every numbers of levels.


import numpy as np
import matplotlib.pyplot as plt
from constants import *
from functions import fermi_stat, density_matrix_photo
from parameter import *

fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(10,7))
fig.tight_layout()



#Measurements conditions
I= 100                           #Light  Intensity [W cm-2]
T=np.array([78, 200, 250, 300])  # Temperature [K]
T_size=np.size(T)
SpectralEnergy= np.arange(95, 170, 1) #Energy spectrum [meV]
En=np.size(SpectralEnergy)   #lenght of the array of the spectral energies



#Parameters of the quantum structure
Period=0
i=0
for i in range(np.size(struc_ini)):
     Period+=struc_ini[i] # period of the structure [Angstrom]

N=np.size(Levels_array,axis=0) #number of levels in the detector 
M=int(N*(N-1)/2)  #number of transitions
n=np.size(Energies[:,0:1]) # size of the electric field
Efermi, PopEquili_0=fermi_stat(Energies, n_s, m_eff, N, n, T) #chemical potential [meV],population at equilibrium normalized to 1
n_LO=1/(np.exp(E_LO/(k_B*T))-1) #Phonon occupation factor



#Density matrix calculation
Population=np.zeros([T_size], dtype=object)
Coherences=np.zeros([T_size], dtype=object)
for t in range(T_size):
    U, Population[t], Coherences[t] = density_matrix_photo(struc_ini, Levels_array, t_dephasing, Energies, Osc_str, SpectralEnergy, Delta_sep, n_r, n_d, n_s, m_eff,  Gamma1, N, M, n, En, n_LO[t], I, z_elem, E_LO, T[t], Times_LO, Efermi[t,:], Times_other, PopEquili_0[t,:,:])


#Photocurrent calculation
j= np.zeros([T_size], dtype=object)
J=np.zeros([T_size], dtype=object)
Signal=np.zeros([T_size], dtype=object)
for t in range(T_size):
    j[t]=np.zeros([n,En], dtype=object)
    J[t]=np.zeros([n,En], dtype=object)
    Signal[t]=np.zeros([n], dtype=object)
    for k in range(n):
        for e in range(En):                
            j[t][k,e]=np.real(-n_d*e_C*1j*U[k,e]*Coherences[t][k,e]*1e12) #[A cm^-2] current density for every levels
            J[t][k,e]=(j[t][k,e][5,0]+j[t][k,e][6,0] +j[t][k,e][9,0])     #total photocurrent from levels with non zero coherences
            Signal[t][k]=np.amax(J[t][k,:])  #maximum photocurrent at spectral peak
            
#Plot
 
    axs[0].plot(Energies[:,0],Signal[t][:]/I,'o')         #Responsivity vs electric field [A/W] for every temperatures
    axs[1].plot(SpectralEnergy,J[t][1,:],'o')             #Photocurrent spectrum for every temperature
    axs[2].plot(T[t],Signal[t][1]/I,'o')                  #Responsivity vs temperature [A/W] for a fixed electric field

    axs[2].set_xlabel('Temperature [K]')
    axs[1].set_xlabel('Spectral Energy [meV]')    
    axs[0].set_xlabel('Electric field [kV/cm]') 

    axs[2].set_ylabel('Responsivity [A/W]')
    axs[1].set_ylabel('Photocurrent density [A cm-2]')
    axs[0].set_ylabel('Responsivity [A/W]')