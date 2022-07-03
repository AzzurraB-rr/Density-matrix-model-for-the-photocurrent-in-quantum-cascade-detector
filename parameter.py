# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 16:52:07 2022

@author: Azzurra
"""


import numpy as np
import math
import itertools
from constants import k_B,e_C,h_bar_SI,h_bar_eV,c,eps_0, Kb_SI
from numpy.linalg import inv




struc_ini = [70,67,18,46,22,38,30,33,43,32]
Levels_array=np.array([1,2,3,4,5,6])
n_d=1e18  #electron density in cm-3
n_s=n_d*(struc_ini[0]-20)*1e-8 #[cm-2] areal electron density
n_r=3.3 #refractive index GaAs
Gamma1=8 #[ meV] linewidth of the optical transition
m_eff=0.068*9.11e-31 #[ kg]
f_21=0.96 #oscillator strenght doped well
E_LO=36 #meV Phonon energy for GaAs
t_dephasing= 0.1 #dephasing time parameter in ps 
Delta_sep=6.76 #[meV] Tunneling coupling between levels 2 and 3. Calculated separately.
Patch_antenna=5.5 #absorption enhancement due to the patch-antenna metamaterial (see references)
#%%
'''
#INPUT: Energies vs Electric field for every levels, 
#       Scattering Times for chosen consecutive transitions,
#       Dipoles calculated in localized basis
txt files input
'''
#------------------------
with open('energies_'+str(struc_ini[2])+'_'+str(struc_ini[4])+'_'+str(struc_ini[6])+'_'+str(struc_ini[8])+'.txt', 'r') as data:
     Energies=np.genfromtxt(data)[:]        #Energies files 0 are Efield, E1,E2,E3,E4,E5,E6, E1'
#Energies[:,1:]=Energies[:,1:]*1e3 #Converting to meV

with open('Times_tot_'+str(struc_ini[2])+'_'+str(struc_ini[4])+'_'+str(struc_ini[6])+'_'+str(struc_ini[8])+'.txt', 'r') as data_t:
      Times_tot=np.genfromtxt(data_t)[:]  #times are in ps #field, t21, t61, t56, t45, t34, '
     
with open('Times_LO_'+str(struc_ini[2])+'_'+str(struc_ini[4])+'_'+str(struc_ini[6])+'_'+str(struc_ini[8])+'.txt', 'r') as data_t:
     Times_LO=np.genfromtxt(data_t)[:]  #times are in ps #field, t21, t61, t56, t45, t34, '

with open('Times_other_'+str(struc_ini[2])+'_'+str(struc_ini[4])+'_'+str(struc_ini[6])+'_'+str(struc_ini[8])+'.txt', 'r') as data_t:
     Times_other=np.genfromtxt(data_t)[:]  #times are in ps #field, t21, t61, t56, t45, t34,
     
with open('z_elem_'+str(struc_ini[2])+'_'+str(struc_ini[4])+'_'+str(struc_ini[6])+'_'+str(struc_ini[8])+'.txt', 'r') as data_z:
     z_elem=1e-8*np.genfromtxt(data_z)[:] #Dipoles are in Angstrom, z21, z,61, z56, z45, z34
     
with open('Osc_str_'+str(struc_ini[2])+'_'+str(struc_ini[4])+'_'+str(struc_ini[6])+'_'+str(struc_ini[8])+'.txt', 'r') as data_osc:
     Osc_str=np.genfromtxt(data_osc)[:] #Oscillator strenght f21,f61, f56, f45, f34
     
with open('Omega_L1437.txt', 'r') as data_om:
     Omega_ext=np.genfromtxt(data_om)[:] #Values of Omega, rabi frequencies in 1/ps, calculated for extended basis, the  couplng are 23 and 24    
