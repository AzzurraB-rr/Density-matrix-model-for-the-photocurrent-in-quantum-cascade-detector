# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 12:00:46 2020

@author: Azzurra
"""

#pROGRAM for simplified densty matrix in localized basis
#form taken from appendix C.2 of Giulia Pegolotti thesis applied to a QCD



import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import inv
import itertools

 
Period=396*1e-8 #period in Angstrom for in cm
e_C=1.6e-19 #in C
h=6.582e-1 #h bar in meV*ps
hbar = 1.0545718e-1
n_d=1e18  #electron density in cm-3
n_s=0.5 #[1e12 cm-2] areal electron density
#Omega=6.76/(2*h) #1/ps
eps_0=8.8542e-12  #F m-1
c=3e8 #m s-1
n_r=3.3 #refractive index GaAs
I=200 #W cm-2
#t_abs=100 #ps
Gamma1=6 #meV
m_eff=0.068*9.11e-31 #kg
f_21=0.96
T= 78#K
k_B=8.62*1e-2 #[meV K]
#Insert levels_struc parameters
E_LO=36 #meV

Efermi=22
#(3.14*np.power(h_bar,2)*n_s*1e4)/(m_eff*e_C*1e-3) # en meV
n_LO=1/(np.exp(E_LO/(k_B*T))-1)
struc_ini = [70,67,18,46,22,38,30,33,43,32] #in Angstrom
Levels=np.array([1,2,3,4,5,6])
N=np.size(Levels,axis=0) #number of levels 
M=int(N*(N-1)/2)  #number of transitions

#%%
'''
Inputs calculated with "Localized_basis_omega_wavefunction_energies_times.py"

#INPUT: Energies vs Electric field for every levels, 
#       Scattering Times for chosen consecutive transitions,
#       Dipoles calculated in localized basis
txt files input
'''
#------------------------
with open('energies_'+str(struc_ini[2])+'_'+str(struc_ini[4])+'_'+str(struc_ini[6])+'_'+str(struc_ini[8])+'.txt', 'r') as data:
     Energies=np.genfromtxt(data)[:]        #Energies files 0 are Efield, E1,E2,E3,E4,E5,E6, E1'
#Energies[:,1:]=Energies[:,1:]*1e3 #Converting to meV
#Energies=Energies[:,0:N+1]    #removes the last values of energies since it s not used
n=np.size(Energies[:,0:1]) #for all the electric field in data1,0 values

with open('Times_tot_'+str(struc_ini[2])+'_'+str(struc_ini[4])+'_'+str(struc_ini[6])+'_'+str(struc_ini[8])+'.txt', 'r') as data_t:
     Times_tot=np.genfromtxt(data_t)[:]  #times are in ps #field, t21, t61, t56, t45, t34, '
     
with open('z_elem_'+str(struc_ini[2])+'_'+str(struc_ini[4])+'_'+str(struc_ini[6])+'_'+str(struc_ini[8])+'.txt', 'r') as data_z:
     z_elem=1e-8*np.genfromtxt(data_z)[:] #Dipoles are in A, z21, z,61, z56, z45, z34
     
with open('Osc_str_'+str(struc_ini[2])+'_'+str(struc_ini[4])+'_'+str(struc_ini[6])+'_'+str(struc_ini[8])+'.txt', 'r') as data_osc:
     Osc_str=np.genfromtxt(data_osc)[:] #Oscillator strenght f21,f61, f56, f45, f34
     
with open('Omega_L1437.txt', 'r') as data_om:
     Omega_ext=np.genfromtxt(data_om)[:] #Values of Omega, rabi frequencies in 1/ps, calculated for extended basis, the activated couplng are 23 and 24

with open('PopulationEquilibrium.txt', 'r') as popula:
     PopEquili=np.genfromtxt(popula)[:,1:]     

Transitions=list(itertools.combinations(Levels, 2)) #combnations of all the transitions
alfa_p=np.zeros([M],int)
alfa_q=np.zeros([M],int)
for p in range(M):
    alfa_p[p]=Transitions[p][0]-1  #it starts from 0 to take account how Python counts in array (from 0)
for q in range(M):
    alfa_q[q]=Transitions[q][1]-1
    

#------------------------------------------------------------------------------------    
SpectralEnergy= np.arange(80, 180, 1) #Energy spectrum in meV
En=np.size(SpectralEnergy)   #lenght of the array of the spectral energies
#-------------------------------------
t_dephasing= np.arange(0.1, 0.2, 0.1) #dephasing time parameter in ps
ndeph=np.size(t_dephasing)   #lenght of the array of the dephasing time

#oscillator strenghts
f_21=np.zeros(n)
f_31=np.zeros(n)  
for k in range(n):
    f_21[k]=Osc_str[k,1]
    f_31[k]=Osc_str[k,2]
    


#----------------------------------------------------------
Omega=np.zeros([n, ndeph, En],dtype=object)  #All the the tunnel couplings
Delta=np.zeros([n, ndeph, En],dtype=object) #All energies differences
Tau_dephasing=np.zeros([n, ndeph,En], dtype=object)
PEqu=np.zeros([n, ndeph,En], dtype=object)  
for l in range(ndeph) :  
    for k in range(n):
        for e in range(En):
            Delta[k,l,e]=np.zeros([N,N])
            Omega[k,l,e]=np.zeros([N,N])
            Tau_dephasing[k,l,e]=np.zeros([N,N])
            PEqu[k,l,e]=np.zeros([N,1],complex)
            for i in range(N):
                for j in range(N):
                    Delta[k,l,e][i,j]=np.abs(Energies[k,i+1]-Energies[k,j+1])
                    Tau_dephasing[k,l,e][i,j]=(i!=j)*0.01
                    PEqu[k,l,e][i]=PopEquili[k,i]
                    if (i==1 and j==2) or (i==2 and j==1):                  #the tunnel couplings different from 0
                         Omega[k,l,e][i,j] = 6.76/(h) # Omega_ext[0]
                         Tau_dephasing[k,l,e][i,j]=t_dephasing[l]     
              #      if (i==1 and j==3) or (i==3 and j==1):
              #           Omega[k,l,e][i,j] = 2.483/(h) #Omega_ext[1] maybe wrong 2.483/(2*h)
              #           Tau_dephasing[k,l,e][i,j]=2.7
 
#Depolarization shift from Schneider book, with epsilon=epsilon_0*n_r^2
E_dep_21=np.zeros(n)
for k in range(n):
    E_dep_21[k]=np.abs(Energies[k,1]-Energies[k,2])+h*h*e_C*e_C*n_d*1e-18*f_21[k]/(2*n_r*n_r*eps_0*m_eff*np.abs(Energies[k,1]-Energies[k,2]))
   

#Lorenztian absorption function
#-------------------------------------------------------------
t_abs_21=np.zeros([n, En],dtype=object)
for i in range(En) :  
    for k in range(n): 
        t_abs_21[k,i]=(SpectralEnergy[i]*1e12*1e12)/(I*6.24*1e18*1e3*0.15*f_21[k])*(((np.abs(E_dep_21[k])-SpectralEnergy[i])**2+Gamma1*Gamma1)/Gamma1)                    

#Single particle matrix Hamiltonian
H_sp=np.zeros([n, ndeph, En],dtype=object)
for l in range(ndeph) :  
    for k in range(n):
        for e in range(En):
            H_sp[k,l,e]=np.zeros([N,N],complex)
            for i in range(N):
                for j in range(N):
                    H_sp[k,l,e][i,j]= (i==j)*Energies[k,j+1]
                    H_sp[k,l,e][i,j]+= (i!=j)*Omega[k,l,e][i,j]


#Dipoles
Z=np.zeros([n, ndeph,En], dtype=object) 
for l in range(ndeph) :  
    for k in range(n):
        for e in range(En):
            Z[k,l,e]=np.zeros([N,N],complex)
            for i in range(N):
                for j in range(N):
                     Z[k,l,e][i,j]+=((i==0 and j==1) or (i==1 and j==0))*-np.abs(z_elem[k,1])
                     Z[k,l,e][i,j]+=((i==2 and j==3) or (i==3 and j==2))*-np.abs(z_elem[k,5])
                     Z[k,l,e][i,j]+=((i==3 and j==4) or (i==4 and j==3))*-np.abs(z_elem[k,4])
                     Z[k,l,e][i,j]+=((i==4 and j==5) or (i==5 and j==4))*-np.abs(z_elem[k,3])
                     Z[k,l,e][i,j]+=((i==5 and j==0) or (i==0 and j==5))*-np.abs(z_elem[k,2])
                     Z[k,l,e][i,j]+=(i==j==2)*-(struc_ini[0]/2+struc_ini[1]+struc_ini[2]/2)*1e-8  
                     Z[k,l,e][i,j]+=(i==j==3)*-((struc_ini[0]/2+struc_ini[1]+struc_ini[2]/2)+(struc_ini[2]/2+struc_ini[3]+struc_ini[4]/2))*1e-8
                     Z[k,l,e][i,j]+=(i==j==4)*-((struc_ini[0]/2+struc_ini[1]+struc_ini[2]/2)+(struc_ini[2]/2+struc_ini[3]+struc_ini[4]/2)+(struc_ini[4]/2+struc_ini[5]+struc_ini[6]/2))*1e-8
                     Z[k,l,e][i,j]+=(i==j==5)*-((struc_ini[0]/2+struc_ini[1]+struc_ini[2]/2)+(struc_ini[2]/2+struc_ini[3]+struc_ini[4]/2)+(struc_ini[4]/2+struc_ini[5]+struc_ini[6]/2)+(struc_ini[6]/2+struc_ini[7]+struc_ini[8]/2))*1e-8
                 
#Writing the matrix to solve the equations:
U_total=np.zeros([n, ndeph,En],dtype=object)
U=np.zeros([n, ndeph,En],dtype=object)
for l in range(ndeph) :  
    for k in range(n):
        for e in range(En):
            U_total[k,l,e]=np.zeros([N,N],complex) 
            U_total[k,l,e]=(1/h)*(H_sp[k,l,e] @ Z[k,l,e]-Z[k,l,e] @ H_sp[k,l,e]) #commutator
            U[k,l,e]=np.zeros([M,1],complex)
            m=0
            for i in range(N):
                for j in range(i+1,N):
                    U[k,l,e][m,0]=U_total[k,l,e][i,j]
                    m=m+1
               
                
C=np.zeros([n, ndeph,En],dtype=object)                    
for l in range(ndeph) :  
    for k in range(n):
        for e in range(En):
            C[k,l,e]=np.zeros([2*M,N],complex)
            for j in range(N):
                for a in range(M):
                    C[k,l,e][a+M,j]=(j==alfa_q[a])*2*1j*Omega[k,l,e][alfa_p[a],alfa_q[a]]
                    C[k,l,e][a+M,j]+=(j==alfa_p[a])*(-2)*1j*(Omega[k,l,e][alfa_p[a],alfa_q[a]])


D=np.zeros([n, ndeph,En],dtype=object)
for l in range(ndeph) :  
    for k in range(n):
        for e in range(En):
            D[k,l,e]=np.zeros([2*M,2*M],complex)
            for a in range(M):
                for b in range(M):
                    D[k,l,e][a,b]=(a==b)*1j*Delta[k,l,e][alfa_p[a],alfa_q[a]]/h
                    D[k,l,e][4,4]=1j*(Energies[k,6]-Energies[k,7])/h
                    D[k,l,e][a,b]+=(alfa_p[a]==alfa_p[b] and alfa_q[a]!=alfa_q[b])*-1j*Omega[k,l,e][alfa_q[a],alfa_q[b]]
                    D[k,l,e][a,b]+=(alfa_p[a]==alfa_q[b] and alfa_q[a]!=alfa_p[b])*1j*Omega[k,l,e][alfa_q[a],alfa_p[b]]
                    D[k,l,e][a,b]+=(alfa_q[a]==alfa_p[b] and alfa_p[a]!=alfa_q[b])*1j*-Omega[k,l,e][alfa_p[a],alfa_q[b]]
                    D[k,l,e][a,b]+=(alfa_q[a]==alfa_q[b] and alfa_p[a]!=alfa_p[b])*1j*Omega[k,l,e][alfa_p[a],alfa_p[b]]
                    D[k,l,e][a+M,b+M]+=(a==b)*1j*Delta[k,l,e][alfa_p[a],alfa_q[a]]/h
                    D[k,l,e][19,19]=1j*(Energies[k,6]-Energies[k,7])/h
                    D[k,l,e][a+M,b+M]+=((alfa_p[a]==alfa_p[b]) and (alfa_q[a]!=alfa_q[b]))*1j*-Omega[k,l,e][alfa_q[a],alfa_q[b]]
                    D[k,l,e][a+M,b+M]+=(alfa_p[a]==alfa_q[b] and alfa_q[a]!=alfa_p[b])*1j*-Omega[k,l,e][alfa_q[a],alfa_p[b]]
                    D[k,l,e][a+M,b+M]+=(alfa_q[a]==alfa_p[b] and alfa_p[a]!=alfa_q[b])*1j*Omega[k,l,e][alfa_p[a],alfa_q[b]]
                    D[k,l,e][a+M,b+M]+=(alfa_q[a]==alfa_q[b] and alfa_p[a]!=alfa_p[b])*1j*Omega[k,l,e][alfa_p[a],alfa_p[b]]
                    D[k,l,e][a,b+M]+=(a==b)*1/(Tau_dephasing[k,l,e][alfa_p[a],alfa_q[a]])
                    D[k,l,e][a+M,b]+=(a==b)*1/(Tau_dephasing[k,l,e][alfa_p[a],alfa_q[a]])

E=np.zeros([n, ndeph,En],dtype=object)
for l in range(ndeph):  
    for k in range(n):
        for e in range(En):
            E[k,l,e]=np.zeros([N,M],complex)
            for i in range(N):
                for a in range(M):
                    E[k,l,e][i,a]=(i==alfa_p[a])*(-1j)*Omega[k,l,e][alfa_p[a],alfa_q[a]]
                    E[k,l,e][i,a]+=(i==alfa_q[a])*(1j)*Omega[k,l,e][alfa_p[a],alfa_q[a]]
                   

                  

F_abs=np.zeros([n,ndeph,En],dtype=object)
for l in range(ndeph):
    for k in range(n):
        for e in range(En):
            F_abs[k,l,e]=np.zeros([N,N], complex)
            for i in range(N):
                for j in range(N):
                    F_abs[k,l,e][i,j]=(i==j==0)*(-1/(t_abs_21[k,e]))
                    F_abs[k,l,e][i,j]+=(i==j==1)*(-1/(t_abs_21[k,e]))
                  #  F_abs[k,l,e][i,j]+=(i==j==2)*((-1/Times_tot[k,5])*(1+n_LO))
                  #  F_abs[k,l,e][i,j]+=(i==j==3)*(-(1/Times_tot[k,4])*(1+n_LO)-(1/Times_tot[k,5])*(n_LO))
                  #  F_abs[k,l,e][i,j]+=(i==j==4)*(-(1/Times_tot[k,3])*(1+n_LO)-(1/Times_tot[k,4])*(n_LO))
                  #  F_abs[k,l,e][i,j]+=(i==j==5)*(-(1/Times_tot[k,2])*(1+n_LO)-(1/Times_tot[k,3])*(n_LO)) 
                    F_abs[k,l,e][i,j]+=(i==0 and j==1)*(+1/(t_abs_21[k,e]))  
                    F_abs[k,l,e][i,j]+=(i==1 and j==0)*(+1/(t_abs_21[k,e]))#+(1/Times_tot[k,1])*(n_LO)) 
                
                    
                    
F=np.zeros([n,ndeph,En],dtype=object)
for l in range(ndeph):
    for k in range(n):
        for e in range(En):
            F[k,l,e]=np.zeros([N,N], complex)
            for i in range(N):
                for j in range(N):
                    F[k,l,e][i,j]=(i==j==0)*(-(1/Times_tot[k,2])*n_LO*(1/(np.exp((-Efermi+Energies[k,6]-Energies[k,7]-E_LO)/(k_B*T))+1))-(1/Times_tot[k,1])*(n_LO)*(1/(1+np.exp((-Efermi+Delta[k,l,e][0,1]-E_LO)/(k_B*T))))-1/(t_abs_21[k,e]))#*(1+np.exp(((-Energies[k,0]*Period*1e-2))/(k_B*T))))
                    F[k,l,e][i,j]+=(i==j==1)*((-1/Times_tot[k,1])*(1+n_LO)-1/(t_abs_21[k,e]))
                    F[k,l,e][i,j]+=(i==j==2)*((-1/Times_tot[k,5])*(1+n_LO))
                    F[k,l,e][i,j]+=(i==j==3)*(-(1/Times_tot[k,4])*(1+n_LO)-(1/Times_tot[k,5])*(n_LO)*(1/(np.exp((Delta[k,l,e][2,3]-E_LO)/(k_B*T))+1)))
                    F[k,l,e][i,j]+=(i==j==4)*(-(1/Times_tot[k,3])*(1+n_LO)-(1/Times_tot[k,4])*(n_LO)*(1/(np.exp((Delta[k,l,e][3,4]-E_LO)/(k_B*T))+1)))
                    F[k,l,e][i,j]+=(i==j==5)*(-(1/Times_tot[k,2])*(1+n_LO)-(1/Times_tot[k,3])*(n_LO)*(1/(np.exp((Delta[k,l,e][4,5]-E_LO)/(k_B*T))+1))) 
                    F[k,l,e][i,j]+=(i==0 and j==1)*((1/Times_tot[k,1])*(1+n_LO)+1/(t_abs_21[k,e]))  
                    F[k,l,e][i,j]+=(i==1 and j==0)*(+1/(t_abs_21[k,e])+(1/Times_tot[k,1])*(n_LO)*(1/(np.exp((-Efermi+Delta[k,l,e][0,1]-E_LO)/(k_B*T))+1))) 
                    F[k,l,e][i,j]+=((i==3 and j==2))*(1/Times_tot[k,5])*(1+n_LO)
                    F[k,l,e][i,j]+=((i==2 and j==3))*(1/Times_tot[k,5])*(n_LO)*(1/(np.exp((Delta[k,l,e][2,3]-E_LO)/(k_B*T))+1))
                    F[k,l,e][i,j]+=((i==4 and j==3))*(1/Times_tot[k,4])*(1+n_LO)
                    F[k,l,e][i,j]+=((i==3 and j==4))*(1/Times_tot[k,4])*(n_LO)*(1/(np.exp((Delta[k,l,e][3,4]-E_LO)/(k_B*T))+1)) 
                    F[k,l,e][i,j]+=((i==5 and j==4))*(1/Times_tot[k,3])*(1+n_LO)
                    F[k,l,e][i,j]+=((i==4 and j==5))*(1/Times_tot[k,3])*(n_LO)*(1/(np.exp((Delta[k,l,e][4,5]-E_LO)/(k_B*T)+1)))
                    F[k,l,e][i,j]+=((i==0 and j==5))*((1/Times_tot[k,2])*(1+n_LO))#*np.exp((-Energies[k,0]*Period*1e-2)/(k_B*T))) 
                    F[k,l,e][i,j]+=((i==5 and j==0))*((1/Times_tot[k,2])*n_LO*(1/(np.exp((-Efermi+Energies[k,6]-Energies[k,7]-E_LO)/(k_B*T))+1)))           
                    
                    
                    
                    
          #          F[k,l,e][i,j]=(i==j==0)*(-1/(t_abs_21[k,e])-(1/Times_tot[k,1])*(n_LO)-(1/Times_tot[k,2])*n_LO)#*(1+np.exp(((-Energies[k,0]*Period*1e-2))/(k_B*T))))
          #          F[k,l,e][i,j]+=(i==j==1)*((-1/Times_tot[k,1])*(1+n_LO)-1/(t_abs_21[k,e]))
          #          F[k,l,e][i,j]+=(i==j==2)*((-1/Times_tot[k,5])*(1+n_LO))
          #          F[k,l,e][i,j]+=(i==j==3)*(-(1/Times_tot[k,4])*(1+n_LO)-(1/Times_tot[k,5])*(n_LO))
          #          F[k,l,e][i,j]+=(i==j==4)*(-(1/Times_tot[k,3])*(1+n_LO)-(1/Times_tot[k,4])*(n_LO))
          #          F[k,l,e][i,j]+=(i==j==5)*(-(1/Times_tot[k,2])*(1+n_LO)-(1/Times_tot[k,3])*(n_LO)) 
          #          F[k,l,e][i,j]+=(i==0 and j==1)*((1/Times_tot[k,1])*(1+n_LO)+1/(t_abs_21[k,e]))  
          #          F[k,l,e][i,j]+=(i==1 and j==0)*(+1/(t_abs_21[k,e])+(1/Times_tot[k,1])*(n_LO)) 
          #          F[k,l,e][i,j]+=((i==3 and j==2))*(1/Times_tot[k,5])*(1+n_LO)
          #          F[k,l,e][i,j]+=((i==2 and j==3))*(1/Times_tot[k,5])*(n_LO)
          #          F[k,l,e][i,j]+=((i==4 and j==3))*(1/Times_tot[k,4])*(1+n_LO)
          #          F[k,l,e][i,j]+=((i==3 and j==4))*(1/Times_tot[k,4])*(n_LO) 
          #          F[k,l,e][i,j]+=((i==5 and j==4))*(1/Times_tot[k,3])*(1+n_LO)
          #          F[k,l,e][i,j]+=((i==4 and j==5))*(1/Times_tot[k,3])*(n_LO)
          #          F[k,l,e][i,j]+=((i==0 and j==5))*((1/Times_tot[k,2])*(1+n_LO))#*np.exp((-Energies[k,0]*Period*1e-2)/(k_B*T))) 
          #          F[k,l,e][i,j]+=((i==5 and j==0))*((1/Times_tot[k,2])*n_LO)


D_inv=np.zeros([n, ndeph,En],dtype=object)                
Y=np.zeros([n, ndeph,En],dtype=object)
Y_red=np.zeros([n, ndeph,En],dtype=object)
G=np.zeros([n, ndeph,En],dtype=object)
G_n=np.zeros([n, ndeph,En],dtype=object)
Population=np.zeros([n, ndeph,En],dtype=object)
PEquilibrium=np.zeros([n, ndeph,En],dtype=object)
Coherences = np.zeros([n, ndeph,En],dtype=object)
Vector= np.zeros([n, ndeph,En],dtype=object)
j= np.zeros([n, ndeph,En], dtype=object)
J=np.zeros([n, ndeph,En])
Signal=np.zeros([n,ndeph])
Pop=np.zeros([n,ndeph])
for l in range(ndeph):
    for k in range(n):
        for e in range(En):                
            Vector[k,l,e]=np.zeros([N,1],complex)
            D_inv[k,l,e]=inv(D[k,l,e])
            Y[k,l,e]=D_inv[k,l,e] @ C[k,l,e]
            Y_red[k,l,e]=Y[k,l,e][0:M,:]
            G[k,l,e]= E[k,l,e] @ Y_red[k,l,e]-(-F[k,l,e])#matrix F has an inverted sign respect to Giulia thesis
            G[k,l,e][(N-1):N,:]=1
            G_n[k,l,e]=-F_abs[k,l,e]-(E[k,l,e] @ Y_red[k,l,e])
            PEquilibrium[k,l,e]=G_n[k,l,e]@ PEqu[k,l,e]
            PEquilibrium[k,l,e][(N-1):N,:]=-1
            Vector[k,l,e][(N-1):N,:]=1
            Population[k,l,e]= inv(G[k,l,e]) @ (Vector[k,l,e]+PEquilibrium[k,l,e])
            Coherences[k,l,e]=Y_red[k,l,e] @ Population[k,l,e]# inv(G[k,l,e]) @ Vector[k,l,e]
      #      Population[k,l,e]=inv(G[k,l,e]) @ Vector[k,l,e]
            j[k,l,e]=np.real(-n_d*e_C*1j*U[k,l,e]*Coherences[k,l,e]*1e12)
            J[k,l,e]=j[k,l,e][5,0]+j[k,l,e][6,0] +j[k,l,e][9,0] #The current from the two levels where coherences are not 0
 #           J[k,l,e]=np.real(Population[k,l,e][4,0]/Times[k,4])
            Signal[k,l]+=J[k,l,e] 
    plt.plot(Energies[:,0],Signal[:,l],'o')
  #      plt.plot(SpectralEnergy,J[k,l,:],'o')


#Output: a file with the current density value, a file with the population for t dephasing 0,1 as funxion of the electric field                 
#Plot the populations
P=np.zeros([n, N],dtype=object)
for k in range(n):
    for i in range(N):
        P[k,i]=np.real(Population[k,0,50][i]) #Population at the peak


   
k=np.column_stack((Energies[:,0],Signal)) #aggiunge a j la colonna di e field
T_deph1=np.insert(t_dephasing,0,0,axis=0) #aggiunge uno 0 a t deph
AllValues=np.insert(k,[0],T_deph1,axis=0)
np.savetxt("IV_integrated_loc.txt",AllValues,fmt='%.10f')     

Popul=np.column_stack((Energies[:,0],P))
np.savetxt("Population.txt", Popul)

Current=np.zeros([n,En])
for k in range(En):
    for s in range(n):
        Current[s,k]=J[s,0,k]
        
Sig=np.column_stack((SpectralEnergy,Current.T))
Eg1=np.insert(Energies[:,0],0,0,axis=0)
Values=np.insert(Sig,[0],Eg1,axis=0)
np.savetxt("Spectrum_loc.txt",Values,fmt='%.9f')    

##op=np.column_stack((Energies[:,0],P))
#np.savetxt("Population.txt",Pop,fmt='%.6f')
