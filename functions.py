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



#Calculate autoconsistenly the chemical potential and level population for every temperature and electric field
def fermi_stat(Energies, n_s, m_eff, N, n, T):
    Ef_0K=(np.pi*np.power(h_bar_SI,2)*n_s*1e4)/(m_eff*e_C*1e-3)  # [meV] 
    L=np.size(T)
    n_j = np.zeros([L,n,N]) 
    n_tot=np.zeros([L,n])
    mu= np.ones([L,n])*Ef_0K   
    
    
    for i in range(n):
        for t in range(L):
           n_j[t,i,:] = 0. #for j in range(N)]
           flag = 1
           while flag == 1:
               n_j[t,i,0] = (m_eff*Kb_SI*T[t]/(np.pi*np.power(h_bar_SI,2)))*np.log(1+np.exp((mu[t,i])/(k_B*T[t])))
               for j in range(1,N):
                   n_j[t,i,j] = (m_eff*Kb_SI*T[t]/(np.pi*np.power(h_bar_SI,2)))*np.log(1+np.exp((mu[t,i]+Energies[i,7]-Energies[i, j+1])/(k_B*T[t])))
               n_tot[t,i] = np.sum(n_j[t,i,:])
               if n_tot[t,i] > n_s*1e4:
                  mu[t,i]-= 0.1
               else:    
                  flag = 0  
           n_j[t,i,:]=n_j[t,i,:]/n_tot[t,i]

    return mu, n_j  #chemical potential, level normalized population


#Implement the density matrix equation. Model from  Pegolotti Giulia PhD thesis, Universite Paris Diderot, 2014  (appendix C)
def density_matrix_photo(struc_ini, Levels_array, t_dephasing, Energies,Osc_str,SpectralEnergy,Delta_sep, n_r, n_d, n_s, m_eff, Gamma1, N, M, n, En, n_LO, I, z_elem, E_LO, T, Times_LO, Efermi, Times_other, PopEquili):

    #Dephasing time
    Tau_dephasing_matrix=np.zeros([N,N])
    for i in range(N):
        for j in range(N):
            Tau_dephasing_matrix[i,j]=(i!=j)*0.08  #it has always to be different from 0. In this particular structure it is only one between level 2 and 3.
            if (i==1 and j==2) or (i==2 and j==1): 
                Tau_dephasing_matrix[i,j]=t_dephasing
            
            
       #combinations of all the transitions            
    Transitions=list(itertools.combinations(Levels_array, 2)) 
    alfa_p=np.zeros([M],int)
    alfa_q=np.zeros([M],int)
    for p in range(M):
        alfa_p[p]=Transitions[p][0]-1  #it starts from 0 to take account how Python counts in array (from 0)
    for q in range(M):
        alfa_q[q]=Transitions[q][1]-1
    
           
      #oscillator strenghts
    f_21=np.zeros(n)
    f_31=np.zeros(n)  
    for k in range(n):
        f_21[k]=Osc_str[k,1]
        f_31[k]=Osc_str[k,2]
    
    Omega=np.zeros([n, En],dtype=object)  #All the the tunnel couplings
    Delta=np.zeros([n, En],dtype=object)  #All energies differences
    PEqu=np.zeros([n, En], dtype=object)  #Populations at thermal equilibrium (wo light)
    for k in range(n):
        for e in range(En):
            Delta[k,e]=np.zeros([N,N])
            Omega[k,e]=np.zeros([N,N])
            PEqu[k,e]=np.zeros([N,1],complex)  
            for i in range(N):
                for j in range(N):
                    PEqu[k,e][i]=PopEquili[k,i]
                    Delta[k,e][i,j]=np.abs(Energies[k,i+1]-Energies[k,j+1])
                    if (i==1 and j==2) or (i==2 and j==1):          #the tunnel couplings different from 0. To be calculated separately and added as input
                         Omega[k,e][i,j] = Delta_sep/(h_bar_eV) # Omega_ext[0]/(h_bar_eV) 
                    if (i==1 and j==3) or (i==3 and j==1):
                         Omega[k,e][i,j] = 2.483/(h_bar_eV) #Omega_ext[1]
                         

     #Depolarization shift for the doped well from Schneider book, with epsilon=epsilon_0*n_r^2
    E_dep_21=np.zeros(n)
    for k in range(n):
        E_dep_21[k]=np.abs(Energies[k,1]-Energies[k,2])+h_bar_eV*h_bar_eV*e_C*e_C*n_d*1e-18*f_21[k]/(2*n_r*n_r*eps_0*m_eff*np.abs(Energies[k,1]-Energies[k,2]))

    #Lorenztian absorption function  
    absorption=np.zeros([n, En],dtype=object)  #dimensionless
    for i in range(En) :  
        for k in range(n):     #1e3 eV to meV , 1e4 m^2 to cm^2, 1.602e-19 J to eV
            absorption[k,i]=((e_C*e_C*h_bar_SI*1e7*n_s*1e12*(f_21[k]))/(2*eps_0*c*n_r*m_eff*e_C))*Gamma1/((np.abs(Energies[k,1]-Energies[k,2])-SpectralEnergy[i])**2+Gamma1*Gamma1)

    #Absorption time
    t_abs_21=np.zeros([n, En],dtype=object)
    for i in range(En) :  
        for k in range(n):   #1e-3 meV to eV, 1.6e-19 eV to J, 1e12 s to ps
            t_abs_21[k,i]=((SpectralEnergy[i]*1e-3*e_C*1e12*n_s*1e12)/(I*absorption[k,i]))   #[ps]                

    #Single particle matrix Hamiltonian
    H_sp=np.zeros([n, En],dtype=object) 
    for k in range(n):
        for e in range(En):
            H_sp[k,e]=np.zeros([N,N],complex)
            for i in range(N):
                for j in range(N):
                    H_sp[k,e][i,j]= (i==j)*Energies[k,j+1]
                    H_sp[k,e][i,j]+= (i!=j)*0.5*h_bar_eV*Omega[k,e][i,j]

     #Dipoles, to adapt to the particular structure 
    Z=np.zeros([n,En], dtype=object) 
    for k in range(n):
        for e in range(En):
            Z[k,e]=np.zeros([N,N],complex)
            for i in range(N):
                for j in range(N):
                     Z[k,e][i,j]+=((i==0 and j==1) or (i==1 and j==0))*-np.abs(z_elem[k,1])
                     Z[k,e][i,j]+=((i==2 and j==3) or (i==3 and j==2))*-np.abs(z_elem[k,5])
                     Z[k,e][i,j]+=((i==3 and j==4) or (i==4 and j==3))*-np.abs(z_elem[k,4])
                     Z[k,e][i,j]+=((i==4 and j==5) or (i==5 and j==4))*-np.abs(z_elem[k,3])
                     Z[k,e][i,j]+=((i==5 and j==0) or (i==0 and j==5))*-np.abs(z_elem[k,2])
                     Z[k,e][i,j]+=(i==j==2)*-(struc_ini[0]/2+struc_ini[1]+struc_ini[2]/2)*1e-8  
                     Z[k,e][i,j]+=(i==j==3)*-((struc_ini[0]/2+struc_ini[1]+struc_ini[2]/2)+(struc_ini[2]/2+struc_ini[3]+struc_ini[4]/2))*1e-8
                     Z[k,e][i,j]+=(i==j==4)*-((struc_ini[0]/2+struc_ini[1]+struc_ini[2]/2)+(struc_ini[2]/2+struc_ini[3]+struc_ini[4]/2)+(struc_ini[4]/2+struc_ini[5]+struc_ini[6]/2))*1e-8
                     Z[k,e][i,j]+=(i==j==5)*-((struc_ini[0]/2+struc_ini[1]+struc_ini[2]/2)+(struc_ini[2]/2+struc_ini[3]+struc_ini[4]/2)+(struc_ini[4]/2+struc_ini[5]+struc_ini[6]/2)+(struc_ini[6]/2+struc_ini[7]+struc_ini[8]/2))*1e-8
                 

    U_total=np.zeros([n, En],dtype=object)
    U=np.zeros([n,En],dtype=object)
    for k in range(n):
        for e in range(En):
            U_total[k,e]=np.zeros([N,N],complex) 
            U_total[k,e]=(1/h_bar_eV)*(H_sp[k,e] @ Z[k,e]-Z[k,e] @ H_sp[k,e]) #commutator
            U[k,e]=np.zeros([M,1],complex)
            m=0
            for i in range(N):
                for j in range(i+1,N):
                    U[k,e][m,0]=U_total[k,e][i,j]
                    m=m+1
               
                
    C=np.zeros([n,En],dtype=object)                    
    for k in range(n):
        for e in range(En):
            C[k,e]=np.zeros([2*M,N],complex)
            for j in range(N):
                for a in range(M):
                    C[k,e][a+M,j]=(j==alfa_q[a])*2j*Omega[k,e][alfa_p[a],alfa_q[a]]
                    C[k,e][a+M,j]+=(j==alfa_p[a])*(-2)*1j*(Omega[k,e][alfa_p[a],alfa_q[a]])

    #To adapt to the particular structure for the last energies
    D=np.zeros([n, En],dtype=object)
    for k in range(n):
        for e in range(En):
            D[k,e]=np.zeros([2*M,2*M],complex)
            for a in range(M):
                for b in range(M):
                    D[k,e][a,b]=(a==b)*1j*Delta[k,e][alfa_p[a],alfa_q[a]]/h_bar_eV
                    D[k,e][4,4]=1j*(Energies[k,6]-Energies[k,7])/h_bar_eV
                    D[k,e][a,b]+=(alfa_p[a]==alfa_p[b] and alfa_q[a]!=alfa_q[b])*-1j*Omega[k,e][alfa_q[a],alfa_q[b]]
                    D[k,e][a,b]+=(alfa_p[a]==alfa_q[b] and alfa_q[a]!=alfa_p[b])*1j*Omega[k,e][alfa_q[a],alfa_p[b]]
                    D[k,e][a,b]+=(alfa_q[a]==alfa_p[b] and alfa_p[a]!=alfa_q[b])*1j*-Omega[k,e][alfa_p[a],alfa_q[b]]
                    D[k,e][a,b]+=(alfa_q[a]==alfa_q[b] and alfa_p[a]!=alfa_p[b])*1j*Omega[k,e][alfa_p[a],alfa_p[b]]
                    D[k,e][a+M,b+M]+=(a==b)*1j*Delta[k,e][alfa_p[a],alfa_q[a]]/h_bar_eV
                    D[k,e][19,19]=1j*(Energies[k,6]-Energies[k,7])/h_bar_eV
                    D[k,e][a+M,b+M]+=((alfa_p[a]==alfa_p[b]) and (alfa_q[a]!=alfa_q[b]))*1j*-Omega[k,e][alfa_q[a],alfa_q[b]]
                    D[k,e][a+M,b+M]+=(alfa_p[a]==alfa_q[b] and alfa_q[a]!=alfa_p[b])*1j*-Omega[k,e][alfa_q[a],alfa_p[b]]
                    D[k,e][a+M,b+M]+=(alfa_q[a]==alfa_p[b] and alfa_p[a]!=alfa_q[b])*1j*Omega[k,e][alfa_p[a],alfa_q[b]]
                    D[k,e][a+M,b+M]+=(alfa_q[a]==alfa_q[b] and alfa_p[a]!=alfa_p[b])*1j*Omega[k,e][alfa_p[a],alfa_p[b]]
                    D[k,e][a,b+M]+=(a==b)*1/(Tau_dephasing_matrix[alfa_p[a],alfa_q[a]])
                    D[k,e][a+M,b]+=(a==b)*1/(Tau_dephasing_matrix[alfa_p[a],alfa_q[a]])

    E=np.zeros([n, En],dtype=object)
    for k in range(n):
        for e in range(En):
            E[k,e]=np.zeros([N,M],complex)
            for i in range(N):
                for a in range(M):
                    E[k,e][i,a]=(i==alfa_p[a])*(-1j)*Omega[k,e][alfa_p[a],alfa_q[a]]
                    E[k,e][i,a]+=(i==alfa_q[a])*(1j)*Omega[k,e][alfa_p[a],alfa_q[a]]
                   

                  
#superoperator with the absorption time. To adapt where to the levels where absorption happens (0, 1 in our case)
    F_abs=np.zeros([n,En],dtype=object)
    for k in range(n):
        for e in range(En):
            F_abs[k,e]=np.zeros([N,N], complex)
            for i in range(N):
                for j in range(N):
                    F_abs[k,e][i,j]=(i==j==0)*(-1/(t_abs_21[k,e]))
                    F_abs[k,e][i,j]+=(i==j==1)*(-1/(t_abs_21[k,e]))
                    F_abs[k,e][i,j]+=(i==0 and j==1)*(+1/(t_abs_21[k,e]))  
                    F_abs[k,e][i,j]+=(i==1 and j==0)*(+1/(t_abs_21[k,e]))
                
                    
#superoperator F with all scattering rates implemented as a detailed balance. To adapt to the structure. i, j levels                
    F=np.zeros([n,En],dtype=object)
    for k in range(n):
        for e in range(En):
            F[k,e]=np.zeros([N,N], complex)
            for i in range(N):
                for j in range(N):
                    F[k,e][i,j]=(i==j==0)*(-1/(t_abs_21[k,e])-(1/Times_LO[k,2])*n_LO*(1/(np.exp((-Efermi[k]+Energies[k,6]-Energies[k,7]-E_LO)/(k_B*T))+1))-(1/Times_LO[k,1])*(n_LO))#*(1+np.exp(((-Energies[k,0]*Period*1e-2))/(k_B*T))))
                    F[k,e][i,j]+=(i==j==1)*((-1/Times_LO[k,1])*(1+n_LO)-1/Times_other[k,1]-1/(t_abs_21[k,e]))
                    F[k,e][i,j]+=(i==j==2)*((-1/Times_LO[k,5])*(1+n_LO)-1/Times_other[k,5])
                    F[k,e][i,j]+=(i==j==3)*(-(1/Times_LO[k,4])*(1+n_LO)-1/Times_other[k,4]-(1/Times_LO[k,5])*(n_LO)*(1/(np.exp((Delta[k,e][2,3]-E_LO)/(k_B*T))+1)))
                    F[k,e][i,j]+=(i==j==4)*(-(1/Times_LO[k,3])*(1+n_LO)-1/Times_other[k,3]-(1/Times_LO[k,4])*(n_LO)*(1/(np.exp((Delta[k,e][3,4]-E_LO)/(k_B*T))+1)))
                    F[k,e][i,j]+=(i==j==5)*(-(1/Times_LO[k,2])*(1+n_LO)-1/Times_other[k,2]-(1/Times_LO[k,3])*(n_LO)*(1/(np.exp((Delta[k,e][4,5]-E_LO)/(k_B*T)+1)))) 
                    F[k,e][i,j]+=(i==0 and j==1)*((1/Times_LO[k,1])*(1+n_LO)+1/Times_other[k,1]+1/(t_abs_21[k,e]))  
                    F[k,e][i,j]+=(i==1 and j==0)*((1/Times_LO[k,1])*(n_LO)+1/(t_abs_21[k,e])) 
                    F[k,e][i,j]+=((i==3 and j==2))*(1/Times_LO[k,5])*(1+n_LO)+1/Times_other[k,5]
                    F[k,e][i,j]+=((i==2 and j==3))*(1/Times_LO[k,5])*(n_LO)*(1/(np.exp((Delta[k,e][2,3]-E_LO)/(k_B*T))+1))
                    F[k,e][i,j]+=((i==4 and j==3))*(1/Times_LO[k,4])*(1+n_LO)+1/Times_other[k,4]
                    F[k,e][i,j]+=((i==3 and j==4))*(1/Times_LO[k,4])*(n_LO)*(1/(np.exp((Delta[k,e][3,4]-E_LO)/(k_B*T))+1))  
                    F[k,e][i,j]+=((i==5 and j==4))*(1/Times_LO[k,3])*(1+n_LO)+1/Times_other[k,3]
                    F[k,e][i,j]+=((i==4 and j==5))*(1/Times_LO[k,3])*(n_LO)*(1/(np.exp((Delta[k,e][4,5]-E_LO)/(k_B*T)+1)))
                    F[k,e][i,j]+=((i==0 and j==5))*((1/Times_LO[k,2])*(1+n_LO)+1/Times_other[k,2])
                    F[k,e][i,j]+=((i==5 and j==0))*((1/Times_LO[k,2])*n_LO*(1/(np.exp((-Efermi[k]+Energies[k,6]-Energies[k,7]-E_LO)/(k_B*T))+1)))        
                    
    D_inv=np.zeros([n, En],dtype=object)                
    Y=np.zeros([n, En],dtype=object)
    Y_red=np.zeros([n, En],dtype=object)
    G=np.zeros([n,En],dtype=object)
    G_n=np.zeros([n, En],dtype=object)
    Population=np.zeros([n, En],dtype=object)
    PEquilibrium=np.zeros([n, En],dtype=object)
    Coherences = np.zeros([n, En],dtype=object)
    Vector= np.zeros([n, En],dtype=object)
    j= np.zeros([n, En], dtype=object)
    for k in range(n):
        for e in range(En):                
            Vector[k,e]=np.zeros([N,1],complex)
            D_inv[k,e]=inv(D[k,e])
            Y[k,e]=D_inv[k,e] @ C[k,e]
            Y_red[k,e]=Y[k,e][0:M,:]
            G[k,e]= E[k,e] @ Y_red[k,e]-(-F[k,e])#matrix F has an inverted sign respect to Pegolotti G thesis
            G[k,e][(N-1):N,:]=1
            
            G_n[k,e]=-F_abs[k,e]
            PEquilibrium[k,e]=G_n[k,e]@ PEqu[k,e]
            PEquilibrium[k,e][(N-1):N,:]=0
            
            Population[k,e]= inv(G[k,e]) @ PEquilibrium[k,e]
            Coherences[k,e]=Y_red[k,e] @ (Population[k,e])# inv(G[k,l,e]) @ Vector[k,l,e]   
            
    #print(G[k,e])            
  #  print(PEquilibrium[k,e])
    return U, Population, Coherences
   
    