import numpy as np
import math
import matplotlib.pyplot as pylab
import sys
import copy
from bonds_n_coords import lat_bonds #triangular lattice rectangular geometery
from bonds_n_coords_diagonal import diagonal_lat_bonds #triangular lattice diagonal geometery
# This is for Hubbard Hamiltoinian
# Mean Field Hamiltonian is H= KT + U sum_i ( <niup>nidn + niup<nidn> - <niup><nidn> - <Splusi>Sminsi - Splusi<Sminusi> + <Splusi><Sminsi> )
# filling = 1 means half filling
if len(sys.argv) < 6:
    print('Provide inputs: Nx Ny filling string geometery float U')
    sys.exit(1)
Nx = int(sys.argv[1])
Ny = int(sys.argv[2])
N = Nx * Ny
filling = float(sys.argv[3])
Np = int(filling * N)
Nh = N - Np
geometery = sys.argv[4]
U = float(sys.argv[5])
t = 1

if geometery == 'diagonal':
    bonds, open_bonds = diagonal_lat_bonds(Nx,Ny)
else:    
    bonds, open_bonds = lat_bonds(Nx,Ny)
#-------------------------------------------------------------------------------------------------------------------------------------------
def fix_bonds(bonds): #this is different from my other fix_bonds function it depends on if it is (1up 2up ... 1dn 2dn ...) or (1up 1dn 2up 2dn ............) this one is for (1up 2up ......)
    for i in range(len(bonds)):
        bonds.append([bonds[i][0]+N,bonds[i][1]+N])
    print("fixed_bonds = " ,bonds)
    return bonds
bonds = fix_bonds(bonds)
#------------------------------------------------------------------------
def is_hermitian(matrix):
        return np.allclose(matrix, matrix.conj().T)
#------------------------------------------------------------------------
def hf (Hopping, nups, ndns, spluss):
    N = int(len(Hopping[0])/2)
    Hamiltonian = copy.deepcopy(Hopping)
    #Make the intercating part-------------------------
    for i in range(N):
          Hamiltonian[i][i] += U * ndns[i] #niup <nidn> 
          Hamiltonian[i+N][i+N] += U * nups[i] #<niup> nidn
          Hamiltonian[i+N][i] += - U * spluss[i] #-<Splusi> Sminsi
          Hamiltonian[i][i+N] += - U * spluss.conj()[i] #-Splusi <Sminsi>
    #print(is_hermitian(Hamiltonian))
    '''
    for row in Hamiltonian:
        formatted_row = " ".join(f"{element:3}" for element in row)
        print(formatted_row)
    '''
    eigenvalues, eigenvectors = np.linalg.eigh(Hamiltonian)
    #print(eigenvectors)
    #print(eigenvectors[0][1])#, "    ", eigenvectors[0,1])
    ef = eigenvalues[Np-1]
    total_energy = 0.0
    for i in range(Np):
        total_energy += eigenvalues[i]
    #------------------------------------Find the new prameters--------------------------------------------
    nups_new = np.zeros(N)
    ndns_new = np.zeros(N)
    spluss_new = np.zeros(N, dtype=complex)
    for i in range(N):
         for j in range(Np):
             if (eigenvectors[i][j] * (eigenvectors[i][j].conj())).imag < 0.0000001 and (eigenvectors[i+N][j] * (eigenvectors[i+N][j].conj())).imag < 0.0000001:
                    nups_new[i] += (eigenvectors[i][j] * (eigenvectors[i][j].conj())).real
                    ndns_new[i] += (eigenvectors[i+N][j] * (eigenvectors[i+N][j].conj())).real
             else:
                    print(" error : nups or udns is not real")
             spluss_new[i] += eigenvectors[i+N][j] * (eigenvectors[i][j].conj())
    for i in range(N):
        total_energy += - U * nups[i] * ndns[i]
        if (spluss[i] * spluss.conj()[i]).imag < 0.0000001:
            total_energy += U * (spluss[i] * spluss.conj()[i]).real
    return total_energy, nups_new, ndns_new, spluss_new         
#-------------------------------------------------------------------------------------------------------------------------------
Hopping = np.zeros((2*N,2*N), dtype=complex)
#Make the non inteacting part----------------------
for i in range(len(bonds)):
    Hopping[bonds[i][0]-1][bonds[i][1]-1] = -t
    Hopping[bonds[i][1]-1][bonds[i][0]-1] = -t
se=1
min_energy = 1000
for seed in range(50):
       np.random.seed(seed)
       #------------------------Making up the random intial obeservables--------
       nups_intial = np.random.rand(N)
       ndns_intial = np.random.rand(N)
       sumups = sum(nups_intial)
       sumdns = sum(ndns_intial)
       for i in range(len(nups_intial)):
            nups_intial[i] = nups_intial[i] *math.ceil((0.5) *Np) / sumups
            ndns_intial[i] = ndns_intial[i] *math.floor((0.5) *Np) / sumdns
       mags_for_spluss = np.random.rand(N)
       angles_for_spluss = np.random.rand(N) * 2 * np.pi # This is just to creat complex spluss with abs value = 1
       spluss_intial = mags_for_spluss * (np.cos(angles_for_spluss) + 1j * np.sin(angles_for_spluss))
       #print('intial up', sum(nups_intial))
       #print('intial dn', sum(ndns_intial))
       #print('intial Np', sum(nups_intial)+sum(ndns_intial))
       nups = copy.deepcopy(nups_intial)
       ndns = copy.deepcopy(ndns_intial)
       spluss = copy.deepcopy(spluss_intial)
       for iteration in range(20):
            nups_old = copy.deepcopy(nups)
            ndns_old = copy.deepcopy(ndns)
            spluss_old = copy.deepcopy(spluss)
            total_energy, nups_new, ndns_new, spluss_new = hf (Hopping, nups, ndns, spluss)
            #print("E = ", total_energy)
            if total_energy < min_energy:
                min_energy = total_energy
            for i in range(N):
                nups[i] = 0.5 * (nups_old[i] + nups_new[i])
                ndns[i] = 0.5 * (ndns_old[i] + ndns_new[i])
                spluss[i] = 0.5 * (spluss_old[i] + spluss_new[i])
       #print('E = ',total_energy)         
print('--------------------------------------------------------')
print('U = ',U)
print('min energy',min_energy)
print("nups = ",nups_new)
print("ndns = ",ndns_new)
print("spluss = ",spluss_new)
print("sminss = ",np.conj(spluss_new))
Nup = sum(nups_new)       
Ndn = sum(ndns_new)  
print("Nup = ",Nup,'   ','Ndn = ',Ndn,'   ','Np = ',Nup+Ndn)

