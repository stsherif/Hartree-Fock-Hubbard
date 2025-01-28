import numpy as np
import math
import matplotlib.pyplot as pylab
import sys
import copy
from bonds_n_coords import lat_bonds #triangular lattice rectangular geometery
from bonds_n_coords_diagonal import diagonal_lat_bonds #triangular lattice diagonal geometery
# This is for slave fermion mean feild Hubbard Hamiltoinian for infinte U
# Mean Field Hamiltonian is H =  
# filling = 1 means half filling
# f = fermions b = bosons

if len(sys.argv) < 5:
    print('Provide inputs: Nx Ny filling string geometery')
    sys.exit(1)
Nx = int(sys.argv[1])
Ny = int(sys.argv[2])
N = Nx * Ny
filling = float(sys.argv[3])
Np = int(filling * N)
Nh = N - Np
geometery = sys.argv[4]
t = 1
#---------------------------------------------------------------------
def print_matrix(matrix, decimal_points):
    for row in matrix:
        formating = f"{{:.{decimal_points}f}}"
        print(" ".join(formating.format(float(elem)) for elem in row))
#--------------------------------------------------------------------------
if geometery == 'diagonal':
    bonds, open_bonds = diagonal_lat_bonds(Nx,Ny)
else:    
    bonds, open_bonds = lat_bonds(Nx,Ny)
#-------------------------------------------------------------------------------------------------------------------------------------------
def fix_bonds(bonds): #this is different from my other fix_bonds function it depends on if it is (1up 2up ... 1dn 2dn ...) or (1up 1dn 2up 2dn ............) this one is for (1up 2up ......)
    fixed_bonds=copy.deepcopy(bonds)
    for i in range(len(bonds)):
        fixed_bonds.append([bonds[i][0]+N,bonds[i][1]+N])
    return fixed_bonds
fixed_bonds = fix_bonds(bonds)
print("bonds = ",bonds)
print("fixed_bonds = ",fixed_bonds)
#print(max(max(bonds)))
#------------------------------------------------------------------------
def is_hermitian(matrix):
        return np.allclose(matrix, matrix.conj().T)
#------------------------------------------------------------------------
def sfmf (bonds, bs, fs):
    fixed_bonds = fix_bonds(bonds) 
    N = max(max(bonds))
    Hamiltonianf = np.zeros((N,N))
    Hamiltonianb = np.zeros((2*N,2*N))
    for i in range(len(bonds)):
         Hamiltonianf[bonds[i][0]-1][bonds[i][1]-1] = t * bs[bonds[i][0]-1][bonds[i][1]-1]
         Hamiltonianf[bonds[i][1]-1][bonds[i][0]-1] = t * bs[bonds[i][1]-1][bonds[i][0]-1]
    for j in range(len(fixed_bonds)):
        Hamiltonianb[fixed_bonds[j][0]-1][fixed_bonds[j][1]-1] = t * fs[fixed_bonds[j][0]-1][fixed_bonds[j][1]-1]
        Hamiltonianb[fixed_bonds[j][1]-1][fixed_bonds[j][0]-1] = t * fs[fixed_bonds[j][0]-1][fixed_bonds[j][1]-1]
    #print("Hf", Hamiltonianf)
    #print("Hb", Hamiltonianb)
    eigenvaluesf, eigenvectorsf = np.linalg.eigh(Hamiltonianf)
    eigenvaluesb, eigenvectorsb = np.linalg.eigh(Hamiltonianb)
    '''
    print('vecf',eigenvectorsf)
    print('valf',eigenvaluesf)
    print('vecb',eigenvectorsb)
    print('valb',eigenvaluesb)
    '''
    #print(eigenvectors[0][1])#, "    ", eigenvectors[0,1])
    total_energy = eigenvaluesb[0] * Np
    for i in range(Nh):
        total_energy += eigenvaluesf[i]
    #------------------------------------Find the new prameters--------------------------------------------
    bs_new = np.zeros((N,N))
    fs_new = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
             for k in range(Nh):
                    fs_new[i][j] += (eigenvectorsf[i][k] * (eigenvectorsf[j][k].conj())).real
    #make fs into (2*N x 2*N)
    fs_new_extented = np.zeros((2*N,2*N))
    fs_new_extented[:N, :N] = fs_new
    fs_new_extented[N:, N:] = fs_new
    for i in range(N):
        for j in range(N):
            bs_new[i][j] = eigenvectorsb[i][0] * (eigenvectorsb[j][0]).conj() + eigenvectorsb[i+N][0] * (eigenvectorsb[j+N][0].conj())

    return total_energy, bs_new, fs_new_extented         
#-------------------------------------------------------------------------------------------------------------------------------
#make symmetric matrix
def make_sym(matrix):
    return (matrix + matrix.T) / 2
#---------------------------------------------------------------------------------
min_energy = 1000
for seed in range(2):
       np.random.seed(seed)
       #------------------------Making up the random intial obeservables--------
       bs_intial = np.random.rand(N,N)
       fs_intial = np.random.rand(N,N)
       fs_intial = make_sym(fs_intial)
       bs_intial = make_sym(bs_intial)
       '''
       sumh = 0
       sump = 0
       for i in range(N):
           sumh += fs_intial[i][i]
           sump += bs_intial[i][i]
       for i in range(N):
           for j in range(N):
               bs_intial[i][j] = bs_intial[i][j] * Np / sump
               fs_intial[i][j] = fs_intial[i][j] * Nh / sumh
       '''
       #print('fsin',fs_intial)
       #print('bsin',bs_intial)
       fs_intial_extented = np.zeros((2*N,2*N))
       fs_intial_extented[:N, :N] = fs_intial
       fs_intial_extented[N:, N:] = fs_intial
       bs = copy.deepcopy(bs_intial)
       fs = copy.deepcopy(fs_intial_extented)
       #print('fsinex',fs)
       print('new seed')
       for iteration in range(10):
            bs_old = copy.deepcopy(bs)
            fs_old = copy.deepcopy(fs)
            total_energy, bs_new, fs_new = sfmf (bonds, bs, fs)
            print("E = ", total_energy)
            if total_energy < min_energy and iteration > 5:
                min_energy = total_energy
            for i in range(N):
                for j in range(N):
                    bs[i][j] = 0.5 * (bs_old[i][j] + bs_new[i][j])
                    fs[i][j] = 0.5 * (fs_old[i][j] + fs_new[i][j])
       print("E = ", total_energy)        
       #print('E = ',total_energy)         
print('--------------------------------------------------------')
print('min energy',min_energy)
print("bs = ")
print_matrix(bs_new,8)
print("fs = ")
print_matrix(fs_new,8)
#Nb = sum(bs_new)       
#Nf = sum(fs_new)  
#print("Nup = ",Nb,'   ','Ndn = ',Nf,'   ','Np = ',Nb+Nf)

