'created by Anne Vera Jeschke January 2023'
import numpy as np
import matplotlib.pyplot as plt
import h5py

"This program creates two squares in the origin which are then shifted to their wanted positions"

#Dim
dim = 2

#side length square
r_outer = 4

#shift of the rings from origin on x-axis
shift = 4.16

#particle spacing
delta_p = 0.1
#projected speed
v_p = 0.059

#want an initial Distribution, set density true
density = True
initial_density = 1 
mass = 0.01

#calc number of particles within square
N_length = int(2*r_outer/delta_p)
N_square = int(N_length**2)

#coordinates of particles in square
r = np.zeros((N_square, dim))

#2D meshgrid
a = np.mgrid[0:N_length, 0:N_length]

#create square 
for i in range(dim):
    k=0
    for j in range(N_square):
        if(j%(N_length)==0 and j>0):
            k+= 1
            k = k%(N_length)
        #print(i, k, j)

        r[j, i] = (a[i, k, j%N_length]-(N_length-1)/2)*delta_p
   
N = N_square

#construct two squares with N particles which then are shifted along the x-axis
r_ring = np.zeros((N, dim+1))
r_ring2 = np.zeros((N, dim+1))
v = np.zeros((N, dim+1))
v2 = np.zeros((N, dim+1))

m = np.ones(2*N)*mass # 2N because of two squares
rho = np.ones(2*N)*initial_density

#create square 1

for i in range(N_square):
    r_ring[i, 0] = r[i,0] - 6.0
    r_ring[i, 1] = r[i, 1]
    r_ring[i, 2] = 0.0
    v[i, 0] = v_p
    v[i, 1] = 0.0
        

#create ring 2

for i in range(N_square):
    
    r_ring2[i, 0] = r[i,0] + 6.0
    r_ring2[i, 1] = r[i, 1]
    r_ring2[i, 2] = 0.0
    v2[i, 0] = -v_p
    v[i, 1] = 0.0
    

# put two squares in one array 
print(r_ring.shape, r_ring2.shape)
r_final = np.concatenate((r_ring, r_ring2))
v_final = np.concatenate((v, v2))
print(r_final.shape)

#check symmetry
print("sum over all x-values: ", np.sum(r_final[:, 0]))
print("sum over all y-values: ", np.sum(r_final[:, 1]))
print("sum over all vx-values: ", np.sum(v_final[:, 0]))
print("sum over all vy-values: ", np.sum(v_final[:, 1]))

h5f = h5py.File("squares.h5", "w")
print("Saving to squares.h5 ...")

#write to hdf5 data set
h5f.create_dataset("r", data=r_final)
h5f.create_dataset("v", data=v_final)
h5f.create_dataset("m", data= m)
if(density):
    h5f.create_dataset("rho", data = rho)

h5f.close()
print("Finished")