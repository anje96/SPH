'created by Anne Vera Jeschke February 2023'
import numpy as np
import matplotlib.pyplot as plt
import h5py

"this program creates two 2-D rings around the origin which are then shifted to their final position"
#Dim of Rings
dim = 2

#ring properties: inner and outer radius
r_inner = 3
r_outer = 4

#shift of the rings from origin on x-axis
shift = 4.16

#particle spacing
delta_p = 0.1
#projected speed
v_p = 0.059

#want an initial Distribution, set density true
density = True
initial_density = 1 #mass of particles is set to 0.01, particle spacing is 0.1, dim = 2, 
mass = 0.01

#create initial distribution through creating a 2d grid with dims 2*r_outer x 2*r_outer, then delete particles which are not on the ring

#calc number of particles within square
N_length = int(2*r_outer/delta_p)
#number of particles in square
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
    


# count particles in one ring
N= 0
arr = np.zeros(N_square) 
for i in range(N_square):
    radius = np.sqrt(r[i,0]**2 + r[i,1]**2)
    if( radius >= r_inner and radius <= r_outer):
        N += 1
        arr[i] = 1

#construct two rings with N particles which then are shifted along the x-axis
r_ring = np.zeros((N, dim+1)) # first ring
r_ring2 = np.zeros((N, dim+1)) # second ring
v = np.zeros((N, dim+1))
v2 = np.zeros((N, dim+1))

m = np.ones(2*N)*mass # 2N because of two rings
rho = np.ones(2*N)*initial_density

#create ring 1
counter = 0
for i in range(N_square):
    if(arr[i]== 1):
        r_ring[counter, 0] = r[i,0] + 12.9
        r_ring[counter, 1] = r[i, 1] + 8
        r_ring[counter, 2] = 0
        v[counter, 0] = -v_p
        counter += 1

#create ring 2
counter = 0
for i in range(N_square):
    if(arr[i]== 1):
        r_ring2[counter, 0] = r[i,0] + 4.1
        r_ring2[counter, 1] = r[i, 1] + 8
        r_ring2[counter, 2] = 0
        v2[counter, 0] = v_p
        counter += 1


# put two rings in one array 
print(r_ring.shape, r_ring2.shape)
r_final = np.concatenate((r_ring, r_ring2))
v_final = np.concatenate((v, v2))
print(r_final.shape)

h5f = h5py.File("rings.h5", "w")
print("Saving to rings.h5 ...")

#write to hdf5 data set
h5f.create_dataset("r", data=r_final)
h5f.create_dataset("v", data=v_final)
h5f.create_dataset("m", data= m)
if(density):
    h5f.create_dataset("rho", data = rho)

h5f.close()
print("Finished")