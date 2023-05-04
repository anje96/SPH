'created by Anne Vera Jeschke January 2023'
#create position, velocity and mass of particles in 3-D Grid
import numpy as np
import matplotlib.pyplot as plt
import h5py

#want an initial Distribution, set density true
density = True


#Dim
dim = 3

#Number of particles N= n^3
n = 5
N = n**dim
#coordinates of particles
r = np.zeros((N, dim))
# velocity of particles
v = np.zeros((N, dim))
# mass of particles
m = np.ones(N)

if(density):
    rho = np.ones(N)

h5f = h5py.File("cubic_N{}.h5".format(N), "w")
print("Saving to cubic_N{}.h5 ...".format(N))



#3D meshgrid
a = np.mgrid[0:n, 0:n, 0:n]
#print(a)
#write coordinates to r only for 3D
for i in range(dim):
    k=0
    l=0
    for j in range(N):
        if(j%(n**2)==0 and j>0):
            k+= 1
            k = k%(n**2)
        if(j%n==0 and j>0):
            l += 1
            l = l%n

        r[j, i] = a[i, k, l, j%n]-(n-1)/2
    
#write Initial Conditions to csv    
#np.savetxt('initialCond.csv', np.c_[r,v,m], fmt="%1.5f", delimiter=",", header="r_x, r_y, r_z, v_x, v_y, v_z, m")

#write to hdf5 data set
h5f.create_dataset("r", data=r)
h5f.create_dataset("v", data=v)
h5f.create_dataset("m", data= m)
if(density):
    h5f.create_dataset("rho", data = rho)

h5f.close()
print("Finished")