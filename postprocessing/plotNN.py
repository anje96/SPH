import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import h5py
import os
import re
import argparse
import glob


" based on https://github.com/jammartin/ParaLoBstar/blob/main/tools/conservation/main.py"
"adapted by Anne Vera Jeschke March 2023"

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Plot NN of particles.")
    parser.add_argument("--data", "-d", metavar="str", type=str, help="input directory",
                        nargs="?", default="output")
    parser.add_argument("--output", "-o", metavar="str", type=str, help="output directory",
                        nargs="?", default="output")
    parser.add_argument("--nearestneighbor", "-n", metavar="str", type=str, help="nearest neighbors",
                        nargs="?", default="5")

    args = parser.parse_args()

    time = []
    positions = []
    density = []
    nn = []

for h5file in sorted(glob.glob(os.path.join(args.data, "*.h5")), key=os.path.basename):
    print("Processing ", h5file, " ...")
    data = h5py.File(h5file, 'r')
    time.append(re.findall(r'(\d+)*.h5', h5file))

    print("...reading positions...")
    positions.append(np.array(data["r"][:]))

    print("...reading density...")
    density.append(np.array(data["rho"][:]))

    print("...reading nearest neighbors...")
    nn.append(np.array(data["NN"][:]))

print("...done.")

numPlots = len(positions)

for i in range(numPlots):
    print("...plotting timestep figure {} / {}....".format(i+1, numPlots))
    col = ["blue"]*len(positions[0])
    particle = int(args.nearestneighbor)
    col[particle] = "yellow"
    NN = nn[i]
    maxNN = len(NN[particle])
    neighbor = NN[particle,0]
    nCounter = 0
    while(neighbor != -1 and nCounter < maxNN):
        col[neighbor] = "red"
        nCounter += 1
        neighbor = NN[particle, nCounter]
    #print(col)   
    
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot() #projection='3d'
    r = positions[i]
    circle1 = plt.Circle((r[particle,0], r[particle,1]), 0.3, color="green", fill=False)
    p = ax.scatter(r[:, 0], r[:, 1], c=col, marker=".", s=2) #, r[:, 2]
    ax.add_patch(circle1)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    #ax.set_zlabel('Z')
    #to compare different timesteps
    #ax.set_xlim3d(-10, 10)
    #ax.set_ylim3d(-10, 10)
    #ax.set_zlim3d(-10, 10)
    ax.set_title('NN: {}'.format(particle))
    #fig.colorbar(p)
    #print("{0:}/NN_of_{1:04d}_at_Timestep{1:04d}.png".format(args.output, particle, i))
    #print(i)
    plt.savefig("{0:}/NN_at_Timestep{1:04d}.png".format(args.output, i))
    plt.close()