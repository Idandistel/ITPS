#main

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
from Honeycomb2 import Honeycomb
from elasto_plasticity import elasto_plasticity
from hyper_new import hyperelasticity
#creating the structure
#parameter defenition

a=40 # face length(middle)
t=1  #  half web thickness
b=10 #isolation thickness
resolution=1

n_times_x=2 # linear pattern over the base cell in x directon
n_times_y=2 # linear pattern over the base cell in y directon

fiber_width=2
fiber_clearance=10

phi=np.zeros([2,4]) # fiber angle matrix
phi[0,0]=20
phi[0,1]=160
phi[0,2]=135
phi[0,3]=60
phi[1,0]=80
phi[1,1]=55
phi[1,2]=30
phi[1,3]=0

elastic=0
#geometry defenition
A=Honeycomb(int(a*resolution),int(t*resolution),int(b*resolution),int(fiber_width*resolution),phi,
int(fiber_clearance*resolution),n_times_x,n_times_y)

[Nx,Ny,Nz]=A.shape
print('Structure dim-'+str(A.shape))
if elastic ==0:
    F = hyperelasticity(A)
if elastic ==1:
    F=elasto_plasticity(A)


#material structure view
fig = plt.figure(5)
ax = fig.add_subplot(111, projection='3d')


x,y,z=np.mgrid[:Nx,:Ny,:Nz]

c=F           # 3D deformation gradient tensor
img = ax.scatter(x, y, z, c=c, cmap=plt.hot())
fig.colorbar(img)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()