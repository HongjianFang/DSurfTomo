#!/usr/bin/env python
# the checkerboard should use the initial model as background model, then add 
# some pertubations 
# remember to move MOD.true to the directory where you want to run DSurfTomo
import numpy as np
import matplotlib.pyplot as plt
#parameters you need to change
#start
nx=18
ny=18
nz=9
minvel=0.9
velgrad=0.6
dep1=[0,0.2,0.4,0.6,0.8,1.1,1.4,1.8,2.5]
anosize=0.5
amplitude=0.4
#end
x=range(1,nx+1)
y=range(1,ny+1)
z=range(1,nz+1)
z=np.ones(nz)
bg=np.zeros((nz,nx,ny))
cross=np.zeros((nz,ny))
vs1=np.zeros(nz)
xy=np.kron(np.sin(anosize*np.array(y)),np.sin(anosize*np.array(x)))
xyz=np.kron(z,xy)
pxy=xyz.reshape(nz,nx,ny)
for k in range(nz):
  for j in range(ny):
    for i in range(nx):
      bg[k,i,j]= minvel+dep1[k]*velgrad
mod=np.zeros((nz*ny,nx))
for k in range(nz):
  for j in range(ny):
    for i in range(nx):
	mod[(k)*ny+j,i]=bg[k,i,j]+pxy[k,i,j]*amplitude
k=5
plt.imshow(mod[k*ny:(k+1)*ny,:],cmap='jet_r')
np.savetxt('MOD.true',mod,fmt='%4.4f')
plt.colorbar()
plt.show()

