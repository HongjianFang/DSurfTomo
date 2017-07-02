#!/usr/bin/env python
# how to run:
# ./GenerateIniMOD.py
# remember to move MOD to the directory where you want to run DSurfTomo
import numpy as np

#parameters need to be changed
#start
nx=75
ny=96
nz=17
minvel=0.8
velgrad=0.5
dep1=1.5+np.array([-1.5, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,10.0,11.0,13.0,16.0,20.0,30.0])
vel1=np.loadtxt('mod.1d')
#end
vs1=np.zeros(nz)
mod=np.zeros((nz*ny,nx))
for k in range(nz):
  for j in range(ny):
    for i in range(nx):
      mod[k*ny+j,i]= vel1[k]/1.75#minvel+dep1[k]*velgrad
with open('MOD','w') as fp:
    for i in range(nz):
        fp.write('%5.1f' % dep1[i])
    fp.write('\n')
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                fp.write('%7.3f' % mod[k*ny+j,i])
            fp.write('\n')
for i in range(nz):
  print dep1[i],
