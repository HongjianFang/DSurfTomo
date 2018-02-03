#!/usr/bin/env python

import numpy as np
from matplotlib import pylab as plt

#parameters you need to change
#begin
numray = 3320 # how many rays do you have (all periods)
rb = 1 # ray index for some period (first one)
re = 50 # ray index for some period (last one)
maxseg = 500 # 
#end

n = 0
numseg = np.zeros(numray,)
raylat = np.zeros((numray,maxseg))
raylon = np.zeros((numray,maxseg))
with open('raypath.out') as ray:
    for line in ray:
        linseg = line.split()
        if linseg[0]=='#':
            n = n+1
            seg = 0
            numseg[n-1] = float(linseg[1])
        else:
            seg = seg+1
            raylat[n-1,seg-1] = float(linseg[0])
            raylon[n-1,seg-1] = float(linseg[1])

for ii in range(rb,re):
    plt.plot(raylon[ii,0:numseg[ii]],raylat[ii,0:numseg[ii]],'k-')

plt.show()
