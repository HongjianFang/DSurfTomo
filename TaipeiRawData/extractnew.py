#!/usr/bin/env python

# a python script to extract data for the direct inversion of 
# surface wave dispersion measurements
# By Hongjian Fang @ USTC(fanghj@mail.ustc.edu.cn) 2014/04/24
# the data format is
# # sta1_lat sta1_lon nt iwave igr
#   sta2_lat sta2_lon phase/group velocity
# ......
import os
import numpy as np
'''how to use it?
python extractSurfTomo.py > surfdata.dat
and move surfdata.dat to the directory where you want to run DSurfTomo
'''

# parameters need to change according to your case
# start
nsrc=20
nrc=20
nf=26
MinP=0.5
MaxP=3.0
interval=0.1
wavetp=2
veltp=0
d = {}
d['TB01']=1
d['TB02']=2
d['TB03']=3
d['TB04']=4
d['TB05']=5
d['TB06']=6
d['TB07']=7
d['TB08']=8
d['TB09']=9
d['TB10']=10
d['TB11']=11
d['TB12']=12
d['TB13']=13
d['TB14']=14
d['TB15']=15
d['TB16']=16
d['TB17']=17
d['TB18']=18
d['TB19']=19
d['TB20']=20
#end

phav=np.zeros((nrc,nsrc,nf))
sta1la=np.zeros((nrc,nsrc,nf))
sta1lo=np.zeros((nrc,nsrc,nf))
sta2la=np.zeros((nrc,nsrc,nf))
sta2lo=np.zeros((nrc,nsrc,nf))
for file in os.listdir('.'):
  if 'CD' in file:
    data=np.loadtxt(file)
    shape=data.shape
# you may need to change the following 2 lines depend on your file name
    source=file.split('-')[0].split('.')[1]
    receiver=file.split('-')[1].split('.')[0]
    for i in range(2,shape[0]):
     if data[i,0]>=MinP and data[i,0]<=MaxP:
                  f=np.int(np.rint((data[i,0]-MinP)/interval+1))
                  rc=d[receiver]-1
                  src=d[source]-1
                  phav[rc,src,f-1]=data[i,1]
                  sta1la[rc,src,f-1]=data[0,1]
                  sta1lo[rc,src,f-1]=data[0,0]
                  sta2la[rc,src,f-1]=data[1,1]
                  sta2lo[rc,src,f-1]=data[1,0]

srclat=-999.9
per=-999
nsrcout=0

for ifr in range(nf):
 for isrc in range(nsrc):
  for irc in range(nrc):
   if phav[irc,isrc,ifr]>0:
    if np.abs(sta1la[irc,isrc,ifr]-srclat)>1e-4 or abs(ifr+1-per)>1e-4:
      nsrcout=nsrcout+1
      print ('%s %10.6f %10.6f %d %d %d') % ('#',sta1la[irc,isrc,ifr],sta1lo[irc,isrc,ifr],ifr+1,wavetp,veltp)
      print ('%10.6f %10.6f %6.4f') % (sta2la[irc,isrc,ifr],sta2lo[irc,isrc,ifr],phav[irc,isrc,ifr])
    else:
      print ('%10.6f %10.6f %6.4f') % (sta2la[irc,isrc,ifr],sta2lo[irc,isrc,ifr],phav[irc,isrc,ifr])
    srclat=sta1la[irc,isrc,ifr]
    per=ifr+1
