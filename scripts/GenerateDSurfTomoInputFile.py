#/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
#
# File Name : GenerateDSurfTomoInputFile.py
#
# Purpose : Generate the Input File for DSurfTomo using the data themselves
#
# Creation Date : 05-07-2023
#
# Last Modified : Tue 18 Jul 2023 08:10:22 PM CST
#
# Created By : Hongjian Fang: fanghj1990@gmail.com
#
#_._._._._._._._._._._._._._._._._._._._._.*/

import os
import numpy as np
import glob
import pandas as pd

#These are for generating the right format
DataPath = 'TaipeiRawData'
DataNames = 'CDisp*.dat'
wavetype = 2
phasetype = 0

# These following can be changed in DSurfTomo.in
nz = 9
sublayers = 3
minvel = 0.5
maxvel = 2.8
noiselevel = 0.02
datathreshold = 1.5
DataFiles = glob.glob(DataPath+'/'+DataNames)

sta1latall = []
sta1lonall = []
sta2latall = []
sta2lonall = []
periodsall = []
dispersionall = []
stationpairID = []

for ifile in DataFiles:
    dispersion = np.loadtxt(ifile)
    pairID = ifile.split('/')[-1]
    sta1lon = dispersion[0,0]
    sta1lat = dispersion[0,1]
    sta2lon = dispersion[1,0]
    sta2lat = dispersion[1,1]
    nperiods,_ = dispersion.shape
    nperiods = nperiods - 2
    periods = np.zeros(nperiods,)
    disper = np.zeros(nperiods,)
    for ii in range(nperiods):
        periods[ii] = dispersion[ii+2,0]
        disper[ii] = dispersion[ii+2,1]
   
    dispersionall = np.hstack([dispersionall,disper])
    periodsall = np.hstack([periodsall,periods])
    sta1latall = np.hstack([sta1latall,np.ones(nperiods,)*sta1lat])
    sta1lonall = np.hstack([sta1lonall,np.ones(nperiods,)*sta1lon])
    sta2latall = np.hstack([sta2latall,np.ones(nperiods,)*sta2lat])
    sta2lonall = np.hstack([sta2lonall,np.ones(nperiods,)*sta2lon])
    stationpairID = stationpairID + [pairID] * nperiods

dataall = pd.DataFrame({'sta1lat':sta1latall, 'sta1lon':sta1lonall, \
                        'sta2lat':sta2latall, 'sta2lon':sta2lonall, \
                        'periods':periodsall, 'dispersion': dispersionall, \
                        'pairid':stationpairID})

dataall = dataall.sort_values(by = ['periods', 'pairid'])
dataall = dataall.set_index('periods')
with open('surfdata.dat','w') as fout:
    UniqPeriods = dataall.index.unique()
    for iperiod,period in enumerate(UniqPeriods):
        datasubperiod = dataall.loc[period]
        if isinstance(datasubperiod,pd.Series):
            continue
        datasubperiod = datasubperiod.reset_index().set_index('sta1lat')
        sta1lat = datasubperiod.index.unique()
        for ista1 in sta1lat:
            datasubstation = datasubperiod.loc[ista1]
            if isinstance(datasubstation,pd.DataFrame):
                fout.write(f'# {datasubstation.index[0]} {datasubstation["sta1lon"].iloc[0]} {iperiod+1} {wavetype} {phasetype}\n')
                for ista2 in range(len(datasubstation)):
                    fout.write(f'{datasubstation["sta2lat"].iloc[ista2]} {datasubstation["sta2lon"].iloc[ista2]} {datasubstation["dispersion"].iloc[ista2]}\n')
            else:
                fout.write(f'# {datasubstation.name} {datasubstation["sta1lon"]} {iperiod+1} {wavetype} {phasetype}\n')
                fout.write(f'{datasubstation["sta2lat"]} {datasubstation["sta2lon"]} {datasubstation["dispersion"]}\n')


print('Finish reformatting disperison data')

# Determin the grid interval, 1/3 of the smallest station distance

# Weight and Damp, you will have to run the code a couple of time to decide
# The default values is 5.0 and 1.0, respectively
weight = 1.0
damp = 1.0

# Sparsity fraction, this parameter represent the sparsity of the sensitivity matrix and is set to be 0.2 for safe keeping, but could be as low as 0.02.it is not important as long as your code does not throw you memory erros
sparsityfraction = 0.1

# Iteration number, the default is 10 time, most often you get stable results after 4-5 iterations.
maxiteration = 10


distall = np.sqrt((dataall['sta2lat']-dataall['sta1lat'])**2+(dataall['sta2lon']-dataall['sta1lon'])**2)
mindist = distall.min()
gridintval = int(np.ceil(mindist/3*1000))/1000
# Determine the origin, 1 grids
originLat = np.max([np.max(dataall.sta1lat.max()),np.max(dataall.sta2lat.max())]) + 1*gridintval
originLon = np.min([np.min(dataall.sta1lon.min()),np.min(dataall.sta2lon.min())]) - 1*gridintval
largestLat = np.max([dataall['sta1lat'].max()-dataall['sta1lat'].min(),dataall['sta2lat'].max()-dataall['sta2lat'].min()])
largestLon = np.max([dataall['sta1lon'].max()-dataall['sta1lon'].min(),dataall['sta2lon'].max()-dataall['sta2lon'].min()])
inputfile = open('DSurfTomo.in','w')
dummy = 'cccc'
inputfile.write(f'{dummy*10}\n')
inputfile.write(f'{dummy*10}\n')
inputfile.write(f'{dummy*10}\n')
inputfile.write('surfdata.dat\n')
nx = int(np.ceil(largestLat/gridintval))+5
ny = int(np.ceil(largestLon/gridintval))+5
nreceivers = len(dataall.sta1lat.unique())+1
inputfile.write(f'{nx} {ny} {nz}\n')
inputfile.write(f'{originLat:<7.3f} {originLon:7.3f}\n')
inputfile.write(f'{gridintval:<7.3f} {gridintval:7.3f}\n')
inputfile.write(f'{nreceivers}\n')
inputfile.write(f'{weight} {damp}\n')
inputfile.write(f'{sublayers}\n')
inputfile.write(f'{minvel} {maxvel}\n')
inputfile.write(f'{maxiteration}\n')
inputfile.write(f'{sparsityfraction}\n')
nperiods = len(UniqPeriods)
inputfile.write(f'{nperiods}\n')
#inputfile.write(f'{*UniqPeriods}\n')
inputfile.write(' '.join([str(iperiod) for iperiod in UniqPeriods])+'\n')
inputfile.write(f'0\n')
inputfile.write(f'0\n')
inputfile.write(f'0\n')
inputfile.write(f'0\n')
inputfile.write(f'{noiselevel}\n')
inputfile.write(f'{datathreshold}\n')
inputfile.close()

print('Finishing generating input file\n')
