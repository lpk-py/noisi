
# coding: utf-8

from obspy import read
from obspy.geodetics import gps2dist_azimuth
import numpy as np
from glob import glob
import os

#######################
# USER INPUT
#######################
synthetics_dir = '.'
data_dir = 'data/'
# coordinates:
coords = {'EQ.000..MXZ':(0.00,0.00),'EQ.003..MXZ':(0.00,45.),'EQ.004..MXZ':(0.00,67.5)}

#######################


traces = glob(os.path.join(synthetics_dir,'*.sac'))
amps = np.random.rand(len(traces)*2)


cnt = 0
for t in traces:
    sta1 = os.path.splitext(os.path.basename(t))[0].split('--')[0]
    sta2 = os.path.splitext(os.path.basename(t))[0].split('--')[1]
    
    datafile = os.path.join(data_dir,os.path.basename(t))
    tr = read(t)[0]
    
    i = tr.stats.npts // 2 + 1 
    tr.data[0:i] *= amps[cnt]
    cnt += 1
    tr.data[i+1:] *= amps[cnt]
    cnt += 1
    
    tr.stats.network = sta1.split('.')[0]
    tr.stats.station = sta1.split('.')[1]
    tr.stats.location = ''
    tr.stats.channel = sta1.split('.')[3]
    tr.stats.sac={}
    tr.stats.sac.kuser0 = sta2.split('.')[0]
    tr.stats.sac.kevnm = sta2.split('.')[1]
    tr.stats.sac.kuser1 = ''
    tr.stats.sac.kuser2 = sta2.split('.')[3]
    tr.stats.sac.stla = coords[sta1][0]
    tr.stats.sac.stlo = coords[sta1][1]
    tr.stats.sac.evla = coords[sta2][0]
    tr.stats.sac.evlo = coords[sta2][1]
    tr.stats.sac.dist = gps2dist_azimuth(coords[sta1][0],coords[sta1][1],coords[sta2][0],coords[sta2][1])[0]
    tr.write(datafile,format='SAC')
   



