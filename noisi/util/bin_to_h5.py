#!/Users/lermert/anaconda2/bin/python
#ToDo Docs
#ToDo build calling into main
#ToDo check whether these are all the metadata we sensibly need.


import h5py
import os
import sys
import numpy as np

#- User input: -----------------------------------------
#-------------------------------------------------------
nbytes_stationname = 512
size_of_float = 4
data_quantity = 'VEL' #'DIS','VEL','ACC'
#-------------------------------------------------------
#-------------------------------------------------------



f_in = sys.argv[1]
nbytes_total = os.path.getsize(f_in)
f_sources = sys.argv[2]
f_out_name = os.path.splitext(f_in)[0]+'.h5'

f_in = open(f_in,'rb')
f_sources = np.load(f_sources)
f_out = h5py.File(f_out_name, "w")

# Get metadata
f_in.seek(nbytes_stationname)
ntimesteps = np.fromfile(f_in,dtype='f'+str(size_of_float),count=1)[0]
ntimesteps = int(ntimesteps)
Fs = round(np.fromfile(f_in,dtype='f'+str(size_of_float),count=1)[0],6)
# Record lengths: station name plus two header values plus length of data array
nbytes_trace = nbytes_stationname + (2 + ntimesteps) * size_of_float 
# Number of records actually contained
ntraces = int(nbytes_total / nbytes_trace)
print 'This file contains %g Traces.' %ntraces
# Reference station: 
# ToDo: May differentiate location?
refstation = os.path.basename(sys.argv[1])
refstation = os.path.splitext(refstation)[0]

# DATASET NR 1: STATS
stats = f_out.create_dataset('stats',data=(0,))
stats.attrs['reference_station'] = refstation
stats.attrs['data_quantity'] = data_quantity
stats.attrs['ntraces'] = ntraces
stats.attrs['Fs'] = Fs
stats.attrs['nt'] = int(ntimesteps)

# DATASET NR 2: Source grid
sources = f_out.create_dataset('sourcegrid',data=f_sources[0:2])

# DATASET Nr 3: Seismograms itself
traces = f_out.create_dataset('data',(ntraces,ntimesteps),dtype=np.float32)


# jump to the beginning of the trace in the binary file
f_in.seek(0)
i = 0
print('Starting to read seismograms from: %s' %sys.argv[1])
while i  < ntraces:
    if i%10000 == 0:
        print('Converted %g of %g traces' %(i,ntraces))
    # read station name, copy to output file
    staname = f_in.read(nbytes_stationname)
    staname = str(staname.decode('utf-8')).strip()
    # These are only read to jump over the entries
    nt_temp = int(np.fromfile(f_in,dtype='f'+str(size_of_float),count=1)[0])
    Fs_temp = np.fromfile(f_in,dtype='f'+str(size_of_float),count=1)[0]
   
    # Get the index of that station -- > This links it with the right source coordinate pair.
    staindex = int(staname.split('.')[1])
    values = np.fromfile(f_in,dtype='f'+str(size_of_float),count=ntimesteps)
    
    # Save in traces array
    traces[staindex,:] += values
    
    i += 1
f_out.close()

