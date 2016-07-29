#!/users/lermert/anaconda2/bin/python
# path to python has to be included here

# A script to filter and downsample binary specfem output using python
# This script is specially adapted to the specfem unformatted fortran binary output.


# Tasks: 
#   - Get filter coefficients
#   - read in trace
#   - apply filter to trace
#   - apply downsampling or interpolation to trace
#   - append trace to a new binary file 
#   - parallelize in a sensible way
#   - include an option to extract only the z component
#   - transfer output to hdf5

from mpi4py import MPI
from scipy.signal import iirfilter, lfilter
from obspy import Trace
import numpy as np
import os
import sys
from warnings import warn
from scipy.signal import cheb2ord, cheby2
from scipy.signal import zpk2sos, sosfilt
#ToDo: build in a more stable filter (?); cut off the first x seconds before zero time; take derivative! 
# ToDo: Write directly to hdf5? (The nice thing about unformatted bin is convenient concatenating)


#- User input: -----------------------------------------
#-------------------------------------------------------
ntimesteps=95400 # nr. of time steps, find in specfem output
duration= 300.430176 * 60. + 100. #duration of synthetics in seconds (look up in specfem output, it's given in minutes there)
offset_seconds =  100. # This is the added time specfem adds before t=0 to account for the source 'rise time'.
ncomponents = 3
nbytes_stationname = 512
size_of_float = 4
dtype_output = 'f4'
output_quantity = 'VEL'
output_directory = '/scratch/daint/lermert/output_decimated/'
channel = 'MXZ' # 'all' or specify channels as in specfem output, e.g. 'MXZ'
freq = 0.05 # Lowpass corner
fs_new = 0.4 # Interpolate to new sampling rate in Hz 

#------------------------------------------------------
#------------------------------------------------------

def cheby2_lowpass(df,freq,maxorder=8):
    # From obspy
    nyquist = df * 0.5
    # rp - maximum ripple of passband, rs - attenuation of stopband
    rp, rs, order = 1, 96, 1e99
    ws = freq / nyquist  # stop band frequency
    wp = ws  # pass band frequency
    # raise for some bad scenarios
    if ws > 1:
        ws = 1.0
        msg = "Selected corner frequency is above Nyquist. " + \
              "Setting Nyquist as high corner."
        warnings.warn(msg)
    while True:
        if order <= maxorder:
            break
        wp = wp * 0.99
        order, wn = cheb2ord(wp, ws, rp, rs, analog=0)
    return cheby2(order, rs, wn, btype='low', analog=0, output='zpk')
    
    

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print("Hello from rank ",rank)
print("Size is ",size)

try:
    f_in = sys.argv[1]
except IndexError:
    msg = 'Usage: python decimate_synthetics.py <path to input file>'
    raise ValueError(msg)
    
print('Setting output file name according to input file name provided on command line:')
station = os.path.splitext(os.path.basename(f_in))[0]
f_out  = "%s..%s.%gHz.bin_%g" %(station,channel,fs_new,rank)
dir_out = os.path.join(output_directory,'{}..{}'.format(station,channel))
f_out = os.path.join(dir_out,f_out)
if rank == 0:
    os.system('mkdir -p '+dir_out)
comm.Barrier()

print(f_out)
# Separate output files for different cores - otherwise would have to keep track of sequence of traces
f_out = open(f_out,'wb')


# Record lengths:
nbytes_trace = nbytes_stationname + 8 + ntimesteps * \
size_of_float * 2 + 16
nbytes_station = nbytes_trace * ncomponents

# Number of records actually contained
nstations = os.path.getsize(f_in) / nbytes_station
print 'Number of station records: ',nstations


# Open files:
f_in = open(f_in,'rb')

# Read an example time axis
f_in.seek(nbytes_stationname+8+4)
examp_time = np.zeros(ntimesteps)
examp_time = np.fromfile(f_in,dtype='f',count=ntimesteps)
   
print 'Prescribed duration: ', duration
print 'Inferred duration: ', examp_time[-1]-examp_time[0]
#dt = np.sum(np.abs(examp_time[:-1])) / (ntimesteps-1)
dt = duration / ntimesteps
dt_test = abs(examp_time[-1]-examp_time[0]) / ntimesteps

if dt_test != dt and rank == 0:
    msg = 'Small discrepancy between inferred and prescribed sampling rate:' 
    warn(msg)
    print dt
    print dt_test
    
# Determine the sampling rate
fs_old = 1./dt
# Get taper

# Get filter coeff
z, p, k = cheby2_lowpass(fs_old,freq)
sos = zpk2sos(z, p, k)

# Determine which channels to save...
if channel is not 'all':
    cha_incr = 3
    if channel == 'MXN':
        cha_offset = 0
    elif channel == 'MXE':
        cha_offset = 1
    elif channel == 'MXZ':
        cha_offset = 2
else:
    cha_incr = 1
    cha_offset = 0        

counter = 0
total_traces = ncomponents * nstations

for ns in range(rank,nstations,size):   
    for nc in range(cha_offset,ncomponents,cha_incr):
        
        if counter%1000 == 0: 
            print('Completed %s out of approx %s traces' %(counter,round(total_traces/size)))
        
        # jump to the beginning of the entry in the binary file
        f_in.seek(ns * nbytes_station + nc * nbytes_trace + 4)
        
        # read station name, copy to output file
        staname=f_in.read(nbytes_stationname)
        if channel != 'all' and channel not in staname.split('.'):
            msg = 'Something went wrong with reading, please double check\
             number of time steps and other user input.'
            raise ValueError(msg)
        
        
        # Numpy arrays are double precision by default, maybe use single?
        values = np.zeros(ntimesteps,dtype=np.dtype(dtype_output))
        infnr = np.zeros(2,dtype=np.dtype(dtype_output))
        
        # Jump over the station name and the time axis record....
        f_in.seek( ns * nbytes_station + nc * nbytes_trace + nbytes_stationname  +\
         8 + ntimesteps * size_of_float + 12 )
        
        #for nt in range(ntimesteps):
        #    values[nt] = np.fromfile(f_in,dtype=dtype_output,count=1)
        values = np.fromfile(f_in,dtype=dtype_output,count=ntimesteps)
        
        tr = Trace(data=values)
        # Remove the extra time that specfem added
        tr.trim(starttime = tr.stats.starttime+offset_seconds)
        
        # Filter and downsample
        # Since the same filter will be applied to all synthetics consistently, non-zero-phase should be okay
        # ToDo: Think about whether zerophase would be better
        
        # taper first
        tr.data *= taper
        tr.data = sosfilt(sos,tr.data)
        tr.stats.sampling_rate = fs_old
        tr.interpolate(fs_new)
        
        # Differentiate
        if output_quantity == 'VEL' or output_quantity == 'ACC':
            tr.differentiate()
            if output_quantity == 'ACC':
                tr.differentiate()
        
        tr.data = tr.data.astype(dtype_output)
        
        
        infnr[0] += tr.stats.npts
        infnr[1] += tr.stats.sampling_rate
            
        f_out.write(staname)
        infnr.tofile(f_out)
        tr.data.tofile(f_out)
        counter +=1
print 'New nr. of time steps after interpolation: ', tr.stats.npts
f_out.close()
        

    

