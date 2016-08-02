from mpi4py import MPI
import os
import numpy as np
import json
from obspy import read
from noisi import WaveField
from noisi.util.windows import my_centered, zero_buddy
from scipy.signal.signaltools import _next_regular
import matplotlib.pyplot as plt


def kern_td(c,adj_src,kern):
    """
    ...
    """
    
    
    adj_src = read(adj_src)[0]
    delta = adj_src.stats.delta
    
    with WaveField(c) as c:

        ntraces = c.stats['ntraces']
        nt = c.stats['nt']
        
        f = my_centered(adj_src.data,nt)
        kernel = np.zeros(ntraces,dtype=np.float32)
        # Very basic frequency integration
        # This vectorized form does seem to be impossible memory-wise
        # kernel = np.dot(c.data[:],f) * delta
        for i in range(ntraces):
        
            kernel[i] += np.dot(c.data[i,:],f)
        
        kernel *= delta
        # Save as npy
        np.save(kern,kernel)
 
def kern_fd(c,adj_src,kern,n):
    """
    ...
    """
    adj_src = read(adj_src)[0]
    
    
    
    with WaveField(c,complex=True) as c:
        ntraces = c.stats['ntraces']
        nf = c.stats['nt']
        delta = adj_src.stats.sampling_rate
    
        k = np.zeros(ntraces,dtype=np.float32)
    
    
        # Zero padding of adjoint source.
        # Problem: If the data / the adjoint source are shorter than the generating wavefield, then how to reconstruct the right frequency axis here???
        
        # n is passed in as argument now n = _next_regular(2*adj_src.stats.npts-1)
        i0 = zero_buddy(adj_src.stats.npts,n,causal_function=False)
        adj = np.zeros(n)
        adj[i0:i0+adj_src.stats.npts] += adj_src.data
        
        
        # FFT
        f = np.fft.rfft(adj)
        plt.plot(np.fft.rfftfreq(n,d=delta),np.imag(f))
        plt.show()
        # multiply and rudimentary integral
        for i in range(ntraces):

            k[i] = np.real(np.dot( (1j*c.data_i[i,:] + c.data_r[i,:]), f))
        
        k *= delta
        
        np.save(kern,k)
    
    
     
 
#def get_kernel_function(ktype):
#    
#    if ktype == 'td':
#        func = kern_td
#        print(func)
#    elif ktype == 'fd':
#        func = kern_fd
#        print(func)
#    else:
#        msg = ('Kernel type {} is not implemented'.format(ktype))
#        raise NotImplementedError(msg)
#    
#    return func
#               
    
def run_kernel(source_configfile,step):
    
    source_config=json.load(open(source_configfile))
    
    step_n = 'step_{}'.format(int(step))
    step_dir = os.path.join(source_config['source_path'],step_n)
    ktype = source_config['ktype']
    # All the adjoint sources available:
    files = [f for f in os.listdir(os.path.join(step_dir,'adjt'))]
    files = [os.path.join(step_dir,'adjt',f) for f in files]
    
    if files == []:
        print("No adjoint sources found.")
        return()


    # very simple embarrassingly parallel run
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    files = files[rank:len(files):size]

     
    #kernel_function = get_kernel_function(source_config['ktype'])
       
    for f in files:
        
        corr = os.path.splitext(os.path.basename(f))[0]
        sta1 = corr.split('--')[0]
        sta2 = corr.split('--')[1]
        
        if ktype == 'fd':
            kfile = os.path.join(source_config['source_path'],'green_c_fd',
            '{}--{}.h5'.format(sta1,sta2))      
        else:                                   
            kfile = os.path.join(source_config['source_path'],'green_c',
            '{}--{}.h5'.format(sta1,sta2))
        print(kfile)
        
        koutfile = os.path.join(step_dir,'grad',
        '{}--{}.npy'.format(sta1,sta2))
        print koutfile
        
        adjfile = os.path.join(step_dir,'adjt',
        '{}--{}.sac'.format(sta1,sta2))
        print koutfile
        
        if ktype == 'td':
            kern_td(kfile,adjfile,koutfile)
        
        elif ktype == 'fd':    
            fh = open(os.path.join(source_config['project_path'],
            source_config['source_name'],'n.txt'))
            n = int(fh.read())
            print(n)
            kern_fd(kfile,adjfile,koutfile,n)
            
        else:
            msg = ('Kernel type {} is not implemented'.format(ktype))
            raise NotImplementedError(msg)