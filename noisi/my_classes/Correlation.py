#!/Users/lermert/anaconda2/bin/python
from __future__ import print_function
import numpy as np
import os
import json
import click
from glob import glob
from scipy.signal.signaltools import fftconvolve,_next_regular
from obspy import Trace, read
from noisi import NoiseSource, WaveField
from noisi.util import geo, natural_keys
from obspy.signal.invsim import cosine_taper
from noisi import filter
from scipy.signal import sosfilt
from noisi.util.windows import my_centered, zero_buddy
from noisi.util.corr_pairs import define_correlationpairs


class Correlation(object):
    
    
    def __init__(station1,station2,config,source_config):
        
        



def paths_input(cp,source_conf,step):
    
    inf1 = cp[0].split()
    inf2 = cp[1].split()
    sta1 = "{}.{}..{}".format(*inf1)
    sta2 = "{}.{}..{}".format(*inf2)
    
    # Wavefield files
    conf = json.load(open(os.path.join(source_conf['project_path'],'config.json')))
    
    wf1 = glob(os.path.join(conf['wavefield_path'],sta1+'.h5*'))[0]
    wf2 = glob(os.path.join(conf['wavefield_path'],sta2+'.h5*'))[0]
    
    # Starting model noise source
    nsrc = os.path.join(source_conf['project_path'],
                     source_conf['source_name'],'step_'+str(step),
                     'starting_model.h5')

    # Adjoint source
    adjt = os.path.join(source_conf['project_path'],
                     source_conf['source_name'],'step_'+str(step),
                     'adjt',"{}--{}.sac".format(sta1,sta2))
       
    return(wf1,wf2,nsrc,adjt)
    
def paths_output(cp,source_conf,step):
    
    inf1 = cp[0].split()
    inf2 = cp[1].split()
    sta1 = "{}.{}..{}".format(*inf1)
    sta2 = "{}.{}..{}".format(*inf2)
    
    # Correlation file
    
    corr_trace_name = "{}--{}.sac".format(sta1,sta2)    
    corr_trace_name =  os.path.join(source_conf['project_path'],
        source_conf['source_name'],'step_'+str(step),'corr',
        corr_trace_name)   

    # Kernel (without measurement) file
    kern_name = "{}--{}.npy".format(sta1,sta2)
    kern_name = os.path.join(source_conf['source_path'],
        'step_'+str(step), 'kern',
        kern_name)

    return (kern_name,corr_trace_name)
    
def g1g2_corr(wf1,wf2,corr_file,kernel,adjt,
    src,source_conf,step,autocorr=False):
    
    #ToDo: Take care of saving metainformation
    #ToDo: Think about how to manage different types of sources (numpy array vs. get from configuration -- i.e. read source from file as option)
    #ToDo: check whether to include autocorrs from user (now hardcoded off)
    #ToDo: Parallel loop(s)
    #ToDo tests

    with WaveField(wf1) as wf1, WaveField(wf2) as wf2:
        
        
        if wf1.stats['Fs'] != wf2.stats['Fs']:
            msg = 'Sampling rates of synthetic green\'s functions must match.'
            raise ValueError(msg)
        
        
        nt = int(wf1.stats['nt'])
        Fs = round(wf1.stats['Fs'],8)
        
        taper = cosine_taper(nt,p=0.05)
        n_lag = int(source_conf['max_lag'] * Fs)
        
        if nt - 2*n_lag <= 0:
            
            
            click.secho('Resetting maximum lag to %g seconds: Synthetics are too\
 short for a maximum lag of %g seconds.' %(nt//2/Fs,n_lag/Fs))
            n_lag = nt // 2
            
        n = _next_regular(nt)  
    
        # initialize new hdf5 files for correlation and green's function correlation
        #with #wf1.copy_setup(corr_file,nt=n_lag) as correl, NoiseSource(src) as nsrc:
        with NoiseSource(src) as nsrc:


            # initialize correlation and kernel
            correlation = np.zeros(n_lag)

            if step > 0:
                kern = np.zeros(wf1.stats['ntraces'])
                f = read(adjt)[0]

            # Loop over source locations
            
            for i in range(wf1.stats['ntraces']):
               
                s1 = np.ascontiguousarray(wf1.data[i,:]*taper)
                s2 = np.ascontiguousarray(wf2.data[i,:]*taper)
                
                spec1 = np.fft.rfft(s1,n)
                spec2 = np.fft.rfft(s2,n)
                
                # ToDo: Review math ask Andreas whether it is correct
                # ToDo check sign convention of correlation a(t+tau)b(t) or a(t)b(t+tau)?
                g1g2_tr = np.multiply(spec1,np.conjugate(spec2))
                #g1g2_tr = fftconvolve(s1,s2[::-1],mode='same')
                # extract Green's function correlation
                corr = my_centered(np.fft.ifftshift(np.fft.irfft(g1g2_tr,n)),n_lag)

                kern[i] = np.dot(corr,f.data) * f.stats.delta

                #correl.data[i,:] = corr.astype(np.float32)
                c = np.multiply(g1g2_tr,nsrc.get_spect(i))                
                correlation += my_centered(np.fft.ifftshift(np.fft.irfft(c,n)),n_lag)
                

            trace = Trace()
            trace.stats.sampling_rate = wf1.stats['Fs']
            trace.data = correlation
            trace.write(filename=corr_file,format='SAC')    
            
            np.save(kernel,kern)
            #correlation = correl.space_integral()
            
        
        # Save correlation
        #model_dir = os.path.join(source_conf['project_path'],\
        #source_conf['source_name'],'corr')
        #corr_file = os.path.join(model_dir, "{}--{}.sac".format(sta1,sta2))
        #correlation.write(filename=corr_file,format='SAC') #to do: Include some metadata 

# def corr(c,src,c_int,source_conf):
#     """
#     Obtain a 'noise correlation' from a correlation of Green's functions by factoring in the space-frequency dependent noise source.
    
#     :type c: String 
#     :param c: Path to the h5-file containing the correlation of Green's functions in time domain
#     :type src: String
#     :param src:  Path to the h5-file containing the time-frequency dependent source model
#     :type c_int: String
#     :param c_int: Path to the output file
#     """
    
    
    
#     with WaveField(c) as c, NoiseSource(src) as src:

#         nt = c.stats['nt']
#         #n_lag = 2 * int(source_conf['max_lag'] * c.stats['Fs']) + 1 
#         correlation = np.zeros(nt)
#         n = _next_regular(2*nt-1)
        
#         i0 = zero_buddy(nt,n,causal_function=False)
        
#         with click.progressbar(range(c.stats['ntraces']),\
#         label='Convolving source...' ) as ind:
#             for i in ind:
                
                
                
#                 #spec = np.fft.rfft(c)
#                 #c = np.multiply(g1g2_tr,nsrc.get_spect(i))                
#                 #                corr = my_centered(np.fft.ifftshift(np.fft.\
#                 #                irfft(g1g2_tr,n)),nt)
#                 #corr = fftconvolve(src.get_spect(i),c.data[i,:],mode='same')
#                 data = np.zeros(n)
#                 data[i0:i0+nt] = c.data[i,:]
#                 spec = np.fft.rfft(data)
                
#                 # also the source has to be zeropadded...it is shorter than here
                
#                 corr = np.multiply(spec,src.get_spect(i))
                
#                 correlation += my_centered(np.fft.ifftshift(np.fft.irfft(corr,n)),nt)
        
#         trace = Trace()
#         trace.stats.sampling_rate = c.stats['Fs']
#         trace.data = correlation
#         trace.write(filename=c_int,format='SAC')

# #ToDo: put somewhere else
# def kern(c,adj_src,kern,source_conf):
#     """
#     ...
#     """
    
    
#     adj_src = read(adj_src)[0]
#     delta = adj_src.stats.delta
    
#     with WaveField(c) as c:

#         kern = np.zeros(c.stats['ntraces'])
#         f = my_centered(adj_src.data,(2*source_conf['max_lag']*c.stats['Fs']+1))
        
#         kernel = np.dot(c.data[:],f) * delta
        
    
        
        
def run_corr(source_configfile,step=0):
    #ToDo think about that configuration decorator
    source_conf=json.load(open(source_configfile))
    #conf = json.load(open(os.path.join(source_conf['project_path'],'config.json')))
    
    
    p = define_correlationpairs(source_conf['project_path'])

    # for each pair:
    
    #TRY
    # get the paths to the wavefield files and the noise source file and the output (preliminary correlation and or integrated correlation)
    # is the 'preliminary run' necessary?
    # combine the preliminary correlation with the source spectrum
    #EXCEPT
    # - files not found?

    for cp in p:
        
        try:
            wf1,wf2,src,adjt = paths_input(cp,source_conf,step)
            
            kernel,corr = paths_output(cp,source_conf,step)
            
        except:
            print('Could not determine correlation for: ')
            print(cp)
            continue
            
        if not os.path.exists(corr):
            g1g2_corr(wf1,wf2,corr,kernel,adjt,src,source_conf,step)
        
        #corr(c,src,c_int)
                
            
            
        #except:
        #    print('Could not determine correlation for: ')
        ##    print(cp)
                
                #nsrc = os.path.join(source_conf['project_path'],
                #    source_conf['source_name'],
                #    'sourcemodel.h5')
                #
                
                
                
                
                
                
                