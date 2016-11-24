#!/Users/lermert/anaconda2/bin/python
from __future__ import print_function
from mpi4py import MPI
import numpy as np
import os
import h5py
import json
import click
from glob import glob
from scipy.signal.signaltools import fftconvolve
from scipy.fftpack import next_fast_len
from obspy import Trace, read
from noisi import NoiseSource, WaveField
from noisi.util import geo, natural_keys
from obspy.signal.invsim import cosine_taper
from noisi import filter
from scipy.signal import sosfilt
from noisi.util.windows import my_centered, zero_buddy
from noisi.util.corr_pairs import define_correlationpairs, rem_fin_prs, rem_no_obs
import matplotlib.pyplot as plt
#ToDo: put in the possibility to run on mixed channel pairs
def paths_input(cp,source_conf,step,kernelrun):
    
    inf1 = cp[0].split()
    inf2 = cp[1].split()
    
    
    # Wavefield files
    conf = json.load(open(os.path.join(source_conf['project_path'],'config.json')))
    channel = source_conf['channel']
    sta1 = "{}.{}..{}".format(*(inf1[0:2]+[channel]))
    sta2 = "{}.{}..{}".format(*(inf2[0:2]+[channel]))
    
    if source_conf['preprocess_do']:
        dir = os.path.join(source_conf['source_path'],'wavefield_processed')
    else:
        dir = conf['wavefield_path']
        
   

   
    wf1 = glob(os.path.join(dir,sta1+'.h5'))[0]
    wf2 = glob(os.path.join(dir,sta2+'.h5'))[0]

    
    
    # Starting model for the noise source
    if kernelrun:
        # this is a bit of an odd construction. The base model contains no spatial weights if the gradient is for spatial weights.
        # If even the spectral weights get updated, then it should not even contain spectral weights. 
        # But so far, so good.
        nsrc = os.path.join(source_conf['project_path'],
                     source_conf['source_name'],'step_'+str(step),
                     'base_model.h5')
    else:
        nsrc = os.path.join(source_conf['project_path'],
                     source_conf['source_name'],'step_'+str(step),
                     'starting_model.h5')

    # Adjoint source
    

    try:
        adjt = os.path.join(source_conf['source_path'],
                     'step_'+str(step),
                     'adjt',"{}--{}.sac".format(sta1,sta2))
        adjt = glob(adjt)[0]
    except IndexError:
        print("No adjoint source found for station pair: {}, {}".format(sta1,sta2))
            # ToDo: this is too horrible, please find another solution.
        adjt = '-'

    
    return(wf1,wf2,nsrc,adjt)
    
    
def paths_output(cp,source_conf,step):
    

    id1 = cp[0].split()[0]+cp[0].split()[1]
    id2 = cp[1].split()[0]+cp[1].split()[1]

    if id1 < id2 :
        inf1 = cp[0].split()
        inf2 = cp[1].split()
    else:
        inf2 = cp[0].split()
        inf1 = cp[1].split()

    channel = source_conf['channel']
    sta1 = "{}.{}..{}".format(*(inf1[0:2]+[channel]))
    sta2 = "{}.{}..{}".format(*(inf2[0:2]+[channel]))
    
    # Correlation file
    #corr_name = "{}--{}.h5".format(sta1,sta2)
     # Kernel (without measurement) file
    kern_name = "{}--{}.npy".format(sta1,sta2)
    kern_name = os.path.join(source_conf['source_path'],
        'step_'+str(step), 'kern',
        kern_name)

   
    
    # ktype = source_conf['ktype']
    # print(ktype)
    # if ktype == 'fd':
    #     corr_name = os.path.join(source_conf['source_path'],
    #     'green_c_fd',
    #     corr_name)
    # elif ktype == 'td':
    #     corr_name = os.path.join(source_conf['source_path'],
    #     'green_c',
    #     corr_name)
    # else:
    #     msg = ('ktype parameter must be fd or td.')
    #     raise NotImplementedError(msg)
        
    corr_trace_name = "{}--{}.sac".format(sta1,sta2)    
    corr_trace_name =  os.path.join(source_conf['source_path'],
        'step_'+str(step),'corr',
        corr_trace_name)   
    return (kern_name,corr_trace_name)
    
def get_ns(wf1,source_conf):
    
    # Nr of time steps in traces
    with WaveField(wf1) as wf1:
        nt = int(wf1.stats['nt'])
        Fs = round(wf1.stats['Fs'],8)

    print(Fs)
    
    # Necessary length of zero padding for carrying out frequency domain correlations/convolutions
    n = next_fast_len(2*nt-1)     
    
    # Number of time steps for synthetic correlation
    n_lag = int(source_conf['max_lag'] * Fs)
    if nt - 2*n_lag <= 0:
        click.secho('Resetting maximum lag to %g seconds: Synthetics are too\
 short for a maximum lag of %g seconds.' %(nt//2/Fs,n_lag/Fs))
        n_lag = nt // 2
        
    n_corr = 2*n_lag + 1
    
    return nt,n,n_corr
        
    


def g1g2_corr(wf1,wf2,corr_file,kernel,adjt,
    src,source_conf,kernelrun):
    
    #ToDo: Take care of saving metainformation
    #ToDo: Think about how to manage different types of sources (numpy array vs. get from configuration -- i.e. read source from file as option)
    #ToDo: check whether to include autocorrs from user (now hardcoded off)
    #ToDo: Parallel loop(s)
    #ToDo tests
    
    ntime, n, n_corr = get_ns(wf1,source_conf)
    
    
    taper = cosine_taper(ntime,p=0.05)
    print(kernelrun)
    with WaveField(wf1) as wf1, WaveField(wf2) as wf2:
        
        
        if wf1.stats['Fs'] != wf2.stats['Fs']:
            msg = 'Sampling rates of synthetic green\'s functions must match.'
            raise ValueError(msg)
        
        # initialize new hdf5 files for correlation and green's function correlation
        #with wf1.copy_setup(corr_file,nt=n_corr) as correl, NoiseSource(src) as nsrc:
        #with wf1.copy_setup(corr_file,nt=n_corr) as correl:
        

        with NoiseSource(src) as nsrc:

            correlation = np.zeros(n_corr)
            
            if kernelrun:
                
                #if not os.path.exists(adjt):
                #    print('Adjoint source %s not found, skipping kernel.')
                #    return()

                kern = np.zeros(wf1.stats['ntraces'])
                f = read(adjt)[0]
                f.data = my_centered(f.data,n_corr)
                
            # Loop over source locations
            #with click.progressbar(range(wf1.stats['ntraces']),\
            #label='Correlating...' ) as ind:
            for i in range(wf1.stats['ntraces']):
               
                s1 = np.ascontiguousarray(wf1.data[i,:]*taper)
                s2 = np.ascontiguousarray(wf2.data[i,:]*taper)
                
                spec1 = np.fft.rfft(s1,n)
                spec2 = np.fft.rfft(s2,n)
                
              
                g1g2_tr = np.multiply(np.conjugate(spec1),spec2)
                c = np.multiply(g1g2_tr,nsrc.get_spect(i))


                if kernelrun:
                    # The frequency spectrum of the noise source is included here
                    corr_temp = my_centered(np.fft.ifftshift(np.fft.irfft(c,n)),n_corr)
                    # A Riemann sum -- one could actually build in a more fancy integration here
                    kern[i] = np.dot(corr_temp,f.data) * f.stats.delta
                    
                
                else:
                                    
                    correlation += my_centered(np.fft.ifftshift(np.fft.irfft(c,n)),n_corr)
                
                if i%50000 == 0:
                    print("Finished {} source locations.".format(i))
        

        

        if kernelrun:
            np.save(kernel,kern) 

        else:
            trace = Trace()
            trace.stats.sampling_rate = wf1.stats['Fs']
            trace.data = correlation
            trace.write(filename=corr_file,format='SAC')
            
        

def g1g2_kern(wf1,wf2,corr_file,kernel,adjt,
    src,source_conf):
    
    #ToDo: Take care of saving metainformation
    #ToDo: Think about how to manage different types of sources (numpy array vs. get from configuration -- i.e. read source from file as option)
    #ToDo: check whether to include autocorrs from user (now hardcoded off)
    #ToDo: Parallel loop(s)
    #ToDo tests
    
    ntime, n, n_corr = get_ns(wf1,source_conf)
    
    
    taper = cosine_taper(ntime,p=0.05)


    with WaveField(wf1) as wf1, WaveField(wf2) as wf2:
        
        
        if wf1.stats['Fs'] != wf2.stats['Fs']:
            msg = 'Sampling rates of synthetic green\'s functions must match.'
            raise ValueError(msg)
        
        # initialize new hdf5 files for correlation and green's function correlation
        #with wf1.copy_setup(corr_file,nt=n_corr) as correl, NoiseSource(src) as nsrc:
        #with wf1.copy_setup(corr_file,nt=n_corr) as correl:
        

        with NoiseSource(src) as nsrc:

            correlation = np.zeros(n_corr)
            

            kern = np.zeros(wf1.stats['ntraces'])

            # Try to use a trick: Use causal and acausal part of f separately.
            f = read(adjt)[0]
            f.data = my_centered(f.data,n_corr)
            
            n_acausal_samples = (f.stats.npts-1)/2
            specf = np.fft.rfft(f[n_acausal_samples:],n)
            # Loop over source locations
            #with click.progressbar(range(wf1.stats['ntraces']),\
            #label='Correlating...' ) as ind:
            for i in range(wf1.stats['ntraces']):

               
                s1 = np.ascontiguousarray(wf1.data[i,:]*taper)
                spec1 = np.fft.rfft(s1,n)
                T = np.multiply(np.conjugate(spec1),nsrc.get_spect(i))
                # plt.plot(np.abs(spec1)/np.max(np.abs(spec1)))
                # plt.plot(np.abs(T)/np.max(np.abs(T)))
                # plt.show()
                T = np.fft.irfft(T,n)[-len(s1):]
                # if i in [1,2,3,4,5]:
                #     plt.plot(T[::-1]/np.max(np.abs(T)))
                #     plt.plot(s1/np.max(np.abs(s1)),'--')
                #     plt.show()                
                # Get s2 in the shape of f
                # we need to add half of the length of f before G to replace the
                # acausal part that G does not have to start with:
                
                #s2 = np.zeros(n)
                s2 = np.ascontiguousarray(wf2.data[i,:]*taper)
                #s2[n_acausal_samples:n_acausal_samples+ntime] = s2_caus
                spec2 = np.fft.rfft(s2,n)
                # plt.plot(s2_caus/np.max(np.abs(s2_caus)))
                # plt.plot(f.data/np.max(np.abs(f.data)))
                # plt.plot(s2/np.max(np.abs(s2)))
                # plt.plot(n_acausal_samples,0.5,'rd')

                # plt.show()

                # transform both f and s2 to fourier d 
                # (again, zeropadding but just to avoid circular convolution)
                
                
                
                # plt.plot(np.abs(spec2)/np.max(np.abs(spec2)))
                # plt.plot(np.abs(specf)/np.max(np.abs(specf)))
                # plt.show()
                g2f_tr = np.multiply(np.conjugate(spec2),specf)
                #plt.plot(n_acausal_samples,0.5,'rd')
                #plt.plot(n,0.5,'gd')
                u_dagger = np.fft.irfft(g2f_tr,n)[-len(s1):]
                #plt.plot(u_dagger/np.max(np.abs(u_dagger)))
                #plt.plot(T[::-1]/np.max(np.abs(T)))

                #plt.show()

            
                # The frequency spectrum of the noise source is included here
               ## corr_temp = my_centered(np.fft.ifftshift(np.fft.irfft(c,n)),n_corr)
                # A Riemann sum -- one could actually build in a more fancy integration here
                kern[i] = np.dot(u_dagger,T) * f.stats.delta
                

                if i%50000 == 0:
                    print("Finished {} source locations.".format(i))

                #np.save(kernel,kern) 
            return(kern)

            

#def corr(c,src,c_int,source_conf):
#    """
#    Obtain a 'noise correlation' from a correlation of Green's functions by factoring in the space-frequency dependent noise source.
#    
#    :type c: String 
#    :param c: Path to the h5-file containing the correlation of Green's functions in time domain
#    :type src: String
#    :param src:  Path to the h5-file containing the time-frequency dependent source model
#    :type c_int: String
#    :param c_int: Path to the output file
#    """
#    
#    
#    
#    with WaveField(c) as c, NoiseSource(src) as src:
#
#        nt = c.stats['nt']
#        #n_lag = 2 * int(source_conf['max_lag'] * c.stats['Fs']) + 1 
#        correlation = np.zeros(nt)
#        n = _next_regular(2*nt-1)
#        
#        i0 = zero_buddy(nt,n,causal_function=False)
#        
#        with click.progressbar(range(c.stats['ntraces']),\
#        label='Convolving source...' ) as ind:
#            for i in ind:
#                
#                
#                
#                #spec = np.fft.rfft(c)
#                #c = np.multiply(g1g2_tr,nsrc.get_spect(i))                
#                #                corr = my_centered(np.fft.ifftshift(np.fft.\
#                #                irfft(g1g2_tr,n)),nt)
#                #corr = fftconvolve(src.get_spect(i),c.data[i,:],mode='same')
#                data = np.zeros(n)
#                data[i0:i0+nt] = c.data[i,:]
#                spec = np.fft.rfft(data)
#                
#                # also the source has to be zeropadded...it is shorter than here
#                
#                corr = np.multiply(spec,src.get_spect(i))
#                
#                correlation += my_centered(np.fft.ifftshift(np.fft.irfft(corr,n)),nt)
#        
#        trace = Trace()
#        trace.stats.sampling_rate = c.stats['Fs']
#        trace.data = correlation
#        trace.write(filename=c_int,format='SAC')




def run_corr(source_configfile,step,kernelrun=False,steplengthrun=False,obs_only=True):

    step = int(step)


    #ToDo think about that configuration decorator
    source_config=json.load(open(source_configfile))
    #conf = json.load(open(os.path.join(source_conf['project_path'],'config.json')))
    
    p = define_correlationpairs(source_config['project_path'])
    
    # Remove pairs for which no observation is available
    if obs_only:
        directory = os.path.join(source_config['source_path'],'observed_correlations')
        p = rem_no_obs(p,source_config,directory=directory)
    if steplengthrun:
        directory = os.path.join(source_config['source_path'],
            'step_'+str(step),'obs_slt')
        p = rem_no_obs(p,source_config,directory=directory)
        

    # Remove pairs that have already been calculated
    p = rem_fin_prs(p,source_config,step,kernelrun)

    # for each pair:
    
    #TRY
    # get the paths to the wavefield files and the noise source file and the output (preliminary correlation and or integrated correlation)
    # is the 'preliminary run' necessary?
    # combine the preliminary correlation with the source spectrum
    #EXCEPT
    # - files not found?


    # simple embarrassingly parallel run:

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    

    # The assignment of station pairs should be such that one core has as many occurrences of the same station as possible; 
    # this will prevent that many processes try to access the same hdf5 file all at once.
    num_pairs = int( round(float(len(p))/float(size)) )
    p_p = p[ rank*num_pairs : rank*num_pairs + num_pairs] 
    
    print('Rank number %g' %rank)
    print('working on pair nr. %g to %g of %g.' %(rank*num_pairs,rank*num_pairs+num_pairs,len(p)))
    
    for cp in p_p:
        
        try:
            wf1,wf2,src,adjt = paths_input(cp,source_config,step,kernelrun)
            print(wf1,wf2,src)
        
            kernel,corr = paths_output(cp,source_config,step)
            
            
        except:
            print('Could not determine correlation for: ')
            print(cp)
            continue
            
        if os.path.exists(corr) and not kernelrun:
            continue

        if os.path.exists(kernel) and kernelrun:
            continue

        if not os.path.exists(adjt) and kernelrun:
            print('No adjoint source found for:')
            print(os.path.basename(adjt))
            continue

        
            #if int(step) == 0:
               # if source_config['ktype'] == 'td':

                    #print('Time domain preliminary kernel...')
        g1g2_corr(wf1,wf2,corr,kernel,adjt,src,source_config,kernelrun=kernelrun)

            #     elif source_config['ktype'] == 'fd':
            #         print('Frequency domain preliminary kernel...')
            #         g1g2_corr_fd(wf1,wf2,c,c_int,src,source_config)
            # #else:
            #     corr(wf1,wf2,c,c_int,src,source_config,kernelrun=kernelrun)
        #
        #corr(c,src,c_int)
                
            
            
        #except:
        #    print('Could not determine correlation for: ')
        ##    print(cp)
                
                #nsrc = os.path.join(source_conf['project_path'],
                #    source_conf['source_name'],
                #    'sourcemodel.h5')
                #
                
                
                
# def g1g2_corr_fd(wf1,wf2,corr_file,corr_int_file,src,source_conf):
    
#     #ToDo: Take care of saving metainformation
#     #ToDo: Think about how to manage different types of sources (numpy array vs. get from configuration -- i.e. read source from file as option)
#     #ToDo: check whether to include autocorrs from user (now hardcoded off)
#     #ToDo: Parallel loop(s)
#     #ToDo tests

#     ntimes, n, n_corr = get_ns(wf1,source_conf)
#     print('------')
#     print('ntimes: '+str(ntimes))
#     if n%2 == 0:
#         n_freq = n / 2 + 1
#     else:
#         n_freq = (n+1) / 2 
#     print('n_freq: '+str(n_freq))
#     print('n: '+str(n))
#     print('------')
#     taper = cosine_taper(ntimes,p=0.05)
    
#     print('got time lengths, taper')
    
    
#     # N needs to be saved for future zero padding.
#     save_n = open(os.path.join(source_conf['project_path'],source_conf['source_name'],
#     'n.txt'),'w')
#     print(str(n),file=save_n)
#     save_n.close()
#     # Zero padding should add an acausal part to the Green's function
#     mid_index = zero_buddy(ntimes,n,causal_function=True)
    
    
    
#     with WaveField(wf1) as wf1, WaveField(wf2) as wf2:
        
#         print('opened wave field files')
        
#         if wf1.stats['Fs'] != wf2.stats['Fs']:
#             msg = 'Sampling rates of synthetic green\'s functions must match.'
#             raise ValueError(msg)
        
#         # initialize new hdf5 files for correlation and green's function correlation
#         #with wf1.copy_setup(corr_file,nt=n_freq,complex=True) as correl,\
#         #NoiseSource(src) as nsrc:
#         correl = wf1.copy_setup(corr_file,nt=n_freq,complex=True)
#         nsrc = NoiseSource(src)
#         print('opened noise source and output file.')
        
#         correlation = np.zeros(n_corr)
#         print('initialized correlations.')
        
        
        
#         # Loop over source locations
#         with click.progressbar(range(wf1.stats['ntraces']),\
#         label='Correlating...' ) as ind:
#             for i in ind:
                
#                 s1 = np.ascontiguousarray(wf1.data[i,:]*taper)
#                 s2 = np.ascontiguousarray(wf2.data[i,:]*taper)
                
#                 #s1[mid_index:mid_index+ntimes] = wf1.data[i,:]*taper
#                 #s2[mid_index:mid_index+ntimes] = wf2.data[i,:]*taper
                
#                 spec1 = np.fft.rfft(s1,n)
#                 spec2 = np.fft.rfft(s2,n)
                
              
#                 g1g2_tr = np.multiply(spec1,np.conjugate(spec2))
                
#                 # extract Green's function correlation here
#                 # Save only as much as the adjoint source will be long.
#                 # This has to be done only once.
#                 correl.data_i[i,:] = np.imag(g1g2_tr).astype(np.float32)
#                 correl.data_r[i,:] = np.real(g1g2_tr).astype(np.float32)
#                 correl.file.flush()
                
#                 c = np.multiply(spec1,np.conjugate(spec2))
#                 c = np.multiply(c,nsrc.get_spect(i))                
#                 correlation += my_centered(np.fft.ifftshift(np.fft.irfft(c,n)),n_corr)
                
#         trace = Trace()
#         trace.stats.sampling_rate = wf1.stats['Fs']
#         trace.data = correlation
#         trace.write(filename=corr_int_file,format='SAC')    
#         correlation = correl.space_integral()
     
                

# def corr(wf1,wf2,corr_file,corr_int_file,src,source_conf):
    
#     #ToDo: Take care of saving metainformation
#     #ToDo: Think about how to manage different types of sources (numpy array vs. get from configuration -- i.e. read source from file as option)
#     #ToDo: check whether to include autocorrs from user (now hardcoded off)
#     #ToDo: Parallel loop(s)
#     #ToDo tests
    
#     nt, n, n_corr = get_ns(wf1,source_conf)
#     if source_config['ktype'] == 'fd':
#         n_corr = nt
#     taper = cosine_taper(nt,p=0.05)
    
#     with WaveField(wf1) as wf1, WaveField(wf2) as wf2:
        
        
#         if wf1.stats['Fs'] != wf2.stats['Fs']:
#             msg = 'Sampling rates of synthetic green\'s functions must match.'
#             raise ValueError(msg)
        
        
#         # initialize new hdf5 files for correlation and green's function correlation
#         with NoiseSource(src) as nsrc:
#             correlation = np.zeros(n_corr)
#             # Loop over source locations
            
#             for i in range(wf1.stats['ntraces']):
               
#                 s1 = np.ascontiguousarray(wf1.data[i,:]*taper)
#                 s2 = np.ascontiguousarray(wf2.data[i,:]*taper)
                
#                 spec1 = np.fft.rfft(s1,n)
#                 spec2 = np.fft.rfft(s2,n)
                
              
#                 g1g2_tr = np.multiply(spec1,np.conjugate(spec2))
#                 c = np.multiply(g1g2_tr,nsrc.get_spect(i))             
#                 correlation += my_centered(np.fft.ifftshift(np.fft.irfft(c,n)),n_corr)
                    
#             trace = Trace()
#             trace.stats.sampling_rate = wf1.stats['Fs']
#             trace.data = correlation
#             trace.write(filename=corr_int_file,format='SAC')    
#             correlation = correl.space_integral()
             