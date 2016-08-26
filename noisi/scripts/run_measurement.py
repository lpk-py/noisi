import os
import numpy as np
import pandas as pd
from math import log, pi
import click
import json
from scipy.signal import hilbert
from glob import glob
from obspy import read, Trace
from obspy.geodetics import gps2dist_azimuth
import matplotlib.pyplot as plt
#ToDo plot if requested.
from noisi.scripts import measurements as rm
#from noisi.scripts import adjnt_functs as af
from noisi.util.windows import get_window, my_centered, snratio
# Get and return measurement as a table or something.


def get_station_info(stats):

    sta1 = '{}.{}.{}.{}'.format(stats.network,stats.station,stats.location,
    stats.channel)
    sta2 = '{}.{}.{}.{}'.format(stats.sac.kuser0.strip(),stats.sac.kevnm.strip(),
    stats.sac.kuser1.strip(),stats.sac.kuser2.strip())
    lat1 = stats.sac.stla
    lon1 = stats.sac.stlo
    lat2 = stats.sac.evla
    lon2 = stats.sac.evlo
    dist = stats.sac.dist
    az = gps2dist_azimuth(lat1,lon1,lat2,lon2)[2]
    
    
    return([sta1,sta2,lat1,lon1,lat2,lon2,dist,az])

def get_synthetics_filename(obs_filename,synth_location='',fileformat='sac',synth_channel_basename='MX'):

    inf = obs_filename.split('--')

    if len(inf) == 1:
        # old station name format
        inf = obs_filename.split('.')
        net1 = inf[0]
        sta1 = inf[1]
        cha1 = inf[3]
        net2 = inf[4]
        sta2 = inf[5]
        cha2 = inf[7]
    elif len(inf) == 2:
        # new station name format
        inf1 = inf[0].split('.')
        inf2 = inf[1].split('.')
        net1 = inf1[0]
        sta1 = inf1[1]
        net2 = inf2[0]
        sta2 = inf2[1]
        cha1 = inf1[3]
        cha2 = inf2[3]


    cha1 = synth_channel_basename + cha1[-1]
    cha2 = synth_channel_basename + cha2[-1]

    synth_filename = '{}.{}.{}.{}--{}.{}.{}.{}.{}'.format(net1,sta1,synth_location,
        cha1,net2,sta2,synth_location,cha2,fileformat)
    print(synth_filename)
    return(synth_filename)





def measurement(source_config,mtype,step,**options):
    
    """
    Get measurements on noise correlation data and synthetics. 
    options: g_speed,window_params (only needed if mtype is ln_energy_ratio or enery_diff)
    """
    step_n = 'step_{}'.format(int(step))
    
    
    step_dir = os.path.join(source_config['source_path'],
    step_n)
    
    files = [f for f in os.listdir(os.path.join(source_config['source_path'],
    'observed_correlations')) ]
    files = [os.path.join(source_config['source_path'],
    'observed_correlations',f) for f in files]
    
    synth_dir = os.path.join(step_dir,'corr')
    
    
    columns = ['sta1','sta2','lat1','lon1','lat2','lon2','dist','az',
    'obs','l2_norm','snr']
    measurements = pd.DataFrame(columns=columns)
    
    
    
    if files == []:
        msg = 'No input found!'
        raise ValueError(msg)
    
    i = 0
    with click.progressbar(files,label='Taking measurements...') as bar:
        
        for f in bar:
            
            try: 
                tr_o = read(f)[0]
            except:
                print('\nCould not read data: '+os.path.basename(f))
                i+=1
                continue
            try:
                synth_filename = get_synthetics_filename(os.path.basename(f))
                tr_s = read(os.path.join(synth_dir,synth_filename))[0]
            except:
                print('\nCould not read synthetics: ' + synth_filename)
                i+=1
                continue
            tr_s.stats.sac = tr_o.stats.sac.copy() #ToDo: Give the stats to this thing before!
            tr_s.data = my_centered(tr_s.data,tr_o.stats.npts)    
            # Get all the necessary information
            info = get_station_info(tr_o.stats)
           
            # Take the measurement
           
            func = rm.get_measure_func(mtype)
            msr_o = func(tr_o,**options)
            msr_s = func(tr_s,**options)
            
            # timeseries-like measurements:
            if mtype in ['envelope','windowed_envelope','waveform',\
            'windowed_waveform']:
                l2_so = np.trapz(0.5*(msr_s-msr_o)**2) * tr_o.stats.delta
                msr = np.nan
                snr = np.nan
            # single value measurements:
            else:
                l2_so = 0.5*(msr_s-msr_o)**2
                msr = msr_o
                snr = snratio(tr_o,**options)
                 
            info.extend([msr,l2_so,snr])
            measurements.loc[i] = info

            # step index
            i+=1
    
    filename = '{}.measurement.csv'.format(mtype)
    measurements.to_csv(os.path.join(step_dir,filename),index=None)

def run_measurement(source_configfile,measr_configfile,step):


    
    source_config=json.load(open(source_configfile))
    measr_config=json.load(open(measr_configfile))
    
    mtype = measr_config['mtype']
    
    # TODo all available misfits --  what parameters do they need (if any.)
    if measr_config['mtype'] in ['ln_energy_ratio','energy_diff']:
        

        g_speed                         =    measr_config['g_speed']
        window_params                   =    {}
        window_params['hw']             =    measr_config['window_params_hw']
        window_params['sep_noise']      =    measr_config['window_params_sep_noise']
        window_params['win_overlap']    =    measr_config['window_params_win_overlap']
        window_params['wtype']          =    measr_config['window_params_wtype']
        window_params['causal_side']    =    measr_config['window_params_causal']
        window_params['plot']           =    measr_config['window_plot_measurements']
    
    measurement(source_config,mtype,step,g_speed=g_speed,window_params=window_params)
    