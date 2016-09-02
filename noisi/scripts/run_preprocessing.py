from mpi4py import MPI
from noisi import WaveField
import os
from glob import glob
import json

def run_preprocessing(source_config):

    configfile = os.path.join(source_config['project_path'],
    'config.json')
    config = json.load(open(configfile))
    
    files = glob(os.path.join(config['wavefield_path'],'*.h5'))
    processed_path = os.path.join(source_config['source_path'],
        'wavefield_processed')
    
    if not os.path.exists(processed_path):
        os.mkdir(processed_path)
        
    # very simple embarrassingly parallel loop
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm. Get_rank()

    files = files[rank::size]
    
    for file in files:

        newfile = os.path.join(processed_path, os.path.basename(file))

        if os.path.exists(newfile):
            print "File {} was already processed, skipping.".format(os.path.basename(file))
        print "Preprocessing {}".format(os.path.basename(file))
        

         
        if source_config['preprocess_truncate_sec'] is not None:
            
            # truncating
            
            with WaveField(file) as wf:
                wf.truncate(newfile,float(source_config['preprocess_truncate_sec']))


        if source_config['preprocess_decimation_factor'] is not None:

            # Already truncated file?
            if os.path.exists(newfile):
                newfile_temp = newfile + '.temp'
                with WaveField(newfile) as wf:
                    wf.decimate(decimation_factor=source_config['preprocess_decimation_factor'],
                                outfile=newfile_temp,
                                taper_width=0.005)
                os.system("mv {} {}".format(newfile_temp,newfile))
            else:
                with WaveField(file) as wf:
                    wf.decimate(decimation_factor=source_config['preprocess_decimation_factor'],
                                outfile=newfile,
                                taper_width=0.005)
                




            
        if source_config['preprocess_filter_kind'] == 'bandpass':

            # The file has been written previously by wavefield.truncate
            if os.path.exists(newfile):
                with WaveField(newfile) as wf:
                    wf.filter_all(
                        source_config['preprocess_filter_kind'],
                        overwrite=True,
                        freqmin=source_config['preprocess_filter_params'][0],
                        freqmax=source_config['preprocess_filter_params'][1],
                        corners=source_config['preprocess_filter_params'][2],
                        zerophase=source_config['preprocess_filter_params'][3])

            else:
                # The file still has to be written
                with WaveField(file) as wf:
                    wf.filter_all(
                        source_config['preprocess_filter_kind'],
                        overwrite=False,
                        freqmin=source_config['preprocess_filter_params'][0],
                        freqmax=source_config['preprocess_filter_params'][1],
                        corners=source_config['preprocess_filter_params'][2],
                        zerophase=source_config['preprocess_filter_params'][3],
                        outfile=newfile)



            # filtering type,overwrite=False,zerophase=True,**kwargs
            #with WaveField(newfile) as wf:
                
                
                
            #    wf.filter_all()
            
        
        
                
    
    