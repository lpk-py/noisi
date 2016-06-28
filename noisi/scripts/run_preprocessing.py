from noisi import WaveField
import os
from glob import glob
import json

def run_preprocessing(source_config):

    configfile = os.path.join(source_config['project_path'],
    'config.json')
    config = json.load(open(configfile))
    
    files = glob(os.path.join(config['wavefield_path'],'*.h5'))
    
    
    # determine filter parameters
    #if source_config['preprocess_filter_kind'] == 'bandpass':
        
    
    
    for file in files:
        
        
        print "Preprocessing {}".format(os.path.basename(file))
        
         
        if source_config['preprocess_truncate_sec'] is not None:
            
            # truncating
            newfile = file+'_proc'
            with WaveField(file) as wf:
                wf.truncate(newfile,float(source_config['preprocess_truncate_sec']))
            
            # filtering type,overwrite=False,zerophase=True,**kwargs
            #with WaveField(newfile) as wf:
                
                
                
            #    wf.filter_all()
            
        
        
                
    
    