import numpy as np
import json
import os
import shutil
import io

def grid_to_specfem_stations(grid,spec_dir):
    """
    Write noisesource grid to disk as specfem compatible station list.
    """
    
    
    fid = open('temp.txt','w')
    for i in range(len(grid[0,:])):
        fid.write('%08g SRC %10.8f  %10.8f 0.0 0.0\n'\
            %(i,grid[1,i],grid[0,i]))
    
    fid.close()
    
    for dir in os.listdir(spec_dir):
        dst = os.path.join(spec_dir,dir,'DATA','STATIONS')
        print dst
        shutil.copy('temp.txt',dst)
    os.remove('temp.txt')
    
def stations_to_cmtsolutions(stationlist,hdur,outdir):
    

    fid = open(stationlist,'r')
    stationlist = fid.read().split('\n')
    
    for i in range(len(stationlist)):
        
        station = stationlist[i]
        if station =='': continue
        
        print station
        info = station.split()
        
        id = info[0].strip() + '.' + info[1].strip()
        os.mkdir(os.path.join(outdir,id))
        os.mkdir(os.path.join(outdir,id,'DATA'))
        os.mkdir(os.path.join(outdir,id,'OUTPUT_FILES'))
        os.mkdir(os.path.join(outdir,id,'DATABASES_MPI'))
        
        eventfid = open(os.path.join(outdir,id,'DATA','CMTSOLUTION'),'w')
        
        eventfid.write('*** 2000  1  1  1 01 01.00 '+info[3].strip()+\
        ' '+info[4].strip()+' '+info[5].strip()+' 0.0 0.0 ***\n')
        eventfid.write('event name:    %s   \n' %id)
        eventfid.write('time shift:    0.0000   \n')
        eventfid.write('half duration:    %s   \n' %str(hdur))
        eventfid.write('latitude:    %s   \n' %str(info[3].strip()))
        eventfid.write('longitude:    %s   \n' %str(info[4].strip()))
        eventfid.write('depth:    %s   \n' %str(info[5].strip()))
        eventfid.write('Mrr:    0.0000000   \n')
        eventfid.write('Mtt:    0.0000000   \n')
        eventfid.write('Mpp:    0.0000000   \n')
        eventfid.write('Mrt:    0.0000000   \n')
        eventfid.write('Mrp:    0.0000000   \n')
        eventfid.write('Mtp:    0.0000000   ') 



def prepare_specfem_input(configfile):
    
    with io.open(configfile,'r') as fh:
        config = json.load(fh)
    
    spec_dir = os.path.join(config['project_path'],'specfem_input') 
    os.mkdir(spec_dir)
    grid = np.load(os.path.join(config['project_path'],'sourcegrid.npy'))
    stations = os.path.join(config['project_path'],'stations.txt')
    
    
    
    stations_to_cmtsolutions(stations,config["hdur_pointsource"],spec_dir)
    grid_to_specfem_stations(grid,spec_dir)
    
