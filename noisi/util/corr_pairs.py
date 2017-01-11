from glob import glob
import os

def define_correlationpairs(proj_dir,auto_corr=False,only_observed=True,channel='*'):
    """
    Match correlation pairs.
    :param proj_dir: Path to project directory, where a stations.txt file has to be located
    
    :param auto_corr: Include or exclude autocorrelations
    :return corr_pairs: A list of correlation pairs in the format
    [['net1 sta1 lat1 lon1','net2 sta2 lat2 lon2'],[...]]
    """
    
    
    stations = open(os.path.join(proj_dir,'stations.txt'))
    
    stations = stations.read().split('\n')
    i = 0
    corr_pairs = []
    
    while i < len(stations):
        sta_0 = stations[i].strip()
        if auto_corr:
            stas_other = stations[i:]
        else:
            stas_other = stations[i+1:]
        i += 1

        for sta in stas_other:

            if '' in [sta_0,sta]:
                continue
            corr_pairs.append([sta_0,sta])
    return corr_pairs.sort()


def rem_no_obs(stapairs,source_conf,directory,ignore_network=True):


    channel = source_conf['channel']
    channel = '??' + channel[-1]


    stapairs_new = []
    for i in range(len(stapairs)):
        # Check if an observation is actually available
        if stapairs[i] == '':
            break
       
        sta1 = '{}.{}.*.{}'.format(*(stapairs[i][0].split()[0:2]+[channel]))
        sta2 = '{}.{}.*.{}'.format(*(stapairs[i][1].split()[0:2]+[channel]))
        p_new = glob_obs_corr(sta1,sta2,directory,ignore_network)
        
        if p_new ==[]:
            continue
        stapairs_new.append(stapairs[i])
    return stapairs_new


def rem_fin_prs(stapairs,source_conf,step,kernelrun):

    """
    Remove those station pairs from the list for which correlation / kernel has already 
    been calculated.
    :param sta_pairs: List of all station pairs
    :param source_conf: source config dictionary
    :param step: step nr
    :param kernelrun: run for calculating kernels or correlations 
    """

    
    channel = source_conf['channel']

    if kernelrun:
        mod_dir = os.path.join(source_conf['source_path'],'step_{}'.format(step),'kern')
    else:
        mod_dir = os.path.join(source_conf['source_path'],'step_{}'.format(step),'corr')


    stapairs_new = []

    for sp in stapairs:
        id1 = sp[0].split()[0]+sp[0].split()[1]
        id2 = sp[1].split()[0]+sp[1].split()[1]

        if id1 < id2 :
            inf1 = sp[0].split()
            inf2 = sp[1].split()
        else:
            inf2 = sp[0].split()
            inf1 = sp[1].split()

        sta1 = "{}.{}..{}".format(*(inf1[0:2]+[channel]))
        sta2 = "{}.{}..{}".format(*(inf2[0:2]+[channel]))
        
        if kernelrun:
            kern_name = "{}--{}.npy".format(sta1,sta2)
            kern_name = os.path.join(mod_dir,kern_name)
            if not os.path.exists(kern_name):
                stapairs_new.append(sp)
        else:
            corr_name = "{}--{}.sac".format(sta1,sta2)    
            corr_name = os.path.join(mod_dir,corr_name)
            if not os.path.exists(corr_name):
                
                stapairs_new.append(sp)

    return stapairs_new

# Find the filename of the synthetic correlation from the one of the observed correlation
def get_synthetics_filename(obs_filename,dir,synth_location='',
    fileformat='sac',synth_channel_basename='??',ignore_network=True):

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


    sfilename = None

    if ignore_network:
        synth_filename1 = '*.{}.{}.{}--*.{}.{}.{}.{}'.format(sta1,synth_location,
        cha1,sta2,synth_location,cha2,fileformat)
        synth_filename2 = '*.{}.{}.{}--*.{}.{}.{}.{}'.format(sta2,synth_location,
        cha2,sta1,synth_location,cha1,fileformat)

        try: 
            sfilename = glob(os.path.join(dir,synth_filename1))[0]      
        except IndexError:
            try:
                sfilename = glob(os.path.join(dir,synth_filename2))[0]
            except IndexError:
                print('No synthetic file found for data:')
                print(obs_filename)

    else:
        synth_filename1 = '{}.{}.{}.{}--{}.{}.{}.{}.{}'.format(net1,sta1,synth_location,
        cha1,net2,sta2,synth_location,cha2,fileformat)

        try:
            sfilename = glob(os.path.join(dir,synth_filename1))[0]  
        
        except IndexError:
            print('No synthetic file found for data:')
            print(obs_filename)
        

    return sfilename


def glob_obs_corr(sta1,sta2,directory,ignore_network):


    inf1 = sta1.split('.')
    inf2 = sta2.split('.')

    sta1 = inf1[1]
    sta2 = inf2[1]

    cha1 = '??' + inf1[3][-1]
    cha2 = '??' + inf2[3][-1]

    if ignore_network:
        net1 = '*'
        net2 = '*'
    else:
        net1 = inf1[0]
        net2 = inf2[0]


    obs_filename1 = os.path.join(directory,'{}.{}.*.{}*{}.{}.*.{}.*'.format(net1,sta1,cha1,net2,sta2,cha2))
    obs_filename2 = os.path.join(directory,'{}.{}.*.{}*{}.{}.*.{}.*'.format(net2,sta2,cha2,net1,sta1,cha1))
    
    if ignore_network:
        obs_files = glob(obs_filename1)
        obs_files.extend(glob(obs_filename2))
    else:
        obs_files = glob(obs_filename1)


    return obs_files
    



# and now...find the observed filename from the synthetic one
# def get_observed_filenames(sta1,sta2,directory,stack=False,
#     t1=None,t2=None,fileformat=None,old_format=False,
#     channel_basename='LH',tag=None,ignore_network=True):

    
#     inf1 = sta1.split('.')
#     inf2 = sta2.split('.')
    
#     sta1 = inf1[1]
#     sta2 = inf2[1]

#     cha1 = channel_basename + inf1[3][-1]
#     cha2 = channel_basename + inf2[3][-1]

#     if not ignore_network:
#         net1 = inf1[0]
#         net2 = inf2[0]
#     else:
#         net1 = '*'
#         net2 = '*' 



#     if not stack:
#         stack = ''
#     else:
#         stack = '.stack'

#     if t1 is None:
#         t1=''
#     else:
#         t1 = '.'+t1
#     if t2 is None:
#         t2 = ''
#     else:
#         t2 = '.'+t2

#     if fileformat is None:
#         fileformat = '*'

#     if tag is None:
#         tag = '*'
    
#     if not old_format:
#         obs_filename = ('{}.{}.*.{}--{}.{}.*.{}{}{}{}.{}'.format(
#             net1,sta1,cha1,net2,sta2,cha2,stack,t1,t2,fileformat
#             ))
#     else:
#         obs_filename = ('{}.{}.*.{}.{}.{}.*.{}.ccc.{}.{}'.format(
#             net1,sta1,cha1,net2,sta2,cha2,tag,fileformat))

#     obs_filename = os.path.join(directory,obs_filename)
#     obs_filename = glob(obs_filename)


#     if ignore_network:
#         obs_filename2 = ('{}.{}.*.{}--{}.{}.*.{}{}{}{}.{}'.format(
#             net2,sta2,cha2,net1,sta1,cha1,stack,t1,t2,fileformat
#             ))
#         obs_filename2 = os.path.join(directory,obs_filename2)
#         obs_filename.extend(glob(obs_filename2))

#     return obs_filename

