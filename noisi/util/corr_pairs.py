import os

def define_correlationpairs(proj_dir,auto_corr=False):
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
    return corr_pairs



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
                print(corr_name)
                stapairs_new.append(sp)

    return stapairs_new

