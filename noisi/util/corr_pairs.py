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
            if '' not in [sta_0,sta]:
                corr_pairs.append([sta_0,sta])
    return corr_pairs