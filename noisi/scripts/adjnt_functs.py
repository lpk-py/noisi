import numpy as np
from math import pi
#from noisi.scripts import measurements as rm
from noisi.util import windows as wn




def log_en_ratio_adj(corr_o,corr_s,g_speed,window_params):

    success = False

    window = wn.get_window(corr_o.stats,g_speed,window_params)
    win = window[0]
    #msr_o = rm.log_en_ratio(corr_o,g_speed,window_params)
    #msr_s = rm.log_en_ratio(corr_s,g_speed,window_params)
    data = wn.my_centered(corr_s.data,corr_o.stats.npts)


    if window[2] == True:
        sig_c = corr_s.data * win
        sig_a = corr_s.data * win[::-1]
        E_plus = np.trapz(np.power(sig_c,2))*corr_s.stats.delta
        E_minus = np.trapz(np.power(sig_a,2))*corr_s.stats.delta
        u_plus = np.multiply(sig_c,win)
        u_minus = np.multiply(sig_a,win[::-1])
        #adjt_src = 2./pi * (msr_s-msr_o) * (u_plus / E_plus - u_minus / E_minus)
        adjt_src = 2./pi * (u_plus / E_plus - u_minus / E_minus)
        success = True
    else:
        adjt_src = win-win+np.nan
    return adjt_src, success
    
def energy(corr_o,corr_s,g_speed,window_params):
    
    #msr_o = rm.energy(corr_o,g_speed,window_params)
    #msr_s = rm.energy(corr_s,g_speed,window_params)
    success = False
    
    window = wn.get_window(corr_o.stats,g_speed,window_params)
    
    if window_params['causal_side']:
        win = window[0]
    else:
        win = window[0][::-1]

    if window[2]:    
        u = np.multiply(np.power(win,2),corr_s.data)
        #adjt_src = 2./(msr_s-msr_o) * u
        adjt_src = u
        success = True
    else:
        adjt_src = win-win+np.nan
    
    return adjt_src, success
    
    
    

def get_adj_func(mtype):
    if mtype == 'ln_energy_ratio':
        func = log_en_ratio_adj
    
    elif mtype == 'energy_diff':
        func = energy
        
    else:
        msg = 'Measurement functional %s not currently implemented.' %mtype
        raise ValueError(msg)
    return func  