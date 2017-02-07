
# coding: utf-8


import numpy as np
from obspy.geodetics import gps2dist_azimuth
from obspy.signal.invsim import cosine_taper
import matplotlib.pyplot as plt
import h5py
from noisi import WaveField
import json
from glob import glob
import os
from scipy.fftpack import next_fast_len
from scipy.signal import hann


##################################################################
# USER INPUT
##################################################################
# path to project
projectpath = '../'
sourcepath = '.'

df = 0.1
nt = 2883

# geography - a sequence of distributions 'homogeneous', 'ocean',
# 'gaussian_blob' in any order. The order has to match with the 
# order od the list of spectra in params_spectra, i.e. the first 
# distribution will be assigned the first spectrum, the second 
# distribution the second spectrum, etc. 
# Similarly, the first 'gaussian_blob' will be assigned the first
# set of parameters in params_gaussian_blobs, and so on.
distributions = ['homogeneous','homogeneous']
coastres = 'c'
#only_ocean = False
#gaussian_blobs = False

# Geographic gaussian blobs. Will be ignored if no 'gaussian_blob'
# is in the list of distributions.
params_gaussian_blobs = [{'center':(-10.,0.),'sigma_radius_m':2000000.,
'rel_weight':2.,'only_ocean':True}]

# spectra: Length N must match the length of distributions. The order must match as well.
params_spectra = [{'central_freq':0.005,'sigma_freq':0.001,'weight':1.},
                           {'central_freq':0.015,'sigma_freq':0.0025,'weight':1.} ]


                          #{'central_freq':0.075,'sigma_freq':0.1,'weight':10.}]

###############################################################################

grd  = np.load(os.path.join(projectpath,'sourcegrid.npy'))
ntraces = np.shape(grd)[-1]
print 'Loaded source grid'

# config = json.load(open(os.path.join(projectpath,'config.json')))
# source_config = json.load(open(os.path.join(sourcepath,'source_config.json')))
# print 'Loaded config files.'

# if source_config['preprocess_do']:
#     ext = '*.h5'
#     wavefield_path = os.path.join(sourcepath,'wavefield_processed')
# else:
#     ext = '*.h5'
#     wavefield_path = config['wavefield_path']


# wfs = glob(os.path.join(wavefield_path,ext))
# if wfs != []:
# 	print 'Found wavefield.'


# with WaveField(wfs[0]) as wf:
#     df = wf.stats['Fs']
    # nt = wf.stats['nt']
    # # The number of points for the fft is larger due to zeropadding --> apparent higher frequency sampling\n",
    # n = next_fast_len(2*nt-1)
n = next_fast_len(2*nt-1)    
freq = np.fft.rfftfreq(n,d=1./df)
taper = cosine_taper(len(freq),0.01)
print 'Determined frequency axis.'

def get_distance(grid,location):
    def f(lat,lon,location):
        return abs(gps2dist_azimuth(lat,lon,location[0],location[1])[0])
    dist = np.array([f(lat,lon,location) for lat,lon in zip(grid[1],grid[0])])
    return dist

# Use Basemap to figure out where ocean is
def get_ocean_mask():
    print 'Getting ocean mask...'
    from mpl_toolkits.basemap import Basemap
    m = Basemap(rsphere=6378137,resolution=coastres,projection='cea',lat_0=0.,
                lon_0=0.,llcrnrlat=-90.,urcrnrlat=90.,llcrnrlon=-180.,urcrnrlon=180.)
    (x,y) = m(grd[0],grd[1])
    ocean_mask = map(lambda (x,y): not m.is_land(x,y),zip(x,y))
    return ocean_mask


   

def get_geodist(disttype,gaussian_params=None):

    if disttype == 'gaussian':
        dist = get_distance(grd,gaussian_params['center'])
        gdist = np.exp(-(dist)**2/(2*gaussian_params['sigma_radius_m']**2))

        if gaussian_params['only_ocean']:
            if not 'ocean_mask' in locals():
                ocean_mask = get_ocean_mask()
            gdist *= ocean_mask

        return gdist

    elif disttype == 'homogeneous':
        return np.ones(ntraces)

    elif disttype == 'ocean':
        if not 'ocean_mask' in locals():
            ocean_mask = get_ocean_mask()
        return ocean_mask


def get_spectrum(sparams):
    spec = taper*np.exp(-(freq-sparams['central_freq'])**2/(2*sparams['sigma_freq']**2))
    return spec / np.max(np.abs(spec))




#########################
# Create the source distr
#########################

# geography
num_bases = len(distributions)
gauss_cnt = 0
#if gaussian_blobs:
#    num_bases += len(params_gaussian_blobs)

basis_geo = np.zeros((num_bases,ntraces))
#weights_geo = np.zeros(np.shape(basis_geo)[0])
basis_spec = np.zeros((num_bases,len(freq)))
#weights_spec = np.ones((np.shape(grd)[-1],np.shape(basis2)[0]))
print 'Filling distribution...'

for i in range(num_bases):

    if distributions[i] =='gaussian':

        gaussparams = params_gaussian_blobs[gauss_cnt]
        gauss_cnt += 1
        basis_geo[i,:] = get_geodist('gaussian',gaussparams)
        basis_spec[i,:] = get_spectrum(params_spectra[i])

    elif distributions[i] in ['ocean','homogeneous']:

        basis_geo[i,:] = get_geodist(distributions[i])
        basis_spec[i,:] = get_spectrum(params_spectra[i])
    else:
        raise NotImplementedError('Unknown geographical distributions. \
            Must be \'gaussian\', \'homogeneous\' or \'ocean\'.')

# # homogeneous layer
# basis1[0,:] = np.ones(ntraces) 

# if only_ocean:
#     ocean_mask = np.array(get_ocean_mask()).astype(int)
#     basis1[0,:] *= ocean_mask

#     # superimposed Gaussian blob(s)
# if gaussian_blobs:
#     i = 1
#     for blob in params_gaussian_blobs:
#         dist = get_distance(grd,blob['center'])
#         basis1[i,:] = np.exp(-(dist)**2/(2*blob['sigma_radius_m']**2))
        
#         if only_ocean:
#             basis1[i,:] *= ocean_mask
#         i+=1


#print 'Filling spectra...'  
# spectra
#basis2 = np.zeros((len(params_gaussian_spectra),len(freq)))
# 'sort of hum gaussian'
#i = 0
#for spec in params_gaussian_spectra:
#    basis2[i,:] = taper*np.exp(-(freq-spec['central_freq'])**2/(2*spec['sigma_freq']**2))
# This normalization means different integrals...
#    basis2[i,:] /= np.max(np.abs(basis2[0,:]))
#    i+=1


######################
# set the weights
#####################
## geography
#weights1 = np.ones(np.shape(basis1)[0])
#
#if gaussian_blobs:
#    i = 1
#    for blob in params_gaussian_blobs:
#        weights1[i] = blob['rel_weight']
#        i+=1
#    if no_background:
#        weights1[0] = 0.
##print weights1
## spectra --- much harder to assign manually, since we need weights for every location. just assigning ones.\n",
#weights2 = np.ones((np.shape(grd)[-1],np.shape(basis2)[0]))
#i=0
#for spec in params_gaussian_spectra:
#    weights2[:,i] *= spec['weight']
#
print 'Plotting...'
from noisi.util import plot
#
#
#
# geograhic_distr = np.sum(basis_geo,axis=0)
# plot.plot_grid(grd[0],grd[1],geograph_distr,
# outfile = os.path.join(sourcepath,'geog_distr_startingmodel.png'))

for i in range(num_bases):
    plot.plot_grid(grd[0],grd[1],basis_geo[i,:],
    outfile = os.path.join(sourcepath,'geog_distr_basis{}.png'.format(i)))

#
#
plt.figure()
#plt.semilogx(freq,np.sum(basis_spec,axis=0),linewidth=2)
for i in range(num_bases):
    plt.semilogx(freq,basis_spec[i,:],'--')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Source power (scaled)')
plt.savefig(os.path.join(sourcepath,'freq_distr_startingmodel.png'))
#
#
# Save to an hdf5 file

with h5py.File(os.path.join(sourcepath,'step_0','starting_model_new.h5'),'w') as fh:
    fh.create_dataset('coordinates',data=grd.astype(np.float32))
    fh.create_dataset('frequencies',data=freq.astype(np.float32))
    fh.create_dataset('geogr_weights',data=basis_geo.astype(np.float32))
#    fh.create_dataset('distr_weights',data=weights1.astype(np.float32))
    fh.create_dataset('spect_basis',data=basis_spec.astype(np.float32))
#    fh.create_dataset('spect_weights',data=weights2.astype(np.float32))
print 'WARNING: All arrays are stored in single precision.'
print 'Done.'

