import numpy as np
import h5py

from scipy.stats import linregress
import os
from noisi.util.plot import plot_grid


class NoiseSource(object):
    """
   'model' of the noise source that comes in terms of a couple of basis 
    functions and associated weights. The NoiseSource object should contain a 
    function to determine weights from a (kernel? source model?), and to expand from weights and basis 
    functions.
    
    """
    
    
    def __init__(self,model):
            
        # Model is an hdf5 file which contains the basis and weights of the source model!
        
       
        try:
            self.model = h5py.File(model,'r+')
            self.src_loc = self.model['coordinates']
            self.freq = self.model['frequencies']
             
            # Presumably, these arrays are small and will be used very often --> good to have in memory.
            self.distr_basis = self.model['distr_basis'][:]
            self.spect_basis = self.model['spect_basis'][:]
            self.distr_weights = self.model['distr_weights'][:]
            self.spect_weights = self.model['spect_weights'][:]
            
            # The 'max amplitude' is probably needed often, and the distribution should not be too heavy to hold in memory, it has dimension 1 x number of sources
            self.source_distribution = self.expand_distr()
            
        except IOError:
            msg = 'Unable to open model file '+model
            raise IOError(msg)


        
    def __enter__(self):
        return self
    
    def __exit__(self,type,value,traceback):
        
        if self.model is not None:
            self.model.close()
            #ToDo: Check what parameters/data should be written before file closed

    def project_gridded(self):
        pass

    def expand_distr(self):
        return np.dot(self.distr_weights,self.distr_basis)
        

    def get_spect(self,iloc):
        # return one spectrum in location with index iloc
        # The reason this function is for one spectrum only is that the entire gridded matrix of spectra by location is most probably pretty big.
        spect = np.dot(self.spect_weights[iloc],self.spect_basis) 
        return self.source_distribution[iloc] * spect
    
    
    def plot(self,**options):
        
        # plot the distribution
        m = self.expand_distr()
        plot_grid(self.src_loc[0],self.src_loc[1],m,**options)


    
    # Note: Inefficient way of doing things! Whichever script needs the noise source field should rather look up things directly in the hdf5 file.
    # But: ToDo: This could be used internally to write to a file, rather than reading from.
    # Although: A problem to think about: noise source should behave the same, whether it is defined by model or by file. So maybe, since model will be the default option anyway, work with this!
    #def get_spectrum(self,iloc):
    #    # Return the source spectrum at location nr. iloc
    #    
    #    #if self.file is not None:
    #    #    return self.sourcedistr[iloc] * self.spectra[iloc,:]
    #    if self.model is not None:
    #        return self.spectr.
    #        # (Expand from basis fct. in model)
    #def get_sourcedistr(self,i):
    #    # Return the spatial distribution of max. PSD
    #    if self.file is not None:
    #        return np.multiply(self.sourcedistr[:],self.spectra[:,i])
    #    else:
    #        raise NotImplementedError
   
