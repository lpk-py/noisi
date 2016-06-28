import numpy as np

# ToDo: Might not be necessary to have a separate class for this, but could be incorporated into the NoiseSource class.
# ToDo: Docs
class MyBasis(object):
    
    def __init__(self,btype,params):
        """
        possible types and parameters:
        'boxcar' (like, linear frequency or period or time bins), params: nr. of bins
        'gaussian' (Gaussians; params: mean and standard deviation as tuple)
        #ToDo: Maybe: Include default parameters
        """


        
        # Is the type needed?
        self.btype = btype
        self.params = params
        

        #if isinstace(params,np.ndarray):
        #    self.params = params
        #elif isinstance(params,str):
        #    try:
        #        self.params = np.load(params)
        #    except:
        #        raise IOError('Basis function file not found.')
        #else:
        #    raise TypeError
        #

    def to_grid(self,grid):
        """
        For basis functions defined in a simple manner -- e.g. 'Five frequency
        bins' return the basis functions on the specified grid. This depends on the type of basis function.
        :type grid: np.ndarray
        :param grid:  The grid on which our simulation lives. This can be a time- or frequency axis, or a coordinate grid. ToDo: Think about how to handle meshgrids vs. coordinate lists (this is only a rearrangement).
        """

        basis_grid = None

        if self.btype == 'boxcar': #This is mostly for understanding the principle.
            # Get the actual value of self.params! self.params is mutable. No modification wanted...
            basis_shape = [self.params[0]] #if self.params[0] > 1. else []
            basis_shape.extend(i for i in np.shape(grid))
            
            basis_grid = np.zeros(basis_shape)
            print np.shape(basis_grid)
            step = len(grid) // basis_shape[0]
            for i in np.arange(basis_shape[0]):
                basis_grid[i,i*step:(i+1)*step] = 1.

        

        #if self.btype == 'gaussian':
        #    basis_shape = [len(self.params)] #if self.params not None else [1.]
        #    basis_shape.extend(i for i in np.shape(grid))
        #    basis_grid = np.zeros(basis_shape)
        #    for i in np.arange(len(self.params)):
        #        mu = self.params[i][0]
        #        sig2 = self.params[i][1]**2
        #        # ToDo: Use normalized Gaussians yes or no?
        #        basis_grid[i,:] += np.exp(-1*(grid-mu)**2/(2*sig2))
        
        if np.shape(basis_grid)[0] > 1.:
            return basis_grid
        else:
            return basis_grid[0]