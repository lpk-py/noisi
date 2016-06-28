# plotting on the map
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import matplotlib.tri as tri    
import numpy as np
    
def plot_grid(map_x,map_y,map_z,stations=[],v=1.2,globe=False,outfile=None):
    
    m = Basemap(rsphere=6378137,resolution='c',projection='cyl',lat_0=0.,           lon_0=0.,llcrnrlat=np.min(map_y),urcrnrlat=np.max(map_y),
    llcrnrlon=np.min(map_x),urcrnrlon=np.max(map_x))
    if globe:
        map_x = np.append(map_x,map_x[0])
        map_y = np.append(map_y,map_y[0])
        map_z = np.append(map_z,map_z[0])
    triangles = tri.Triangulation(map_x,map_y)
    # tripcolor plot.
    plt.figure()
    plt.subplot(111)
    plt.gca().set_aspect('equal')
   
    plt.tripcolor(triangles, map_z/np.max(np.abs(map_z)),                         shading='flat', vmin=-v,vmax=v, cmap=plt.cm.bwr)
    m.colorbar(location='bottom')
    m.drawcoastlines(linewidth=0.5)
    m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
    m.drawmeridians(np.arange(-180,210,60.),labels=[0,0,0,1]) # draw meridians
    #draw station locations
    for sta in stations:
        m.plot(sta[0],sta[1],'rv',latlon=True)
    
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile,format='png')
    
def plot_sourcegrid(gridpoints,**kwargs):

    plt.figure()
    plt.subplot(111)
    m = Basemap(rsphere=6378137,**kwargs)
    m.drawcoastlines()
    
    m.plot(gridpoints[0],gridpoints[1],'+',markersize=10.,latlon=True)
    plt.show()
    

def plot_window(correlation, window, measurement):
    
    
    maxlag = correlation.stats.npts * correlation.stats.delta
    lag = np.linspace(-maxlag,maxlag,correlation.stats.npts)
    
    plt.plot(lag,correlation.data/np.max(np.abs(correlation.data)))
    plt.plot(lag,window/np.max(np.abs(window)),'--')
    plt.title(correlation.id)
    plt.text(0,-0.75,'Measurement value: %g' %measurement)
    plt.xlabel('Correlation Lag in seconds.')
    plt.ylabel('Normalized correlation and window.')
    
    plt.show()