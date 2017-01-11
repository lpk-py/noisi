# plotting on the map
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import matplotlib.tri as tri    
import numpy as np
    
def plot_grid(map_x,map_y,map_z,stations=[],v=None,globe=False,
    outfile=None,title=None,shade='flat',cmap=None,
    sequential=False,v_min=None,normalize=True,coastres='c',proj='cyl',
    lonmin=None,lonmax=None,latmin=None,latmax=None,mode='interp'):
    
    lat_0  = 0.5*(map_y.max()-map_y.min())

    if lonmin == None:
        lonmin = np.min(map_x)
    if lonmax == None:
        lonmax = np.max(map_x)
    if latmax == None:
        latmax = np.max(map_y)
    if latmin == None:
        latmin = np.min(map_y)



    m = Basemap(rsphere=6378137,resolution=coastres,
    projection=proj,lat_0=0.,
    lon_0=0.,llcrnrlat=latmin,urcrnrlat=latmax,
    llcrnrlon=lonmin,urcrnrlon=lonmax)


    plt.figure()
    plt.subplot(111)
    plt.gca().set_aspect('equal')
    if title is not None:
        plt.title(title)

    

    if normalize:
        map_z /= np.max(np.abs(map_z))
    
    if v is None:
        v = np.max(map_z)

    if sequential:
        cm = plt.cm.magma
        if v_min == None:
            v_min = 0.
    else:
        cm = plt.cm.bwr
        v_min =-v

    if cmap is not None:
        cm = cmap
    
    print('max. value on map: %g' %map_z.max())
    if mode == 'interp':
        triangles = tri.Triangulation(map_x,map_y)
        # tripcolor plot.
        plt.tripcolor(triangles,map_z,shading=shade, vmin=v_min,vmax=v,cmap=cm)
        m.colorbar(location='bottom',pad=0.4)
    elif mode == 'srclocs':
        plt.scatter(map_x,map_y,marker='o',c='white')
    elif mode == 'srcdots':
        
        
        colors = cm(map_z)

        plt.scatter(map_x,map_y,marker='o',c=colors)
        m.colorbar(location='bottom',pad=0.4)
    
    if globe:
        m.drawcoastlines(linewidth=0.5)
    else:
        m.drawcoastlines(linewidth=2.5)
    if globe:
        m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
        m.drawmeridians(np.arange(-180,210,60.),labels=[0,0,0,1]) # draw meridians
    else:
        d_lon = round(abs(lonmax-lonmin) / 5.)
        d_lat = round(abs(latmax-latmin) / 5.)
        parallels = np.arange(latmin,latmax,d_lat).astype(int)
        meridians = np.arange(lonmin,lonmax,d_lon).astype(int)
        m.drawparallels(parallels,labels=[1,0,0,0]) # draw parallels
        m.drawmeridians(meridians,labels=[0,0,0,1])

    #draw station locations
    for sta in stations:
        m.plot(sta[0],sta[1],'kv',markersize=10,latlon=True)
    plt.show()
    if outfile is None:
        plt.show()
    else:
        plt.savefig(outfile,format='png')
    
def plot_sourcegrid(gridpoints,**kwargs):

    plt.figure()
    plt.subplot(111)
    m = Basemap(rsphere=6378137,**kwargs)
    m.drawcoastlines()
    
    m.plot(gridpoints[0],gridpoints[1],marker='+',markersize=10.,latlon=True)
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