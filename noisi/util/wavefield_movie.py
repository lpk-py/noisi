# This example uses a MovieWriter directly to grab individual frames and
# write them to a file. This avoids any event loop integration, but has
# the advantage of working with even the Agg backend. This is not recommended
# for use in an interactive setting.
# -*- noplot -*-

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
from noisi import WaveField
import sys
from mpl_toolkits.basemap import Basemap
from matplotlib.mlab import griddata
import matplotlib.tri as tri     


#################################
v = 1.
stations = [(0.,0.)]
lonmin=-120.
lonmax=120.
latmin=-60.
latmax=60.
latc=0.0
lonc=0.0
resolution = 4
fps = 0.5

wf = WaveField(sys.argv[1])
t_min = float(sys.argv[2])
t_max = float(sys.argv[3])
t_step = float(sys.argv[4])
filename = sys.argv[5]
#################################

FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Wavefield', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=fps, metadata=metadata)

fig = plt.figure()
plt.subplot(111)

map_x = wf.sourcegrid[0]
map_x = map_x[0::resolution]
map_y = wf.sourcegrid[1]
map_y = map_y[0::resolution]
triangles = tri.Triangulation(map_x,map_y)
m = Basemap(rsphere=6378137,resolution='c',projection='cyl',lat_0=latc,           lon_0=lonc,llcrnrlat=latmin,urcrnrlat=latmax,
    llcrnrlon=lonmin,urcrnrlon=lonmax)
    
m.drawcoastlines(linewidth=0.5)
m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0]) # draw parallels
m.drawmeridians(np.arange(-180,210,60.),labels=[0,0,0,1]) # draw meridians
for sta in stations:
    m.plot(sta[0],sta[1],'rv',markersize=10,latlon=True)


with writer.saving(fig, filename, 100):
    for t in np.arange(t_min,t_max,t_step):
        print t
        map_z = wf.get_snapshot(t,resolution=resolution)
        #if globe:
        #    map_z = np.append(map_z,map_z[0])
        plt.tripcolor(triangles, map_z/np.max(np.abs(map_z)),                         shading='flat', vmin=-v,vmax=v, cmap=plt.cm.bwr)
        writer.grab_frame()
    

    
                             
    
        
