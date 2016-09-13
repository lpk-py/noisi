import numpy as np
import pandas as pd
import os
import json
from math import isnan
from noisi.util.plot import plot_grid


def assemble_descent_dir(source_model,step,snr_min):



# where is the measurement database located?
	source_config=json.load(open(source_model))
	datadir = os.path.join(source_config['source_path'],'step_' + str(step))
	msrfile = os.path.join(datadir,"{}.measurement.csv".format(source_config['mtype']))

# Read in the csv files of measurement.

	data = pd.read_csv(msrfile)
	print(data)

# allocate the kernel array
	grd = np.load(os.path.join(source_config['project_path'],'sourcegrid.npy'))

	gradient = np.zeros(np.shape(grd)[1])




# loop over stationpairs
	n = len(data)
	for i in range(n):

		if data.at[i,'snr'] < snr_min:
			continue
# ToDo: deal with station pairs with several measurements (with different instruments)
# (At the moment, just all added. Probably fine on this large scale)
# find kernel file
		sta1 = data.at[i,'sta1']
		sta2 = data.at[i,'sta2']
	
		if sta1.split('.')[-1][-1] in ['E','N','T','R']:
			msg = "Cannot handle horizontal components"
			raise NotImplementedError(msg)
		if sta2.split('.')[-1][-1] in ['E','N','T','R']:
			msg = "Cannot handle horizontal components"
			raise NotImplementedError(msg)
	
	
# ToDo !!! Replace this by a decent formulation, where the channel is properly set !!! No error for E, R, T, N
		sta1 = "{}.{}..MXZ".format(*sta1.split('.')[0:2])
		sta2 = "{}.{}..MXZ".format(*sta2.split('.')[0:2])
	
		kernelfile = os.path.join(datadir,'kern',"{}--{}.npy".format(sta1,sta2))
		if not os.path.exists(kernelfile):
			print("File does not exist:")
			print(os.path.basename(kernelfile))
			continue


# load kernel
		kernel = np.load(kernelfile)
		if True in np.isnan(kernel):
			print("kernel contains nan, skipping")
			print(os.path.basename(kernelfile))
			continue


# multiply kernel and measurement, add to descent dir. Skip if entry is nan	
		if isnan(data.at[i,'obs']):
			print("No measurement in dataset for:")
			print(os.path.basename(kernelfile))
			continue
		else:
			kernel *= data.at[i,'obs']


		gradient += kernel

# save
		kernelfile = os.path.join(datadir,'grad',os.path.basename(kernelfile))
		np.save(kernelfile, kernel)

		del kernel

# plot
	
	kernelfile = os.path.join(datadir,'grad','grad_all.npy')
	np.save(kernelfile,gradient)
	#plotfile = os.path.join(datadir,'step_'+step,'grad_all.png')

	#plot_grid(grd[0],grd[1],gradient,outfile=plotfile)
	