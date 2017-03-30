import numpy as np
import pandas as pd
import os
import json
from glob import glob
from math import isnan
from noisi.util.plot import plot_grid


def assemble_ascent_dir(source_model,step,snr_min,n_min,save_all=False):




# where is the measurement database located?
	source_config=json.load(open(source_model))
	datadir = os.path.join(source_config['source_path'],'step_' + str(step))
	msrfile = os.path.join(datadir,"{}.measurement.csv".format(source_config['mtype']))



# Read in the csv files of measurement.

	data = pd.read_csv(msrfile)


# allocate the kernel array
	grd = np.load(os.path.join(source_config['project_path'],'sourcegrid.npy'))

	gradient = np.zeros(np.shape(grd)[1])




# loop over stationpairs
	cnt_success = 0
	cnt_lowsnr = 0
	cnt_lown = 0
	cnt_overlap = 0
	cnt_unavail = 0
	n = len(data)
	for i in range(n):

		if data.at[i,'snr'] < snr_min and data.at[i,'snr_a'] < snr_min:
			cnt_lowsnr += 1
			continue

		if data.at[i,'nstack'] < n_min:
			cnt_lown += 1
			continue


# ToDo: deal with station pairs with several measurements (with different instruments)
# (At the moment, just all added. Probably fine on this large scale)
# find kernel file
		sta1 = data.at[i,'sta1']
		sta2 = data.at[i,'sta2']
	
		#if sta1.split('.')[-1][-1] in ['E','N','T','R']:
	#		msg = "Cannot yet handle horizontal components"
#			raise NotImplementedError(msg)
#		if sta2.split('.')[-1][-1] in ['E','N','T','R']:
#			msg = "Cannot yet handle horizontal components"
#			raise NotImplementedError(msg)
	
	
# ToDo !!! Replace this by a decent formulation, where the channel is properly set !!! No error for E, R, T, N
		sta1 = "*.{}..{}".format(sta1.split('.')[1],source_config['channel']) # ignoring network: IRIS has sometimes several network codes at same station
		sta2 = "*.{}..{}".format(sta2.split('.')[1],source_config['channel']) # ignoring network: IRIS has sometimes several network codes at same station
	
		kernelfile1 = os.path.join(datadir,'kern',"{}--{}.npy".format(sta1,sta2))
		kernelfile2 = os.path.join(datadir,'kern',"{}--{}.npy".format(sta2,sta1))
		# Same problem with different network codes.
		# Due to station pairs being in alphabetic order of network.station.loc.cha, different network
		# codes also lead to different ordering.
		try:
			kernelfile = glob(kernelfile1)[0]
		except IndexError:
			try: 
				kernelfile = glob(kernelfile2)[0]
			except IndexError:
				kernelfile = kernelfile1 # this file does not actually exist, but there might be a good reason.
				# Check that first, and then complain.




# Skip if entry is nan: This is most likely due to no measurement taken because station distance too short	
		if isnan(data.at[i,'obs']):
			print("No measurement in dataset for:")
			print(os.path.basename(kernelfile))
			cnt_overlap += 1
			continue

# ...unless somehow the kernel went missing (undesirable case!)

		if not os.path.exists(kernelfile):
			print("File does not exist:")
			print(os.path.basename(kernelfile))
			cnt_unavail += 1
			continue


# load kernel
		kernel = np.load(kernelfile)
		if True in np.isnan(kernel):
			print("kernel contains nan, skipping")
			print(os.path.basename(kernelfile))
			continue


# multiply kernel and measurement, add to descent dir. 
	
		else:
			if kernel.shape()[-1] == 1:
				kernel *= (data.at[i,'syn'] - data.at[i,'obs'])
			else:
				kernel[:,0] *= (data.at[i,'syn'] - data.at[i,'obs'])
				kernel[:,1] *= (data.at[i,'syn_a'] - data.at[i,'obs_a'])
				kernel = kernel[:,0] + kernel[:,1]
			cnt_success += 1 # yuhu


		gradient += kernel

# save
		if save_all:
			kernelfile = os.path.join(datadir,'grad',os.path.basename(kernelfile))
			np.save(kernelfile, kernel)

		del kernel


	
	kernelfile = os.path.join(datadir,'grad','grad_all.npy')
	np.save(kernelfile,gradient)


# output metadata
	outfile	= os.path.join(datadir,'grad','grad_info.txt')
	with open(outfile,'w') as fh:

		# print info
		fh.write('Analyzed %g station pairs of %g successfully.\n' %(cnt_success,n))
		fh.write('No data found for %g station pairs.\n' %cnt_unavail)
		fh.write('No measurement taken for %g station pairs due to short interstation distance.\n' %cnt_overlap) 
		fh.write('Signal to noise ratio below threshold for %g station pairs.\n' %cnt_lowsnr)
		fh.write('Number of staacked windows below threshold for %g station pairs.\n' %cnt_lown)
		fh.write('\nParameters:==============================================================\n')
		fh.write('Source dir: %s \n' %source_model)
		fh.write('Step: %g' %int(step))
		fh.write('Minimum SNR: %g' %snr_min)
		fh.write('Minimum stack length: %g' %int(n_min))
		fh.write('Save all interstation gradients: %s' %str(save_all))
		fh.write('\n=========================================================================\n')
		fh.write('Project:\n')
		# append configurations
		cfg = open(os.path.join(source_config['project_path'],'config.json')).read()
		fh.write(cfg)
		fh.write('\n=========================================================================\n')
		fh.write('Source model:\n')
		fh.write(json.dumps(source_config))
		fh.write('\n=========================================================================\n')
		fh.write('Measurement:\n')
		cfg = open(os.path.join(source_config['source_path'],'measr_config.json')).read()
		fh.write(cfg)
# plot
	#plotfile = os.path.join(datadir,'step_'+step,'grad_all.png')

	#plot_grid(grd[0],grd[1],gradient,outfile=plotfile)
	