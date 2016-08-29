from __future__ import print_function
import numpy as np
import os
from glob import glob
from obspy import read

#ToDo pass in all the now default options
def run_preprocess_data(source,step,bandpass=None,decimator=None,Fs_new=None,overwrite=False,fmt='sac'):

	data = os.path.join(source,'step_'+str(step),'observed_correlations','*.'+fmt)
	data = glob(data)

	if data == []:
		print('No data found.')
		return()

	if not overwrite:
		outdir = os.path.join(source,'step_'+str(step),'processed_correlations')
		os.mkdir(outdir)
	else:
		outdir = os.path.join(source,'step_'+str(step),'observed_correlations')


	for f in data:

		try:
			tr = read(f)[0]
		except:
			print("Could not read file:")
			print(f)	
			continue

		if bandpass is not None:

			# Using zerophase is essential for correlation
			tr.filter('bandpass',
				freqmin = bandpass[0],
				freqmax = bandpass[1],
				corners = bandpass[2],
				zerophase=True)


		if decimator is not None:

			tr.decimate(decimator)


		if Fs_new is not None:

			if Fs_new < tr.stats.Fs:
				print('HAVE YOU filtered?')

			tr.interpolate(Fs_new, method='lanczos')


		tr.write(os.path.join(outdir,f),format=fmt)
