from __future__ import print_function
import numpy as np
import os
from glob import glob
from obspy import read
from noisi.util.windows import my_centered
from obspy.signal.interpolation import plot_lanczos_windows


#ToDo pass in all the now default options
def run_preprocess_data(source,bandpass=None,decimator=None,Fs_new=None,overwrite=False,fmt='sac'):

	datalst = os.path.join(source,'observed_correlations','*.'+fmt.lower())
	data = glob(datalst)
	datalst = os.path.join(source,'observed_correlations','*.'+fmt.upper())
	data.extend(glob(datalst))

	if data == []:
		print('No data found.')
		return()

	if not overwrite:
		outdir = os.path.join(source,'processed_correlations')
		os.mkdir(outdir)
	else:
		outdir = os.path.join(source,'observed_correlations')


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

			if Fs_new < tr.stats.sampling_rate:
				print('HAVE YOU filtered?')
			plot_lanczos_windows(a=40,filename='lanczos_response.eps')
			tr.interpolate(Fs_new, method='lanczos',a=40)

		if tr.stats.npts % 2 == 0:
			tr.data = my_centered(tr.data,tr.stats.npts-1)

		tr.write(os.path.join(outdir,os.path.basename(f)),format=fmt)
