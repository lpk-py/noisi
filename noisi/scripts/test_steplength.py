import numpy as np
import pandas as pd
import os
import json

from glob import glob
from math import isnan
import sys
from noisi.util.corr_pairs import glob_obs_corr
from noisi.my_classes.noisesource import NoiseSource

####################################
# ToDo: more fancy and more secure with click or argparse
source_model = sys.argv[1]
oldstep = sys.argv[2]
min_snr = 0.0
min_stck = 640.
nr_msr = 20
step_length = 10
perc_step_length = None
# Only if the following is set to True, a small subset (nr_msr) 
# of data will be selected and copied and their misfit evaluated
# for a step length test. Otherwise, only the update of the source model
# is performed. 
prepare_test_steplength = True
####################################


def _update_steepestdesc(old_model,
	neg_grad,
	step_length=None,
	perc_step_length=None,
	project=False,
	smooth=False):

	if step_length is not None and perc_step_length is not None:
		raise ValueError('Only one of step_length and perc_step_length can\
			be specified.')

# just in case:
	os.system('cp {} {}'.format(old_model,old_model+'.bak'))
# This model will be overwritten, need to be careful with that
# read the model
	src_model = NoiseSource(old_model)

# smooth the model
	if smooth:
		raise NotImplementedError('Sorry, not implemented yet.')

# project the model
	if project:
		raise NotImplementedError('Sorry, not implemented yet.')
		# if implemented, here should be a projection of the new kernel
		# onto the distr_basis functions, thereby yielding the new distr_weights

	else:
	# the spectrum remains unchanged.
	# assuming there is one basis only, this will be updated with the new kernel
		if step_length is not None:	
			src_model.model['distr_basis'][:] += neg_grad * step_length
		elif perc_step_length is not None:
			src_model.model['distr_basis'][:] += neg_grad/np.max(np.abs(neg_grad)) * perc_step_length

# write to file
# close the underlying h5 file	
	src_model.model.close()

	return()



############ Preparation procedure #################################################

# where is the measurement database located?
source_model = os.path.join(source_model,'source_config.json')
source_config=json.load(open(source_model))
datadir = os.path.join(source_config['source_path'],'step_' + str(oldstep))
msrfile = os.path.join(datadir,"{}.measurement.csv".format(source_config['mtype']))

# Initialize the new step directory
newstep = int(oldstep) + 1
newdir = os.path.join(source_config['source_path'],'step_' + str(newstep))

if not os.path.exist(newdir):
	newdir = os.path.join(source_config['source_path'],'step_' + str(newstep))
	os.mkdir(newdir)
	os.mkdir(os.path.join(newdir,'obs_slt'))
	os.mkdir(os.path.join(newdir,'corr'))
	os.mkdir(os.path.join(newdir,'adjt'))
	os.mkdir(os.path.join(newdir,'grad'))
	os.mkdir(os.path.join(newdir,'kern'))

os.system('cp {} {}'.format(os.path.join(datadir,'base_model.h5'),newdir))
os.system('cp {} {}'.format(os.path.join(datadir,'starting_model.h5'),newdir))




if prepare_test_steplength:

	# Read in the csv files of measurement.
	data = pd.read_csv(msrfile)
	# Get a set of n randomly chosen station pairs. Criteria: minimum SNR, 
	# ---> prelim_stations.txt
	data_accept = data[(data.snr > min_snr)]
	data_accept = data_accept[(data_accept.nstack > min_stck)]

	data_select = data_accept.sample(n=nr_msr)

	stafile = open(os.path.join(newdir,'stations_slt.txt'),'w')
	stafile.write("Station pairs to be used for step lenght test:\n")

	inffile = open(os.path.join(newdir,'step_length_test_info.txt'),'w')
	inffile.write('Parameters:\n')
	inffile.write('source_model: %s\n' %source_model)
	inffile.write('old step: %s\n' %oldstep)
	inffile.write('min_snr %g\n' %min_snr)
	inffile.write('min_stck %g\n' %min_stck)
	inffile.write('step_length %g\n' %step_length)
	inffile.write('-'*40)
	inffile.write("\nStation pairs to be used for step lenght test:\n")

	cum_misf = 0.0
	# Take care of the test set for the step length test
	
	for i in range(nr_msr):

		sta1 = data.at[i,'sta1']
		sta2 = data.at[i,'sta2']
		
		lat1 = data.at[i,'lat1']
		lat2 = data.at[i,'lat2']
		lon1 = data.at[i,'lon1']
		lon2 = data.at[i,'lon2']

		misf = data.at[i,'l2_norm']
		cum_misf += misf
		# synthetics in the old directory?
		#synth_filename = os.path.join(datadir,'corr','{}--{}.sac'.format(sta1,sta2))
		#print(synth_filename)
		# copy the relevant observed correlation, oh my
		obs_dir = os.path.join(source_config['source_path'],'observed_correlations')
		obs_correlations = glob_obs_corr(sta1,sta2,obs_dir,ignore_network=True)
		
		if len(obs_correlations) > 0:
			
			sta1 = sta1.split('.')
			sta2 = sta2.split('.')
			stafile.write('{} {} {} {}\n'.format(*(sta1[0:2]+[lat1]+[lon1])))
			stafile.write('{} {} {} {}\n'.format(*(sta2[0:2]+[lat2]+[lon2])))

			inffile.write('{} {}, {} {} L2 misfit: {}\n'.format(*(sta1[0:2]+sta2[0:2]+[misf]))
			

			for corrs in obs_correlations:

				os.system('cp {} {}'.format(corrs,os.path.join(newdir,'obs_slt')))
			


	inffile.write('-'*40)
	inffile.write('\nCumulative misfit: %g\n' %cum_misf)
	inffile.write('-'*40)
	inffile.close()
	stafile.close()
else: 
	pass

# Set up a prelim_sourcemodel.h5: 
# Contains starting model + step length * (-grad) for steepest descent
# This would be the point to project to some lovely basis functions..
grad = os.path.join(datadir,'grad','grad_all.npy')
neg_grad = -1. * np.load(grad)
new_sourcemodel = os.path.join(newdir,'starting_model.h5')

_update_steepestdesc(new_sourcemodel,neg_grad,step_length=step_length,
	perc_step_length=perc_step_length,project=False,smooth=False)
# (outside of this script) forward model selected correlations
# (outside of this script) evaluate misfit for selected correlations