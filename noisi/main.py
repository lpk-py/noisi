from __future__ import print_function
import sys
import os
import io
import click
import json
import time

from noisi.scripts.source_grid import setup_sourcegrid as setup_sgrid
from noisi.scripts.run_correlation import run_corr
from noisi.util.prepare_sem_input import prepare_specfem_input
from noisi.scripts.run_measurement import run_measurement
from noisi.scripts.run_adjointsrcs import run_adjointsrcs
from noisi.scripts.run_kernel import run_kernel
from noisi.scripts.run_preprocessing import run_preprocessing

@click.group()
def run():
    """
    Main routine for noise correlation modeling and noise source inversion.
    """
    pass
    

@run.command(help='Initialize a new project.')
@click.argument('project_name')
def setup_project(project_name):
    if os.path.exists(project_name):
        click.echo('Project exists already, must give it a new name.')
        exit()   
    os.makedirs(os.path.join(project_name,'observed_correlations'))
    os.mkdir(os.path.join(project_name,'green_c'))
    #os.mkdir(os.path.join(project_name,'green_c_fd'))
    
    from . import _ROOT
    with io.open(os.path.join(_ROOT,'config','config.json'),'rb+') as fh:
        conf = json.loads(fh.read())
        
    conf['date_created'] = time.strftime("%Y.%m.%d")
    conf['project_name'] = project_name
    conf['project_path'] = os.path.abspath(project_name)
    
    with io.open(os.path.join(project_name,'config.json'),'wb') as fh:
        cf = json.dumps(conf,sort_keys=True, indent=4, separators=(",", ": "))
        fh.write(cf)
    
    click.secho("Copied default config.json to project directory, please edit.")
    
@run.command(help='Determine the source grid and get specfem STATIONS file.')
@click.argument('project_path')
def setup_sourcegrid(project_path):
    setup_sgrid(os.path.join(project_path,'config.json'))
    
@run.command(help='Prepare specfem input.')
@click.argument('project_path')
def specfem_input(project_path):
    prepare_specfem_input(os.path.join(project_path,'config.json'))
    
    
@run.command(help='Initialize a new source model.')
@click.argument('source_model')
def setup_source(source_model):
    if os.path.exists(source_model):
        click.echo('Source exists already, must give it a new name.')
        exit()
        
    if not os.path.exists('config.json'):
        click.echo('No config file for project found \
        (detailing e.g. source grid). Run setup_project first.')
        exit()
    
    os.makedirs(os.path.join(source_model,'step_0'))
    for d in ['adjt','grad','corr']:
        os.mkdir(os.path.join(source_model,'step_0',d))

    from . import _ROOT
    
    with io.open(os.path.join(_ROOT,'config','source_config.json'),'r') as fh:
        conf = json.loads(fh.read())
        conf['date_created'] = unicode(time.strftime("%Y.%m.%d"))
        conf['project_name'] = os.path.basename(os.getcwd())
        conf['project_path'] = os.getcwd()
        conf['source_name'] = source_model
        conf['source_path'] = os.path.abspath(source_model)
        
        
    with io.open(os.path.join(source_model,'source_config.json'),'wb') as fh:
        cf = json.dumps(conf,sort_keys=True, indent=4, separators=(",", ": "))
        fh.write(cf)
    
    with io.open(os.path.join(_ROOT,'config','measr_config.json'),'r') as fh:
        conf = json.loads(fh.read())
        conf['date_created'] = unicode(time.strftime("%Y.%m.%d"))
       
    with io.open(os.path.join(source_model,'measr_config.json'),'wb') as fh:
        cf = json.dumps(conf,sort_keys=True, indent=4, separators=(",", ": "))
        fh.write(cf)
    
    from . import _ROOT
    os.system('cp {} {}'.format(os.path.join(_ROOT,'setup_noisesource.ipynb'),
    source_model))  
    os.system('cp {} {}'.format(os.path.join(_ROOT,'setup_noisesource.py'),
    source_model))  
    
    click.secho("Copied default source_config.json and measr_config.json to source model directory, please edit. \
Please run setup_noisesource.ipynb or setup_noisesource.py after editing to create starting model.")
    
@run.command(help='Filter & truncate synthetics before correlation.')
@click.argument('source_model')
def preprocessing(source_model):
    source_model = os.path.join(source_model,'source_config.json')
    source_config = json.load(open(source_model))
    if source_config['preprocess_do']:
        
        
        dir = os.path.join(source_config['project_path'],'wavefield_processed')
        
        try:
            os.mkdir(dir)
        except:
            pass    
        run_preprocessing(source_config)

@run.command(help='Calculate correlations for selected source model.')
@click.argument('source_model')
@click.argument('step')
def correlation(source_model,step):
    source_model = os.path.join(source_model,'source_config.json')
    run_corr(source_model,step)
    
#@run.command(help='Calculate correlations for selected source model,\nstoring the preliminary kernels in frequency domain.')
#@click.argument('source_model')
#def correlation(source_model):
#    source_model = os.path.join(source_model,'source_config.json')
#    run_corr_fd(source_model)
    

@run.command(help='Run measurement and adjoint sources.')
@click.argument('source_model')
# To do: Include a --test option that produces only plots 
# To do: include a get_parameters_options or something, so that there is no extra step necessary in run_measurement
@click.argument('step')
def measurement(source_model,step):
    
    measr_config = os.path.join(source_model,'measr_config.json')
    source_model = os.path.join(source_model,'source_config.json')
    
    run_measurement(source_model,measr_config,int(step))
    run_adjointsrcs(source_model,measr_config,int(step))


@run.command(help='Calculate gradients.')
@click.argument('source_model')
@click.argument('step')

def kernel(source_model,step):
    source_model = os.path.join(source_model,'source_config.json')
    run_kernel(source_model,step)
    
    
#import config.configure as conf
#import config.configure_source as cs
#import config.source_grid as sg
#@run.command(help='Step-by-step configuration of noise modeling project.')
#def configure():
#    out = conf.get_config_from_user()
#    print(out)
#
#@run.command(help='Create starting model for time-space dep. noise sources.')
#def add_sourcemodel():
#    sourcemodel_name = raw_input('Source model name: ')
#    proj_dir = raw_input('Project directory [%s]: ' %os.getcwd())
#    if proj_dir == '': proj_dir = os.getcwd()
#    out = cs.configure_noisesource_model(proj_dir,sourcemodel_name)
#    print(out)
#
#@run.command(help='Set up the 2-D coordinate grid of source locations.')
#def setup_sourcegrid():
#    proj_dir = raw_input('Project directory [%s]: ' %os.getcwd())
#    if proj_dir == '': proj_dir = os.getcwd()
#    format = raw_input('Output format needed [specfem]: ') 
#    if format == '': 
#        format='specfem'
#    
#    config = os.path.join(proj_dir,'config.txt')
#    sg.setup_sourcegrid(config,format)
#

##ToDo: Set this up properly, esp. looking up the directories!
#@run.command(help='Run the forward model: Create a correlation dataset.')
#def get_forward():
#    proj_dir = raw_input('Project directory [%s]: ' %os.getcwd())
#    sourcemodel_name = raw_input('Source model name: ')
#    out = cs.configure_noisesource_model(proj_dir,sourcemodel_name)
#    print(out)
## ToDo: Extend list of misfits    
#@run.command(help='Run forward model, obtain kernel,\
# update.\n Misfit functions are: energy_difference, log_energy_ratio.')
#@click.argument('misfit', type=click.Choice(['energy_difference',\
#'log_energy_ratio']))
#@click.option('--stop_at',type=click.Choice(['forward','measure','kernels']),\
#help='Use this option to interrupt after part of the workflow.')    
#def iteration(misfit,stop_at):
#    print('We\'ll get there!')
