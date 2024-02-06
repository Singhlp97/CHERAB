#!/bin/env python

import os
import json
import argparse
import numpy as np
import pickle
import time

import matplotlib.pyplot as plt

from raysect.core.workflow import MulticoreEngine
from raysect.optical import World, translate, rotate, rotate_z, Point3D, Vector3D, rotate_basis
from raysect.optical import AbsorbingSurface, SpectralFunction
from raysect.optical.material import Conductor, RoughConductor
from raysect.primitive import Mesh
from raysect.optical import InterpolatedSF
    
from raysect.primitive.mesh.vtk import import_vtk, export_vtk
from raysect.primitive.mesh.stl import import_stl

from load_plasma_radiation import load_plasma_radiation

from create_smart_pixels import CreateSmartPixels

###################################################################################################################

def load_mesh(filename = None):
    mesh_type = os.path.splitext(filename)[1]
    if mesh_type == ".vtk":
        return import_vtk(filename)
    elif mesh_type == ".obj":
        return import_obj(filename)
    else:
        raise ValueError("Unrecognised mesh file type '{}'.".format(mesh_file))

###################################################################################################################

def setup_radiation_load(cfg = None, mesh = None, world = None):

	from raysect.optical.observer import MeshCamera, PowerPipeline1D, MonoAdaptiveSampler1D

	pipeline = PowerPipeline1D()

	frame_sampler = MonoAdaptiveSampler1D(pipeline,
    									  fraction    = cfg['raytracing']['frame_sampler']['fraction'],
    									  ratio       = cfg['raytracing']['frame_sampler']['ratio'],
    									  min_samples = cfg['raytracing']['frame_sampler']['min_samples'],
    									  cutoff      = cfg['raytracing']['frame_sampler']['cutoff'])

	camera = MeshCamera(mesh,
    					parent         = world,
    					pipelines      = [pipeline],
    					frame_sampler  = frame_sampler,
    					surface_offset = cfg['raytracing']['observer']['surface_offset'])

	return (camera, pipeline)

###################################################################################################################

def setup_h_alpha_camera(cfg = None, world = None): #mesh = None, 

	from raysect.optical.observer import PowerPipeline2D, MonoAdaptiveSampler2D
	if   cfg['raytracing']['observer']['specs']['type'] == "PinholeCamera":
		from raysect.optical.observer import PinholeCamera
	elif cfg['raytracing']['observer']['specs']['type'] == "VectorCamera":
		from raysect.optical.observer import VectorCamera

	#########

	offset = cfg['raytracing']['observer']['surface_offset'] 

	# CAUTION.
	#
	# shift_direction with z coordinate = 0 is OK only IF camera image plane is vertical!
	
	# MMM 2022.05.09

	if cfg['raytracing']['observer']['specs']['calcam']['use_calcam_calibration'] is True:

		directory = cfg['raytracing']['observer']['specs']['calcam']['calcam_directory']

		pixel_pupil = pickle.load(open(os.path.join(directory, 'pixel_pupil.pkl'), 'rb'))
		pixel_directions = pickle.load(open(os.path.join(directory, 'pixel_directions.pkl'), 'rb'))

		observation_point = Point3D(pixel_pupil[0], pixel_pupil[1], pixel_pupil[2])
		observation_line = Vector3D(0,0,0)

		# observation_line = arithmetic_average{pixel_directions}
		for i in range(pixel_directions.shape[0]):
			for j in range(pixel_directions.shape[1]):
				observation_line += pixel_directions[i,j]

		observation_line = (observation_line / (pixel_directions.shape[0] * pixel_directions.shape[1])).normalise()
		observation_point += observation_line * offset

	else:

		pixel_origins    = None
		pixel_directions = None

		#triangles = mesh.data.triangles
		#vertices = mesh.data.vertices

		#baricentre = Vector3D(vertices[:,0].mean(), vertices[:,1].mean()*0.0, vertices[:,2].mean())
		#shift_direction = Vector3D(baricentre.x, baricentre.y, 0).normalise()
		observation_point = Vector3D(-0.49223175, 1.19159549, 0.10287001) #baricentre - shift_direction * offset
		observation_line = Vector3D(+0.41332112712126, -0.9104318228272708, -0.016719505334317833)

	#########

	pipeline = PowerPipeline2D(display_progress = False)

	frame_sampler = MonoAdaptiveSampler2D(pipeline,
								          fraction    = cfg['raytracing']['frame_sampler']['fraction'],
								          ratio       = cfg['raytracing']['frame_sampler']['ratio'],
								          min_samples = cfg['raytracing']['frame_sampler']['min_samples'],
								          cutoff      = cfg['raytracing']['frame_sampler']['cutoff'])

	if cfg['raytracing']['observer']['specs']['type'] == "PinholeCamera":

		pixels = (cfg['raytracing']['observer']['specs']['nx_pixels'],
			  	  cfg['raytracing']['observer']['specs']['ny_pixels'])

		camera = PinholeCamera(pixels,
							   parent        = world,
							   pipelines     = [pipeline],
							   frame_sampler = frame_sampler,
						       fov           = cfg['raytracing']['observer']['specs']['fov'])
		
		camera.transform = translate(observation_point.x, observation_point.y, observation_point.z) * \
		 				   rotate_basis(observation_line, observation_line.cross(Vector3D(0, 0, 1)))

	elif cfg['raytracing']['observer']['specs']['type'] == "VectorCamera":

		if cfg['raytracing']['observer']['specs']['smart_pixelling']['use_smart_pixelling'] is True:

			(pixel_origins, pixel_directions) = CreateSmartPixels(cfg = cfg,
																  observation_line  = observation_line,
																  observation_point = observation_point)
		
		elif cfg['raytracing']['observer']['specs']['calcam']['use_calcam_calibration'] is True:

			nx = pixel_directions.shape[0]
			ny = pixel_directions.shape[1]
			pixel_origins = np.ndarray((nx, ny), dtype = object)
			print(pixel_origins.shape)
			pixel_pupil = Point3D(pixel_pupil[0], pixel_pupil[1], pixel_pupil[2])

			print(nx, ny)
			# from pupil to pixel_origins on a plane shifted towards source of tiny amount "offset"
			for i in range(nx):
				for j in range(ny):
					pixel_direction = pixel_directions[i,j].normalise()
					shift = offset * np.abs(observation_line.dot(pixel_direction))
					pixel_origins[i,j] = pixel_pupil + pixel_direction * shift

		camera = VectorCamera(parent           = world,
							  pixel_origins    = pixel_origins,
							  pixel_directions = pixel_directions,
							  pipelines        = [pipeline],
							  frame_sampler    = frame_sampler)

	return (camera, pipeline, pixel_origins)

###################################################################################################################
###################################################################################################################
###################################################################################################################

parser = argparse.ArgumentParser('Simulate Halpha camera view.')
parser.add_argument('configFile', type = str, help = "A JSON format task configuration file.")
args = parser.parse_args()

with open(args.configFile, 'r') as f:
    cfg = json.load(f)

# input

run_directory   = cfg['run']
input_directory = os.path.join(cfg['baserun'], 'input', run_directory)
print(input_directory)
wall_directory  = os.path.join(input_directory, cfg['mesh']['wall_directory'])
observer_directory  = os.path.join(input_directory, cfg['mesh']['observer_directory']) 

try: reflecting_directory = os.path.join(input_directory, cfg['mesh']['reflecting_directory']) 
except: pass

# mesh

periodicity = cfg['mesh']['periodicity']
periodicity_limiters = cfg['mesh']['periodicity_limiters']

wall_mesh_files = os.listdir(wall_directory)
observer_mesh_files = os.listdir(observer_directory)
try: reflecting_mesh_files = os.listdir(reflecting_directory)
except: pass

# output

type_radiation = cfg['plasma']['type_radiation']

use_step   = str(cfg['raytracing']['sampling']['use_step_function'])
use_abs    = str(cfg['plasma']['absorption']['use_absorption_function'])
use_scat   = str(cfg['plasma']['scattering']['use_scattering_function'])

if cfg['plasma']['homemade']['use_homemade_emission'] is False:
	use_astra  = str(cfg['plasma']['ASTRA']['use_ASTRA_emission'])
	use_b2     = str(cfg['plasma']['SOLPS']['use_B2_emission'])
	use_eirene = str(cfg['plasma']['SOLPS']['use_EIRENE_emission'])
	use_extra  = str(cfg['plasma']['SOLPS']['use_extra_emission'])
	details = ["step="   + use_step,
			   "abs="    + use_abs, 
			   "scat="   + use_scat, 
			   "astra="  + use_astra,
			   "b2="     + use_b2,
			   "eirene=" + use_eirene,
			   "extra="  + use_extra,
			   cfg['plasma']['Mapper']]
else:
	use_homemade = str(True)
	details = ["step="   + use_step,
			   "abs="    + use_abs,
			   "scat="   + use_scat,
			   "homemade="  + use_homemade,
			   cfg['plasma']['Mapper']]

details = ".".join(details)

# - check if any existing directory
# - assign serial number accordingly
# - create output directory 
# - copy configFile.json in output directory to have run info available

output_directory = os.path.join(cfg['baserun'], 'output', run_directory, type_radiation)
os.system("mkdir -p " + output_directory)
cases = os.listdir(output_directory)

num_case = 0

for case in cases:

	in_case = os.listdir(os.path.join(output_directory, case))

	# rm directory if empty
	if in_case == [] or in_case == ['input']:
		os.system("rm -rf " + os.path.join(output_directory, case))

	# skip num_case index at beginning of string
	elif ".".join(case.split(".")[1:]) == details:
		if int(case.split(".")[0]) > num_case:
			num_case = int(case.split(".")[0])

details = ".".join([str(num_case + 1), details])

output_directory = os.path.join(output_directory, details)

os.system("mkdir -p " + os.path.join(output_directory, "input/"))
os.system("cp configFile.json " + os.path.join(output_directory, "input/"))
if cfg['plotting']['save_figures'] is True:
	cfg['output_directory_extended'] = os.path.join(output_directory, "input/")

###################################
# OVERALL RAD_FUNCTION GENERATION #
###################################

world = World()

plasma = load_plasma_radiation(cfg = cfg, parent = world)

pkl_dir = os.path.join(cfg['baserun'],
                       'input',
                       cfg['run'],
                       cfg['plasma']['SOLPS']['SOLPS_directory'])
os.system("rm -rf " + pkl_dir + "/*.pkl")

# if images forgotten open, closed automatically after 100s for simulation to proceed
plt.show()

# ciao

# #############################
# # CREATION OFABSORBING WALL #
# #############################

for wall_mesh_file in wall_mesh_files:

	try:
		rescale = cfg['raytracing']['sampling']['resize_speedup_scaling']
	except:
		rescale = 1E+00

	try:
		wall_mesh = import_vtk(os.path.join(wall_directory, wall_mesh_file), scaling = rescale)
		for i in range(int(360/periodicity)):
			wall_mesh.instance(parent = world,
		                  	   material = AbsorbingSurface(),
		                  	   transform = rotate_z(i * periodicity))

	except:
		wall_mesh = import_stl(os.path.join(wall_directory, wall_mesh_file), scaling = 1E-03)
		for i in range(int(360/periodicity)):
			wall_mesh.instance(parent = world,
		                  	   material = AbsorbingSurface(),
		                  	   transform = rotate_z(i * periodicity))

# ###############################
# # CREATION OF REFLECTING WALL #
# ###############################

# print('Importing refecting mesh...\n')

# wavelength = np.array([656])
# index_C = InterpolatedSF(wavelength, np.array([2.8679]))  # at 656 nm
# extinction_C = InterpolatedSF(wavelength, np.array([1.6197]))  # at 656 nm 1.6197
# index_Fe = InterpolatedSF(wavelength, np.array([2.9171]))  # at 656 nm
# extinction_Fe = InterpolatedSF(wavelength, np.array([3.0964]))  # at 656 nm

# for reflecting_mesh_file in reflecting_mesh_files:

# 	try:
# 		rescale = cfg['plasma']['sampling']['resize_speedup_scaling']
# 	except:
# 		rescale = 1E+00

# 	reflecting_mesh = import_vtk(reflecting_directory + reflecting_mesh_file,
# 								 scaling = rescale)
# 	reflecting_mesh.instance(parent = world,
# 	                         material = AbsorbingSurface())    

#########################
# OBSERVATION PROCEDURE #
#########################

#############################################

for observer_mesh_file in observer_mesh_files:

	output_basename = os.path.splitext(observer_mesh_file)[0]

	mesh = import_vtk(observer_directory + observer_mesh_file)

	#############################################

#	if cfg['baserun'] == "h_alpha_camera":
#		(camera, pipeline, pixel_origins) = setup_h_alpha_camera(cfg = cfg, world = world)#, mesh = mesh)
#	elif cfg['baserun'] == "radiation_load":
	(camera, pipeline) = setup_radiation_load(cfg = cfg, world = world, mesh = mesh)

	camera.ray_max_depth  = cfg['raytracing']['max_ray_depth']
	camera.spectral_bins  = cfg['raytracing']['observer']['spectral_bins']
	camera.min_wavelength = cfg['raytracing']['observer']['min_wavelength']
	camera.max_wavelength = cfg['raytracing']['observer']['max_wavelength']

	camera.pixel_samples  = cfg['raytracing']['observer']['pixel_samples']
	camera.render_engine  = MulticoreEngine(cfg['raytracing']['number_CPUs'])

	#ciao

    #############################################

	render_pass = 0

	while (not camera.render_complete) and (render_pass < cfg['raytracing']['max_render_passes']):

		render_pass += 1
		print('Render pass {}:'.format(render_pass))
		
		camera.observe()		

		from write_results_radiation_load import write_results
		write_results(cfg = cfg, name = output_basename, camera = camera, pipeline = pipeline)

