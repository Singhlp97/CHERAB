import os
import csv
import math
import numpy as np
import matplotlib.pyplot as plt

from decimal import Decimal

from raysect.core import Point2D, Point3D, translate, rotate_z
from raysect.optical.material import VolumeTransform
from raysect.primitive import Cylinder, Subtract

from cherab.tools.equilibrium import import_eqdsk
from cherab.core.math import Interpolate1DLinear, Interpolate2DLinear
from cherab.core.math.mappers import AxisymmetricMapper, DiscreteToroidalMapper
from cherab.tools.emitters import RadiationFunction
from cherab.core.math import PolygonMask2D

from quantum_data import FractionExcitedHydrogenGas

from create_step_function import CreateStepFunction
from create_absorption_function import CreateAbsorptionFunction

from solps_python_scripts.reactions.read_reactions import read_one_reaction
from solps_python_scripts.reactions.compute_rates  import evaluate_rate

from homemade_distributions import make_homemade_emission
from accessories import plot_1d, plot_2d, make_chamber_wall

##############################################################################################
##############################################################################################
##############################################################################################

def make_homemade_emitter(cfg = None, parent = None):

	"""

	MODIFYYYYYYYYYYY

	Non-spectral emitter with the emissivity defined as SOLPSFunction2D.
	:param SOLPSMesh solps_mesh: SOLPS simulation mesh.
	:param SOLPSFunction2D radiation_function: Emissivity in W m-3.
	:param Node parent: parent node in the scenegraph, e.g. a World object.
	:param float step: Volume integration step in meters.
	:rtype: Primitive

	TODO mmm:
	- Check whether more efficient to have N 3D emitters or to have N 2D emitters
	  then combined in 1 overall 3D emitter
	- Add possibility to choose among different type of radiation emission when
	  Axisymmetric source is used
	"""

	RAD2DEG = 180 / np.pi

	(_, _, (Da_f2d, _)) = make_homemade_emission(cfg = cfg, verbose = True)

	location = os.path.join(cfg['baserun'],
							'input',
							cfg['run'],
							cfg['plasma']['ASTRA']['ASTRA_directory'])

	eqFile = cfg['plasma']['ASTRA']['geqdsk_name']
	eqTime = cfg['plasma']['ASTRA']['geqdsk_time']

	equilibrium = import_eqdsk(os.path.join(location, eqFile + '_' + eqTime + '.geqdsk'))

	#################################################
	################## ABSORPTION ###################
	#################################################

	use_absorption_function = cfg['plasma']['absorption']['use_absorption_function']

	if use_absorption_function is True:
	    (absorption_function_2d, sn_max) = CreateAbsorptionFunction(cfg = cfg, solps_simulation = None)
	    absorption_function_3d = 0.0
	else:
	    absorption_function_2d = None
	    absorption_function_3d = None
	    sn_max = math.inf

	#################################################
	################## SCATTERING ###################
	#################################################

	use_scattering_function = cfg['plasma']['scattering']['use_scattering_function']
	collisions_max = cfg['plasma']['scattering']['collisions_max']

	if use_scattering_function is True:
	    (scattering_function_2d, sn_max) = CreateAbsorptionFunction(cfg = cfg, solps_simulation = None) # CreateScatteringFunction to be createddddddd       !!!!!
	    scattering_function_3d = 0.0
	else:
	    scattering_function_2d = None
	    scattering_function_3d = None
	    sn_max = math.inf

	#####################################################
	################## TYPE OF MAPPER ###################
	#####################################################

	if cfg['plasma']['Mapper'] == "DiscreteToroidalMapper":
	    
	    periodicity = 2 * np.pi / int(cfg['limiters']['number'])
	    where_non_zero = cfg['limiters']['angular_width_deg'] / RAD2DEG

	    radiation_function_3d = DiscreteToroidalMapper(Da_f2d,
	                                                   periodicity = periodicity,
	                                                   where_non_zero = where_non_zero)

	    if use_absorption_function is True:
	    	absorption_function_3d = DiscreteToroidalMapper(absorption_function_2d,
	    													periodicity = periodicity,
	    													where_non_zero = where_non_zero)

	    if use_scattering_function is True:
	    	scattering_function_3d = DiscreteToroidalMapper(scattering_function_2d,
	    													periodicity = periodicity,
	    													where_non_zero = where_non_zero)

	elif cfg['plasma']['Mapper'] == "AxisymmetricMapper":
	    
	    radiation_function_3d = AxisymmetricMapper(Da_f2d)

	    if use_absorption_function is True:
	    	absorption_function_3d = AxisymmetricMapper(absorption_function_2d)

	#######################################################
	################## TYPE OF SAMPLING ###################
	#######################################################

	use_step_function = cfg['raytracing']['sampling']['use_step_function']

	if use_step_function is True:
	    step_function_3d = CreateStepFunction(cfg = cfg, solps_simulation = None)
	else:
	    step_function_3d = None

	#########################################################
	################## RADIATIVE EMISSION ###################
	#########################################################

	emitter = RadiationFunction(radiation_function =      radiation_function_3d,

	                            use_step_function =       use_step_function,
	                            step_function_3d =        step_function_3d,
	                            step_max =                cfg['raytracing']['sampling']['step_max'],

	                            use_scattering_function = use_scattering_function, # test
	                            scattering_function_3d =  scattering_function_3d,  # test: scattering
	                            collisions_max =          collisions_max,          # test: 100 by default

	                            use_absorption_function = use_absorption_function,
	                            absorption_function_3d =  absorption_function_3d,
	                          #  sn_max =                  sn_max,

	                            step = cfg['raytracing']['sampling']['step_uniform'])

	material = VolumeTransform(emitter, transform = translate(0, 0, 0))

	##########################################################
	################## ENCLOSING PRIMITIVE ###################
	##########################################################

	padding = 1E-04

	max_r = np.max([equilibrium.r_data.max(), 0.0]) + padding
	min_r = equilibrium.r_data.min() - padding
	max_z = equilibrium.z_data.max() + padding
	min_z = equilibrium.z_data.min() - padding

	# hollow cylinder

	height = max_z - min_z
	hollow_cylinder = Subtract(Cylinder(max_r, height),
	                           Cylinder(min_r, height),
	                           transform = translate(0, 0, +min_z))

	######################################################
	################## CREATING SOURCE ###################
	######################################################   

	if cfg['plasma']['Mapper'] == "DiscreteToroidalMapper":
	    align = - 0.5 * RAD2DEG * (periodicity - where_non_zero)
	else:
	    align = 0.0

	plasma_volume = hollow_cylinder
	plasma_volume.parent = parent
	if cfg['plasma']['Mapper'] == "DiscreteToroidalMapper":
		material.transform  = translate(0, 0, -min_z) * rotate_z(0.5*periodicity*RAD2DEG-13/7*2.05*where_non_zero*RAD2DEG)# * rotate_z(+30-11.25+align)
	else:
		material.transform  = translate(0, 0, -min_z) * rotate_z(align)
	plasma_volume.material = material

	######################################################
	################## PLOTTING SOURCE ###################
	######################################################

	if cfg['plotting']['plot_total_emission'] is True:

		plot_2d(ri = np.array([min_r, max_r]), zi = np.array([min_z, max_z]), f2d = Da_f2d, scale = "log", vmin = 1E+18, eq = equilibrium, title = r'$D_{\alpha}$ emission $[ph/m^3/s]$', cmap_label = '$log_{10}$')

		if cfg['plotting']['save_figures'] is True:
			if os.path.exists(os.path.join(cfg['output_directory_extended'], 'emission.png')) is True:
				plt.savefig(os.path.join(cfg['output_directory_extended'], 'emission_extra.png'))
		else:
			plt.savefig(os.path.join(cfg['output_directory_extended'], 'emission.png'))

	#######################

	return emitter
