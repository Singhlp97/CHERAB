import os
import csv
import numpy as np
import matplotlib.pyplot as plt

from raysect.core.math.function.float import Discrete2DMesh

##############################################################################################
##############################################################################################
##############################################################################################

"""
	TODO mmm:

	- create DiscreteMesh2D (~ KdTree) once and for all
	  and then use instance() method
	- put Discrete2D as above but in a separate method, NOT repeat N times
	- add properties in SOLPSSimulation.eirene
	- add velocities (see read_ft46.py)
"""

class EIRENESimulation:

	###################

	def __init__(self, mesh):

		self._mesh = mesh

		self._initial_setup()

	def _initial_setup(self):

		# mesh quantities

		self._vertex_coords = self._mesh['nodes']
		self._triangles = self._mesh['cells']

		# densities [m^{-3}]

		self._pdena = None
		self._pdena_f2d = None

		self._pdenm = None
		self._pdenm_f2d = None

		self._pdeni = None
		self._pdeni_f2d = None

		# energy densities ~ pressures [Pa]

		self._edena = None
		self._edena_f2d = None

		self._edenm = None
		self._edenm_f2d = None

		self._edeni = None
		self._edeni_f2d = None

		# temperatures [eV]

		self._tdena = None
		self._tdena_f2d = None

		self._tdenm = None
		self._tdenm_f2d = None

		self._tdeni = None
		self._tdeni_f2d = None

		# Halpha emission

		# [s^{-1}]: Einstein coefficient for H*(3)->H*(2) transition
		self._Ha_Einstein_coeff = 4.41E+07
		self._Ha_emission = None
		self._Ha_emission_f2d = None

		# excited state fractional abundance [-]

		self._fexct = None
		self._fexct_2d = None

		# absorption macroscopic cross-section [m^{-1}]

		self._SigmaAbs = None
		self._SigmaAbs_f2d = None

	###################

	def _clean(self, value):
		value[np.isnan(value)] = 0.0
		return value

	###################

	@property
	def mesh(self):
		return self._mesh

	###################

	@property
	def vertex_coords(self):
		return self._vertex_coords

	###################

	@property
	def triangles(self):
		return self._triangles

	###################

	@property
	def pdena(self):
		return self._pdena

	@property
	def pdena_f2d(self):
		return self._pdena_f2d

	@pdena.setter
	def pdena(self, value):
		value = np.array(value, dtype = np.float64, copy = False)
		self._pdena = self._clean(value)
		self._pdena_f2d = Discrete2DMesh(vertex_coords = self._vertex_coords,
										 triangles = self._triangles,
										 triangle_data = np.reshape(self._pdena, (self._mesh['cells'].shape[0])),
										 limit = False,
										 default_value = 0.0)	

	###################

	@property
	def pdenm(self):
		return self._pdenm

	@property
	def pdenm_f2d(self):
		return self._pdenm_f2d

	@pdenm.setter
	def pdenm(self, value):
		value = np.array(value, dtype = np.float64, copy = False)
		self._pdenm = self._clean(value)
		self._pdenm_f2d = Discrete2DMesh(vertex_coords = self._vertex_coords,
										 triangles = self._triangles,
										 triangle_data = np.reshape(self._pdenm, (self._mesh['cells'].shape[0])),
										 limit = False,
										 default_value = 0.0)

	###################

	@property
	def pdeni(self):
		return self._pdeni

	@property
	def pdeni_f2d(self):
		return self._pdeni_f2d

	@pdeni.setter
	def pdeni(self, value):
		value = np.array(value, dtype = np.float64, copy = False)
		self._pdeni = self._clean(value)
		self._pdeni_f2d = Discrete2DMesh(vertex_coords = self._vertex_coords,
										 triangles = self._triangles,
										 triangle_data = np.reshape(self._pdeni, (self._mesh['cells'].shape[0])),
										 limit = False,
										 default_value = 0.0)	

	###################

	@property
	def edena(self):
		return self._edena

	@property
	def edena_f2d(self):
		return self._edena_f2d

	@edena.setter
	def edena(self, value):
		value = np.array(value, dtype = np.float64, copy = False)
		self._edena = self._clean(value)
		self._edena_f2d = Discrete2DMesh(vertex_coords = self._vertex_coords,
										 triangles = self._triangles,
										 triangle_data = np.reshape(self._edena, (self._mesh['cells'].shape[0])),
										 limit = False,
										 default_value = 0.0)	

	###################

	@property
	def edenm(self):
		return self._edenm

	@property
	def edenm_f2d(self):
		return self._edenm_f2d

	@edenm.setter
	def edenm(self, value):
		value = np.array(value, dtype = np.float64, copy = False)
		self._edenm = self._clean(value)
		self._edenm_f2d = Discrete2DMesh(vertex_coords = self._vertex_coords,
										 triangles = self._triangles,
										 triangle_data = np.reshape(self._edenm, (self._mesh['cells'].shape[0])),
										 limit = False,
										 default_value = 0.0)

	###################

	@property
	def edeni(self):
		return self._edeni

	@property
	def edeni_f2d(self):
		return self._edeni_f2d

	@edeni.setter
	def edeni(self, value):
		value = np.array(value, dtype = np.float64, copy = False)
		self._edeni = self._clean(value)
		self._edeni_f2d = Discrete2DMesh(vertex_coords = self._vertex_coords,
										 triangles = self._triangles,
										 triangle_data = np.reshape(self._edeni, (self._mesh['cells'].shape[0])),
										 limit = False,
										 default_value = 0.0)	

	###################

	@property
	def tdena(self):
		return self._tdena

	@property
	def tdena_f2d(self):
		return self._tdena_f2d

	@tdena.setter
	def tdena(self, value):
		value = np.array(value, dtype = np.float64, copy = False)
		self._tdena = self._clean(value)
		self._tdena_f2d = Discrete2DMesh(vertex_coords = self._vertex_coords,
										 triangles = self._triangles,
										 triangle_data = np.reshape(self._tdena, (self._mesh['cells'].shape[0])),
										 limit = False,
										 default_value = 0.0)	

	###################

	@property
	def tdenm(self):
		return self._tdenm

	@property
	def tdenm_f2d(self):
		return self._tdenm_f2d

	@tdenm.setter
	def tdenm(self, value):
		value = np.array(value, dtype = np.float64, copy = False)
		self._tdenm = self._clean(value)
		self._tdenm_f2d = Discrete2DMesh(vertex_coords = self._vertex_coords,
										 triangles = self._triangles,
										 triangle_data = np.reshape(self._tdenm, (self._mesh['cells'].shape[0])),
										 limit = False,
										 default_value = 0.0)	

	###################

	@property
	def tdeni(self):
		return self._tdeni

	@property
	def tdeni_f2d(self):
		return self._tdeni_f2d

	@tdeni.setter
	def tdeni(self, value):
		value = np.array(value, dtype = np.float64, copy = False)
		self._tdeni = self._clean(value)
		self._tdeni_f2d = Discrete2DMesh(vertex_coords = self._vertex_coords,
										 triangles = self._triangles,
										 triangle_data = np.reshape(self._tdeni, (self._mesh['cells'].shape[0])),
										 limit = False,
										 default_value = 0.0)	

	###################

	@property
	def Ha_Einstein_coeff(self):
		return self._Ha_Einstein_coeff

	@property
	def Ha_emission(self):
		return self._Ha_emission

	@property
	def Ha_emission_f2d(self):
		return self._Ha_emission_f2d

	@Ha_emission.setter
	def Ha_emission(self, value):
		value = np.array(value, dtype = np.float64, copy = False)
		self._Ha_emission = self._clean(value)
		self._Ha_emission_f2d = Discrete2DMesh(vertex_coords = self._vertex_coords,
										       triangles = self._triangles,
										 	   triangle_data = np.reshape(self._Ha_emission, \
										 			            (self._mesh['cells'].shape[0])),
										 	   limit = False,
										 	   default_value = 0.0)

	###################

	@property
	def f_exct(self):
		return self._f_exct

	@property
	def f_exct_f2d(self):
		return self._f_exct_f2d

	@f_exct.setter
	def f_exct(self, value):
		value = np.array(value, dtype = np.float64, copy = False)
		self._f_exct = self._clean(value)
		self._f_exct_f2d = Discrete2DMesh(vertex_coords = self._vertex_coords,
										  triangles = self._triangles,
										  triangle_data = np.reshape(self._f_exct, \
										  				  (self._mesh['cells'].shape[0])),
										  limit = False,
										  default_value = 0.0)

	###################

	@property
	def SigmaAbs(self):
		return self._SigmaAbs

	@property
	def SigmaAbs_f2d(self):
		return self._SigmaAbs_f2d

	@SigmaAbs.setter
	def SigmaAbs(self, value):
		value = np.array(value, dtype = np.float64, copy = False)
		self._SigmaAbs = self._clean(value)
		self._SigmaAbs_f2d = Discrete2DMesh(vertex_coords = self._vertex_coords,
										    triangles = self._triangles,
										    triangle_data = np.reshape(self._SigmaAbs, \
										  				    (self._mesh['cells'].shape[0])),
										    limit = False,
										    default_value = 0.0)	

	###################

	def __getstate__(self):
		state = {
			'mesh': self._mesh,
			'vertex_coords': self._vertex_coords,
			'triangles': self._triangles,
			'pdena': self._pdena,
			'pdenm': self._pdenm,
			'pdeni': self._pdeni,
			'edena': self._edena,
			'edenm': self._edenm,
			'edeni': self._edeni,
			'tdena': self._tdena,
			'tdenm': self._tdenm,
			'tdeni': self._tdeni,
			'Ha_emission': self._Ha_emission
		}
		return state

	###################

	def __setstate__(self, state):
		self._mesh = state['mesh']
		self._vertex_coords = state['vertex_coords']
		self._triangles = state['triangles']
		self._pdena = state['pdena']
		self._pdenm = state['pdenm']
		self._pdeni = state['pdeni']
		self._edena = state['edena']
		self._edenm = state['edenm']
		self._edeni = state['edeni']
		self._tdena = state['tdena']
		self._tdenm = state['tdenm']
		self._tdeni = state['tdeni']
		self._Ha_emission = state['Ha_emission']