import os
import csv
import numpy as np
import matplotlib.pyplot as plt

from cherab.tools.equilibrium import import_eqdsk
from cherab.core.math import Interpolate1DLinear, Interpolate2DLinear

##############################################################################################
##############################################################################################
##############################################################################################

"""
	TODO mmm:

	- add automatic interpolation to give *_f2d
"""

class ASTRASimulation:

	###################

	def __init__(self):

		self._equilibrium = None

		self._Ha_emission_raw = None
		self._Ha_emission_f2d = None

		self._electron_density_raw = None
		self._electron_density_f2d = None

		self._electron_temperature_raw = None
		self._electron_temperature_f2d = None

		self._neutral_density_raw = None
		self._neutral_density_f2d = None

		self._neutral_temperature_raw = None
		self._neutral_temperature_f2d = None

	###################

	@property
	def equilibrium(self):
		return self._equilibrium
	@equilibrium.setter
	def equilibrium(self, value):
		self._equilibrium = value

	###################

	@property
	def Ha_emission_raw(self):
		return self._Ha_emission_raw
	@Ha_emission_raw.setter
	def Ha_emission_raw(self, value):
		self._Ha_emission_raw = value

	###################

	@property
	def Ha_emission_f2d(self):
		return self._Ha_emission_f2d
	@Ha_emission_f2d.setter
	def Ha_emission_f2d(self, value):
		self._Ha_emission_f2d = value

	###################

	@property
	def electron_density_f2d(self):
		return self._electron_density_f2d
	@electron_density_f2d.setter
	def electron_density_f2d(self, value):
		self._electron_density_f2d = value

	###################

	@property
	def electron_temperature_f2d(self):
		return self._electron_temperature_f2d
	@electron_temperature_f2d.setter
	def electron_temperature_f2d(self, value):
		self._electron_temperature_f2d = value

	###################

	@property
	def neutral_density_f2d(self):
		return self._neutral_density_f2d
	@neutral_density_f2d.setter
	def neutral_density_f2d(self, value):
		self._neutral_density_f2d = value

	###################

	@property
	def neutral_temperature_f2d(self):
		return self._neutral_temperature_f2d
	@neutral_temperature_f2d.setter
	def neutral_temperature_f2d(self, value):
		self._neutral_temperature_f2d = value

	###################

	def __getstate__(self):
		state = {
			'equilibrium': self._equilibrium,
			'Ha_emission_raw': self._Ha_emission_raw,
			'Ha_emission_f2d': self._Ha_emission_f2d,
			'electron_density_raw': self._electron_density_raw,
			'electron_temperature_raw': self._electron_temperature_raw,
			'neutral_density_raw': self._neutral_density_raw,
			'neutral_temperature_raw': self._neutral_temperature_raw
		}
		return state

	###################

	def __setstate__(self, state):
		self._equilibrium = state['equilibrium']
		self._Ha_emission_raw = state['Ha_emission_raw']
		self._Ha_emission_f2d = state['Ha_emission_f2d']
		self._electron_density_f2d = state['electron_density_f2d']
		self._electron_temperature_f2d = state['electron_temperature_f2d']
		self._neutral_density_f2d = state['neutral_density_f2d']
		self._neutral_temperature_f2d = state['neutral_temperature_f2d']