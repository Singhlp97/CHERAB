import os
import sys
import numpy as np
import matplotlib.pyplot as plt

# Author: Matteo Moscheni
# E-mail: matteo.moscheni@tokamakenergy.co.uk
# February 2022

#try:

sys.path.append('/home/lovepreet/CHERAB_files/Fabio/')

	# Tokamak Energy server access required
try:
		from st40_phys_viewer import Get_data
		from solps_python_scripts.plot_solps_TE import plot_vs_ASTRA
except:
		print()
		print('st40_phys_viewer can not be loaded from outside Tokamak Energy!')
		print()

from solps_python_scripts.plot_solps       import plot_1d, plot_2d
from solps_python_scripts.plot_solps       import plot_mesh, plot_wall_loads
from solps_python_scripts.plot_solps       import plot_sv, plot_rates, plot_mfps
from solps_python_scripts.plot_time_traces import plot_time_traces

from solps_python_scripts.read_b2mn          import read_b2mn
from solps_python_scripts.read_b2fgmtry      import read_b2fgmtry
from solps_python_scripts.read_b2fstate      import read_b2fstate
from solps_python_scripts.read_b2fplasmf     import read_b2fplasmf
from solps_python_scripts.read_ft44          import read_ft44
from solps_python_scripts.read_ft46          import read_ft46
from solps_python_scripts.read_triangle_mesh import read_triangle_mesh

from solps_python_scripts.reactions.read_reactions import read_one_reaction, read_all_reactions
from solps_python_scripts.reactions.compute_rates import compute_rates

from solps_python_scripts.extract_vessel_from_eqdsk import extract_vessel_from_eqdsk
from solps_python_scripts.input_dat_to_mesh_extra import input_dat_to_mesh_extra

from solps_python_scripts.utilities.squeue  import squeue
from solps_python_scripts.utilities.lq      import compute_lq, compute_ln
from solps_python_scripts.utilities.fx      import compute_fx
from solps_python_scripts.utilities.last10s import read_last10s	

#except:

	#print()
	#print('Please, set the proper path in solps_python_scripts/setup.py (line 7) :)')
	#print('If does NOT work, de-comment try/except to see actual error...')
	#print()
