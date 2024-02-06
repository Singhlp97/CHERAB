import numpy as np
import matplotlib.pyplot as plt

from accessories import SoundSpeed

from solps_python_scripts.utilities.last10s import read_last10s

mp  = 1.67E-27 # [kg]
mD  = 2 * mp
mD2 = 2 * mD

def compute_neutrals():

	where = '/ancalagon_data/physics/matteo.moscheni/solps-iter/runs/playAround_9229/equ@65ms_20220513/Freeze.20220512/'

	last10s = read_last10s(where = where)

	tw3dr = 1E-01 # [eV]: wall temperature

	ti3dr = last10s['ti3dr'][:,1]
	te3dr = last10s['te3dr'][:,1]
	ni3dr = last10s['ne3dr'][:,1]
	fn3dr = last10s['fn3dr'][:,1]

	sint3dr = last10s['ft3dr'][:,1] / last10s['ft3drP'][:,1]

	Pfast  = 9.5E-01
	Ptherm = 1.0E+00 - Pfast

	# tg3dr = Pfast * ti3dr + Ptherm * tw3dr
	tg3dr = 3E+00 * np.ones(ti3dr.shape)
	mg    = Pfast * mD    + Ptherm * mD2

	# 3 eV = frank-condon dissociation energy if D2 break-up occurssssssssssssssss
	# is it a dominant process? or just fortuitous

	ng3dr = ni3dr * sint3dr * np.sqrt(ti3dr / mD * mg / tg3dr)
	# ng3dr_fn = fn3dr / np.sqrt(tg3dr * 1.6E-19 / mg)

	plt.figure()
	plt.semilogy(last10s['dsr'], ng3dr, label = 'ng3dr@3eV')
	# plt.semilogy(last10s['dsr'], ng3dr / sint3dr, label = 'ng3dr / sint3dr')
	# plt.semilogy(last10s['dsr'], ng3dr_fn, label = 'ng3dr from fn3dr')
	plt.semilogy(last10s['dsr'], ni3dr, label = 'ni3dr')
	plt.semilogy(last10s['dnb23dr'][:,0], last10s['dnb23dr'][:,1], label = 'dnb23dr')
	plt.legend(loc = 'best')
	plt.show()

