import os
import numpy as np

from cherab.tools.equilibrium import import_eqdsk

from homemade_distributions import make_homemade_emission

################################################################################################

def homemade_to_solps_format(cfg = None):

	location = os.path.join(cfg['baserun'],
							'input',
							cfg['run'],
							cfg['plasma']['ASTRA']['ASTRA_directory'])
	eqFile = cfg['plasma']['ASTRA']['geqdsk_name']
	eqTime = cfg['plasma']['ASTRA']['geqdsk_time']
	equilibrium = import_eqdsk(os.path.join(location, eqFile + '_' + eqTime + '.geqdsk'))

	(_, _, (Da_f2d, _)) = make_homemade_emission(cfg = cfg, verbose = False)

	nx = cfg['plasma']['homemade']['nx']
	ny = cfg['plasma']['homemade']['ny']
	x_range = equilibrium.r_range
	y_range = equilibrium.z_range
	x = np.linspace(x_range[0], x_range[1], nx)
	y = np.linspace(y_range[0], y_range[1], ny)
	x_centres = np.zeros((nx-1))
	y_centres = np.zeros((ny-1))

	cx_solps = np.zeros((nx-1, ny-1))
	cy_solps = np.zeros((nx-1, ny-1))
	leftix   = np.zeros((nx-1, ny-1))
	leftiy   = np.zeros((nx-1, ny-1))
	rightix  = np.zeros((nx-1, ny-1))
	rightiy  = np.zeros((nx-1, ny-1))
	bottomix = np.zeros((nx-1, ny-1))
	bottomiy = np.zeros((nx-1, ny-1))
	topix    = np.zeros((nx-1, ny-1))
	topiy    = np.zeros((nx-1, ny-1))
	Da_solps = np.zeros((nx-1, ny-1))

	# Convention:
	#
	# - c*.shape = (nx,ny) > cell centres
	# - values.shape = (nx,ny) in cell centres

	for ix in range(nx-1):
		for iy in range(ny-1):

			cx_solps[ix,iy] = np.mean([x[ix], x[ix+1]])
			cy_solps[ix,iy] = np.mean([y[iy], y[iy+1]])

			x_centres[ix] = cx_solps[ix,iy]
			y_centres[iy] = cy_solps[ix,iy]

			leftix[ix,iy]   = int(ix-1)
			leftiy[ix,iy]   = int(iy)
			rightix[ix,iy]  = int(ix+1)
			rightiy[ix,iy]  = int(iy)
			bottomix[ix,iy] = int(ix)
			bottomiy[ix,iy] = int(iy-1)
			topix[ix,iy]    = int(ix)
			topiy[ix,iy]    = int(iy+1)

			Da_solps[ix,iy] = Da_f2d(cx_solps[ix,iy], cy_solps[ix,iy])

	homemade_b2fgmtry = {

		"nx": nx-1,
		"ny": ny-1,

		"cx": cx_solps,
		"cy": cy_solps,
		"x_centres": x_centres,
		"y_centres": y_centres,

		"leftcut": None,

		"leftix":   leftix.astype(int), 
		"leftiy":   leftiy.astype(int),
		"rightix":  rightix.astype(int), 
		"rightiy":  rightiy.astype(int),
		"bottomix": bottomix.astype(int), 
		"bottomiy": bottomiy.astype(int),
		"topix":    topix.astype(int), 
		"topiy":    topiy.astype(int)
	}

	homemade_SOLPSsim = {"halpha_total_emission": Da_solps}

	return (homemade_b2fgmtry, homemade_SOLPSsim)