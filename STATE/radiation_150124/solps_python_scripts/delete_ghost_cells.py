import numpy as np

def delete_ghost_cells(dct = None):

	for key in dct.keys():
		try:
			if   len(dct[key].shape) == 2:
				dct[key] = dct[key][1:-1,1:-1]
			elif len(dct[key].shape) == 3:
				dct[key] = dct[key][1:-1,1:-1,:]
			elif len(dct[key].shape) == 5:
				dct[key] = dct[key][1:-1,1:-1,:,:]
			else:
				pass
		except: pass

	return dct