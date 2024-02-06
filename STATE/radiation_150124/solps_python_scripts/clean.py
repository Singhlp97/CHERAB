import numpy as np

def clean(dct = None):

	for key in dct.keys():
		try: dct[key][np.isnan(dct[key])] = 0
		except: pass

	return dct