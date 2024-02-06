import os
import numpy as np
import math

# from my_functions.probability_density_function import probability_density_function

#######################################################################################################################


#location = '/home/moscheni/DISRUPTIONS/EUROFusion/from_Matthew/from_zip/results_example/'
#location = './TEST_with_reflection/data_1_10cm_10000_off1e-6_vpcen/'
location = './radiation_load/output/Fabio/total_radiation/2.step=True.abs=False.scat=False.astra=False.b2=True.eirene=False.extra=False.AxisymmetricMapper/'

run_files = os.listdir(location)

rel_errors = []

for file in run_files:

	filename, ext = os.path.splitext(file)

	if ext == '.csv':

		print("importing", file)

		powers = np.loadtxt(open(os.path.join(location + '/', filename + '.csv'), "r"), delimiter=",", skiprows=1)

		num_detectors = np.size(powers, 0)
		
		#triangles = powers[:, 0]
		powers_m2 = powers[:, 1]
		power_m2_errors = powers[:, 2]

		for i in range(num_detectors):
			if powers_m2[i] > 0:
				rel_errors += [power_m2_errors[i] / powers_m2[i]]

AVG_error = np.mean(rel_errors)

print('\n\nAverage error: {:.8G} %'.format(AVG_error * 100))
print()