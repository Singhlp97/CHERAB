import os
import csv
import numpy as np

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'figure.autolayout': True})

from raysect.optical.observer import MeshCamera
from raysect.primitive.mesh.vtk import import_vtk, export_vtk


from raysect.optical.observer import MeshCamera

import pylab as plb
from scipy.optimize import curve_fit
from scipy import asarray as ar, exp 


#####################################################################################################################
########################## RAW DATA REPRESENTATION##############################################################################################
#####################################################################################################################

# locations = ["./output/H_alpha/hal.step=0.abs=0.astra=0/", "./output/H_alpha/hal.step=1.abs=0.astra=0/"]
locations = ["output/9229_equ@16ms_HW/halpha_total_radiation/step=True.abs=False.astra=False.eirene=True.lens=False.DiscreteToroidalMapper.1/",
		     "output/9229_equ@16ms_HW/halpha_total_radiation/step=True.abs=False.astra=False.eirene=True.lens=True.DiscreteToroidalMapper.1/"]
# labels = ["Uniform $\Delta s$", "Smart $\Delta s$"]
# labels = ["Transparent", "Absorbing"]
labels = ["Standard", "With lens"]

colors = ["k", "b", "m", "c"]
ic = 0

for index in range(len(locations)):

	location = locations[index]
	label = labels[index]

	basename = "camera_shift"

	dq = np.loadtxt(location + basename + '_errors.csv', delimiter = ',', skiprows = 0)
	q = np.loadtxt(location + basename + '.csv', delimiter = ',', skiprows = 0)

	nx = dq.shape[0]
	ny = dq.shape[1]
	data = np.zeros((nx * ny, 3))

	for i in range(nx):
		for j in range(ny):
			k = j + i * nx
			data[k, 1] = q[i, j] 
			data[k, 2] = dq[i, j] 

	num_data = np.size(data, 0)

	result = data[:,1]
	result_errors = data[:,2]
	result = result_errors / result
	result[np.isnan(result)] = 0

	weighted_avg = data[:,2].sum() / data[:,1].sum()

	print()
	print('##########')
	print(location)
	print('Average error: {:.4G}%'.format(weighted_avg * 1E+02))
	print('##########')
	print()

	# look at the result of the equatorial limiter with the 1x1 source: there are no zeros
	# on the front face => the number of zeros is the one of the only shadowed edge

	#####################################################################################################################

	min_result = min(result[result > 0])
	max_result = max(result)

	num_bins = int(300)

	bins = np.zeros((num_bins))

	# (num_bins) NOT (+1) because the last one is added separately
	#bins_power_range = np.linspace(min_result, max_result, num_bins)
	bins_power_range = np.logspace(np.log10(min_result), np.log10(max_result), num_bins)
	# more easily manipulated
	bins_power_range = list(bins_power_range)
	delta_bins_range = bins_power_range[1] - bins_power_range[0]
	# you need a further one for the maximum power
	bins_power_range += [max_result + delta_bins_range]

	# watch out the power range starts at the left edge of the bin

	#####################################################################################################################

	for i in range(num_data): 

		pw = result[i]

		for j in range(num_bins):

			# watch out where = and where not
			# faster if < is tested at first
			if j < (num_bins - 1) and pw < bins_power_range[j + 1] and pw >= bins_power_range[j] and pw > 0:
				bins[j] += 1
				# you can exit the loop and check the following power
				break
		
			#elif j == (num_bins - 1) and pw > 0:
				# it is for sure in the last one, i.e. the maximum power
				#bins[-1] += 1

	#####################################################################################################################

	# normalization
	bins = np.array(bins)
	bins_power_range = np.array(bins_power_range)
	bins[0] = 0 # removing the meaningless zeros

	integral = sum(bins * (bins_power_range[1:] - bins_power_range[:-1]))

	bins /= integral

	CONVERSION = 1

	bins_power_range /= CONVERSION # from W/m^2 to MW/m^2

	fig1 = plt.figure(1)
	plt.semilogx(bins_power_range[0:-1], bins * CONVERSION, '.', color = colors[ic], label = label)
	
	fig2 = plt.figure(2)
	pw = data[:,1]
	ic += 1
	plt.semilogy(pw, result, '.', color = colors[ic], label = label)
	ic += 1
	
plt.figure(1)
plt.xlabel('$\epsilon\;[-]$')
plt.ylabel('$f\:(\epsilon)\;[-]$')
plt.title('PDF: relative error')
plt.legend(loc = 'best', markerscale = 5)
plt.grid(True)

fig1.savefig(location + basename + '_errors_PDF.png')

plt.figure(2)
plt.xlabel('$\Gamma$ [photons/s]')
plt.ylabel('$\epsilon\;[-]$')
plt.title('Value-error scatter plot')
plt.legend(loc = 'best', markerscale = 5)
plt.grid(True)

fig2.savefig(location + basename + '_errors_corr.png')


plt.show()
