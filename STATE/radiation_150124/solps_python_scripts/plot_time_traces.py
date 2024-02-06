import os
import time
import numpy as np
import matplotlib.pyplot as plt

def plot_time_traces(where = ".", what = None, routine = "2dt"):

    # Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	if not isinstance(where, list): where = [where]

	variables = ["nesepm", "tesepm", "tisepm",
				 "nemxap", "temxap", "timxap",
				 "nemxip", "temxip", "timxip",
				 "pwmxap", "pwmxip", 
				 "tmne",   "tmte",   "tmti"]

	if what is not None and isinstance(what, str):
		what = [what]

	variables = what or variables

	for here in where:

		original = os.getcwd()
		try: os.chdir(here)
		except: pass

		for variable in variables:

			command = routine + " " + variable + " &"
			os.system(command)
			time.sleep(1.0)

			data = np.loadtxt("./gnuplot.data")

			try:

				threshold = 1E+36
				index = np.where(data > threshold)[0][0]
				data = np.concatenate((data[:index,:], data[index+1:,:]), axis = 0)
				plt.figure()
				if data.shape[1] == 2:
					plt.plot(data[:,0], data[:,1], 'ko-')
				elif data.shape[1] == 3:
					plt.plot(data[:,0], data[:,1], 'ko-', label = "bottom")
					plt.plot(data[:,0], data[:,2], 'mo-', label = "top")
					plt.legend(loc = 'best')
				plt.xlabel('Time [s]')
				plt.title(variable)

			except:
				pass

		try: os.chdir(original)
		except: pass

	plt.show()

	return