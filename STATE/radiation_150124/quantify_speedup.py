import os
import numpy as np
import matplotlib.pyplot as plt

location = "./output/"
basename = "H_alpha"
outdirs = [basename, basename + '_x02', basename + '_x05', basename + '_x10', basename + '_x20']

dic = {}

for outdir in outdirs:

	dic[outdir] = {}

	runlogs = os.listdir(os.path.join(location, outdir))

	for runlog in runlogs:

		details = runlog.split(".")

		if details[0] == "run":

			details = details[2]

			dic[outdir][details] = {}

			f = open(os.path.join(location, outdir, runlog), "r")

			# go to last line to get details

			for line in f:
				if line != "\n" and line.split()[0] == "Render":
					line_prev = line
			line = line_prev
			line = line.split()

			# normalise everything to 1 min run time

			elaps_time = float(line[2].split("s")[0])
			percentage = float(line[3].split("(")[1].split("%")[0]) * 1E-02

			elaps_time_norm = 60.0
			normalisation = elaps_time_norm / elaps_time
			percentage_norm = percentage * normalisation

			dic[outdir][details]['elaps_time'] = elaps_time
			dic[outdir][details]['percentage'] = percentage

			dic[outdir][details]['elaps_time_norm'] = elaps_time_norm
			dic[outdir][details]['percentage_norm'] = percentage_norm

# re-organise data for plotting

num_cases = len(outdirs)

time_0 = [] # step=1E-04_1E-02
time_1 = [] # step=1E-04
time_2 = [] # step=1E-05_1E-03
time_3 = [] # step=1E-05

for outdir in dic:
	for details in dic[outdir]:

		time = dic[outdir][details]['elaps_time_norm']
		perc = dic[outdir][details]['percentage_norm']

		# total time that would have taken to complete the simulation

		if details == "step=1E-04_1E-02": time_0 += [time / perc] 
		if details == "step=1E-04":       time_1 += [time / perc] 
		if details == "step=1E-05_1E-03": time_2 += [time / perc] 
		if details == "step=1E-05":       time_3 += [time / perc] 

increase_factor = [1, 2, 5, 10, 20]

time_0 = np.array(time_0)
time_1 = np.array(time_1)
time_2 = np.array(time_2)
time_3 = np.array(time_3)

fig = plt.figure()

plt.semilogy(increase_factor, time_0, 'ko-', label = "step=1E-04_1E-02")
plt.semilogy(increase_factor, time_1, 'bo-', label = "step=1E-04")
plt.semilogy(increase_factor, time_2, 'co-', label = "step=1E-05_1E-03")
plt.semilogy(increase_factor, time_3, 'mo-', label = "step=1E-05")

plt.xlabel("size increase factor $[-]$")
plt.ylabel("estimated computational time $[s]$")
plt.legend()

fig = plt.figure()

plt.plot(increase_factor, time_1 / time_0, 'ko-', label = "step=1E-04_1E-02")
plt.plot(increase_factor, time_3 / time_2, 'co-', label = "step=1E-05_1E-03")

plt.xlabel("size increase factor $[-]$")
plt.ylabel("absolute speed-up factor $[-]$")
plt.legend()

fig = plt.figure()

r = time_1 / time_0
plt.plot(increase_factor[1:], r[1:]/r[0], 'ko-', label = "step=1E-04_1E-02")
r = time_3 / time_2
plt.plot(increase_factor[1:], r[1:]/r[0], 'co-', label = "step=1E-05_1E-03")

plt.xlabel("size increase factor $[-]$")
plt.ylabel("relative speed-up factor $[-]$")
plt.legend()

plt.show()