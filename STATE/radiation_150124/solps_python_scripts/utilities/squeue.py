import os

def squeue(user = "matteo.moscheni"):

	# Author: Matteo Moscheni
    # E-mail: matteo.moscheni@tokamakenergy.co.uk
    # February 2022

	os.system("squeue -u " + user + " > squeue.tmp")

	f = open('squeue.tmp', 'r')

	JOBIDs = []

	line = f.readline()

	for line in f: JOBIDs.append(line.split()[0])
	f.close()
	os.system("rm squeue.tmp")

	print()

	os.system('squeue -u matteo.moscheni')

	for JOBID in JOBIDs:
		print()
		print("   JOBID: {}".format(JOBID))
		command = 'scontrol show JOBID ' + JOBID + ' | grep -e Work'
		os.system(command)

	print()
	return