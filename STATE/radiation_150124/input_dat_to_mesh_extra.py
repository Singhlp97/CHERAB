import os
import csv
import numpy as np
import matplotlib.pyplot as plt

directory = "put_right_name"
file = open(os.path.join(directory, 'input.dat'), 'r')

start = False
end = False
mesh_extra = []

while end is False:

	line = file.readline() # lines[i]

	if "*** 3b. Data for additional surfaces" in line: start = True
	if "**  Reactions" in line: end = True

	is_SURFMOD = line.split('_')

	if start is True and is_SURFMOD[0] == 'SURFMOD':
		
		line = file.readline()

		while line.split()[0] == '*' and line.split()[2] == ':' and len(line.split()) == 4:

			line = file.readline()
			line = file.readline()
			
			line = file.readline()[:24]
			mesh_extra.append([float(line[:12]), float(line[12:])])

			line = file.readline()
			line = file.readline()

file.close()

mesh_extra = np.array(mesh_extra) * 1E-02
mesh_extra = np.vstack((mesh_extra, mesh_extra[0,:]))

file = open(os.path.join(directory, 'mesh.extra'), 'w+')
for data in mesh_extra:
	file.write(str(data[0]))
	file.write(" ")
	file.write(str(data[1]))
	file.write("\n")

file.close()

plt.figure()
plt.plot(mesh_extra[:,0], mesh_extra[:,1], 'ko-')
plt.axis('equal')
plt.show()
