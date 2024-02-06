import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt("collisions.txt")

collisions_max = 10
bins = np.zeros((collisions_max,1))
for n in range(collisions_max):
	bins[n] = np.sum(data[:,0] == n+1)
bins = bins / data.shape[0]

plt.figure()
plt.semilogy(np.linspace(1,collisions_max,collisions_max), bins, 'ks-')
plt.xlabel('Number of collisions')
plt.ylabel('Probability [-]')
plt.ylim(top = 1.0)

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')

# camera position [m]
start = [1.0, 0.0, 0.0]

i = 0
flag = 0

while i < data.shape[0]:

	j = 0
	points = []

	# at least one collision has occurred
	if data[i,0] > 1: 

		# store camera location and starting point
		points += [[start[0], start[1], start[2]]]
		points += [[data[i-1,1], data[i-1,2], data[i-1,3]]]

		while (i + j) < data.shape[0] and data[i+j,0] > 1:

			points += [[data[i+j,1], data[i+j,2], data[i+j,3]]]
			j += 1

	if points and np.size(points,0) == 10 and flag < 1000:
		flag += 1
		p3d = np.zeros((np.size(points,0), 3))
		for k in range(np.size(points,0)):
			p3d[k,0] = points[k][0]
			p3d[k,1] = points[k][1]
			p3d[k,2] = points[k][2]
		ax.plot(p3d[:,0], p3d[:,1], p3d[:,2])
		

	i += j + 1

ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('z [m]')
plt.show()



