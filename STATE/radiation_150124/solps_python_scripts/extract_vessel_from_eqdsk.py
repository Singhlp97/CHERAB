import os
import numpy as np
import matplotlib.pyplot as plt
from cherab.tools.equilibrium import import_eqdsk

def extract_vessel_from_eqdsk(where = ".", filename = None):

    if '.eqdsk' not in filename and 'geqdsk' not in filename:
        filename = '.'.join([filename, '.eqdsk'])

    eqdsk = import_eqdsk(os.path.join(where, filename))
    wall = eqdsk.limiter_polygon

    # add 0th point to the end and convert to [mm]
    wall = np.concatenate((wall, np.array([wall[0,:]])), axis = 0) * 1E+03

    f = open('.'.join([filename.split('.')[0], 'ogr']), 'w')
    for i in range(wall.shape[0]):
        f.write('{} {}\n'.format(wall[i,0], wall[i,1]))
    f.close()

    fig, ax = plt.subplots()
    ax.plot(wall[:,0], wall[:,1], 'ko-', linewidth = 2)
    ax.set_xlabel('$R \\; [m]$')
    ax.set_ylabel('$Z \\; [m]$')
    ax.set_aspect('equal')
    ax.set_title('Wall profile')

    plt.show()

    return

