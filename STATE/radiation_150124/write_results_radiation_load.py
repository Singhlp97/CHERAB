import os
import csv
import numpy as np

from raysect.primitive.mesh.vtk import export_vtk

def write_results(cfg = None, name = None, camera = None, pipeline = None):

    # calculate power density
    frame = pipeline.frame
    power_density = frame.mean / camera.collection_areas
    power_error = frame.errors() / camera.collection_areas
    triangle_data = {'q [W/m^2]': power_density}

    basename_dir = os.path.join(cfg['output_directory_extended'], '..')

    with open(os.path.join(basename_dir, name + ".csv"), "w", newline = "") as f:
        writer = csv.writer(f)
        writer.writerow(("Triangle Index", "Power Density [W/m^2]", "Error [W/m^2]"))
        for index in range(frame.length):
            writer.writerow((index, power_density[index], power_error[index]))

    export_vtk(camera.mesh, os.path.join(basename_dir, name + '.vtk'), triangle_data = triangle_data, mode="ASCII")


