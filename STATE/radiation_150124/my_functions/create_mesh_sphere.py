import meshzoo
from raysect.primitive.mesh.mesh import Mesh
from raysect.primitive.mesh.vtk import export_vtk


# generating triangulation of spherical surface
# source: https://fenicsproject.discourse.group/t/how-to-mesh-a-sphere-in-pymesh-or-gmsh/1662/6
# unit radius and centre in the origin assumed
vertices, triangles = meshzoo.icosa_sphere(4)
mesh = Mesh(vertices = vertices, triangles = triangles, smoothing = False, flip_normals = True)
export_vtk(mesh, "../mesh_sphere.vtk")

