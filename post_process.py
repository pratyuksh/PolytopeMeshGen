#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import io_mesh as meshio
# import numpy as np

baseMeshDir = "/userdata/prbansal/Documents/projectsOngoing/meshes/"
mfemMeshDir = "/userdata/prbansal/Documents/projectsOngoing/meshes_MFEM/"


# converts the mesh from .h5 to paraview format
def convert_meshes_h5_to_pvd():
    # meshDir = baseMeshDir + "prism_gr/weight3E-1/"
    meshDir = baseMeshDir + "ficheraCube_gr/weight3E-1/"

    L = 4
    for k in range(0, L):
        h5MeshFile = meshDir + "mesh_l%d.h5" % k
        pvd_mesh_file = meshDir + "mesh_l%d.pvd" % k
        meshio.convert_mesh_h5_to_pvd(h5MeshFile, pvd_mesh_file)


# converts the mesh from .h5 to .vtk format
def convert_meshes_h5_to_vtk():
    # inMeshDir = baseMeshDir + "prism_gr/weight5E-1/"
    # outMeshDir = mfemMeshDir + "prism_gr/weight5E-1/"

    inMeshDir = baseMeshDir + "fichera_cube_gr/weight3E-1/"
    outMeshDir = mfemMeshDir + "fichera_cube_gr/weight3E-1/"

    l0 = 6
    L = 7
    for k in range(l0, L):
        h5MeshFile = inMeshDir + "mesh_l%d.h5" % k
        vtkMeshFile = outMeshDir + "mesh_l%d.vtk" % k
        meshio.convert_mesh_h5_to_vtk(h5MeshFile, vtkMeshFile)


    meshio.plot_mesh_file(xmlFile)


if __name__ == "__main__":
    # convert_meshes_h5_to_pvd()
    convert_meshes_h5_to_vtk()

# End of file
