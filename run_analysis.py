#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import generators.io_mesh as meshio
import generators.tetrahedron_properties as tetprop


if __name__=="__main__":

    mesh_dir = "./meshes/prism_gr/weight3E-1/"
    
    h5_mesh_file = mesh_dir+"mesh_l3.h5"
    mesh = meshio.read_mesh_h5(h5_mesh_file)
    vertices = mesh.coordinates()
    cells = mesh.cells()
    
    tet_prop = tetprop.TetrahedronProperties()
    
    mesh_skewness = 0
    mesh_aspect_ratio = 1
    for j in range(0,mesh.num_cells()):
        
        tet = np.zeros((4,3))
        for l in range(0,4):
            tet[l,:] = vertices[cells[j,l]]
        
        mesh_aspect_ratio = min(mesh_aspect_ratio, tet_prop.aspect_ratio(tet))
        mesh_skewness = max(mesh_skewness, tet_prop.skewness(tet))
        
    print(mesh_aspect_ratio)
    print(mesh_skewness)
