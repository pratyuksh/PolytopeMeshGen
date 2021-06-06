#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import converters.geometry as gm
import converters.converter_hdf5_to_mfem as mconv
import generators.io_mesh as meshio


def call_mesh_hdf5_to_mfem_3d():
    
    xFem_type = ["linear", "quadratic", "cubic"]
    graded_type = ["3e-1", "2e-1", "1e-1"]
    
    l = 7
    idx = 2
    
    geom_file = "geometry_files/prism.txt"
    h5_mesh_file = "meshes/prism_graded/weight_5e-1/mesh_l%d.h5"%l
    mfem_mesh_file = "meshes_MFEM/prism_graded/weight_5e-1/mesh_l%d.mesh"%l
    
    #geom_file = "geometry_files/prism.txt"
    #h5_mesh_file = "meshes/prism_graded/weight_"+graded_type[idx]+"/mesh_l%d.h5"%l
    #mfem_mesh_file = "meshes_MFEM/prism_graded/weight_"+graded_type[idx]+"/mesh_l%d.mesh"%l
    
    print(geom_file)
    print(h5_mesh_file)
    print(mfem_mesh_file)
    
    geom = gm.Geometry3d(geom_file)
    preproc_mesh3d = Preprocess_mesh(geom)
    preproc_mesh3d.convert_mesh_hdf5_to_mfem(h5_mesh_file, mfem_mesh_file)
    
   
if __name__=='__main__':
    
    call_mesh_hdf5_to_mfem_3d()


# End of file
