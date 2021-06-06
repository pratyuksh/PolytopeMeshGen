#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import h5py
import numpy as np
import dolfin as df
import generators.io_mesh as meshio


class Hdf5ToMfemConverter:

    def __init__(self, geom):
        
        self.ndim = geom.ndim
        self.get_bdry_attrib_id = geom.get_bdry_attrib_id
        self.get_num_bdry_attribs = geom.num_bdry_facets
        
        # number of cell nodes
        if self.ndim == 2:
            self.nn = 3
        elif self.ndim == 3:
            self.nn = 4
        
    
    # extracts boundary facets from geometry and mesh elements
    def get_mfem_bdry_attribs_2d(self, bdryNodes, bdryCells, dofMap_bdryToDomain):
        
        num_bdryCells = bdryCells.shape[0]
        
        bdry_attribs = []
        for k in range(0,num_bdryCells):
            bdrEl = bdryNodes[bdryCells[k,:]]
        
            e1 = bdryNodes[bdryCells[k,0]]
            e2 = bdryNodes[bdryCells[k,1]]
            
            id1 = dofMap_bdryToDomain[bdryCells[k,0]]
            id2 = dofMap_bdryToDomain[bdryCells[k,1]]
            
            bdry_attrib_id = self.get_bdry_attrib_id(bdrEl)
            bdry_attribs.append([ bdry_attrib_id, id1, id2 ])
        
        return np.asarray(bdry_attribs, dtype=np.int32)
    
    
    # extracts boundary facets from geometry and mesh elements
    def get_mfem_bdry_attribs_3d(self, bdryNodes, bdryCells, dofMap_bdryToDomain):
        
        num_bdryCells = bdryCells.shape[0]
        counter = np.zeros(self.get_num_bdry_attribs, dtype=np.int32)
        
        bdry_attribs = []
        for k in range(0,num_bdryCells):
            bdrEl = bdryNodes[bdryCells[k,:]]
            
            id1 = dofMap_bdryToDomain[bdryCells[k,0]]
            id2 = dofMap_bdryToDomain[bdryCells[k,1]]
            id3 = dofMap_bdryToDomain[bdryCells[k,2]]
            
            bdry_attrib_id = self.get_bdry_attrib_id(bdrEl)
            bdry_attribs.append([ bdry_attrib_id, id1, id2, id3 ])
            if (bdry_attrib_id == 0):
                print(k, bdrEl)
                exit()
            
            counter[bdry_attrib_id-1] = 1 + counter[bdry_attrib_id-1] # update counter
        
        print(counter)
        
        return np.asarray(bdry_attribs, dtype=np.int32)
    
    
    # converts from HDF5 to MFEM mesh format
    def convert_mesh_hdf5_to_mfem(self, h5_mesh_file, mfem_mesh_file):
        
        use_partition_from_file = False
        mesh = df.Mesh()
        
        inFile = df.HDF5File(mesh.mpi_comm(), h5_mesh_file, 'r')
        inFile.read(mesh, 'mesh', use_partition_from_file)
        inFile.close()
        meshio.plot_mesh(mesh)
        
        bdryMesh = df.BoundaryMesh(mesh, 'exterior')
        bdryCells = bdryMesh.cells()
        bdryNodes = bdryMesh.coordinates()
        dofMap_bdryToDomain = bdryMesh.entity_map(0).array()
        #meshio.plot_mesh(bdryMesh)
        
        # get boundary attributes
        # required for MFEM mesh format
        if self.ndim == 2:
            mfem_bdry_attribs = self.get_mfem_bdry_attribs_2d(bdryNodes, bdryCells, dofMap_bdryToDomain)
        elif self.ndim == 3:
            mfem_bdry_attribs = self.get_mfem_bdry_attribs_3d(bdryNodes, bdryCells, dofMap_bdryToDomain)    
        #print(mfem_bdry_attribs.shape)
        #print(mfem_bdry_attribs)
        
        meshio.writeMesh_mfem(mesh.coordinates(), mesh.cells(), mfem_bdry_attribs, mfem_mesh_file)


# End of file
