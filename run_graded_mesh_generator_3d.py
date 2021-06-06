#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import h5py
import numpy as np
import time

import generators.graded_meshes_3d as mgen
import generators.geometry3d as gm
import generators.io_mesh as meshio

# Base mesh directory
base_mesh_dir = "./meshes/"


def gen_initial_mesh():
    mesh_dir = base_mesh_dir+"cube_quasiUniform/"
    # mesh_dir = base_mesh_dir+"prism_graded/"
    # mesh_dir = base_mesh_dir+"fichera_cube_graded/"

    dat_mesh_file = mesh_dir + "mesh_l0.dat"
    vertices, cells = meshio.readMesh3d_dat(dat_mesh_file)

    h5_mesh_file = mesh_dir + "mesh_l0.h5"
    meshio.writeMesh_h5(vertices, cells, h5_mesh_file)


def gen_refined_mesh(mesh_dir, mesh_gen, k0, num_refinements, bool_quasiuniform, bool_mesh_plot, bool_mesh_write):
    
    init_mesh_file = mesh_dir + 'mesh_l%d.h5' % k0
    mesh = meshio.read_mesh_h5(init_mesh_file)
    if bool_mesh_plot:
        meshio.plot_mesh(mesh)
    
    for k in range(0, num_refinements):
        
        print('Running mesh refinement ...')
        
        print('\tTotal number of vertices in Initial Mesh: %d' % mesh.num_vertices())
        print('\tTotal number of cells in Initial Mesh: %d' % mesh.num_cells())
        
        vertices, cells = mesh_gen.refine(mesh, bool_quasiuniform)
        
        print('\n\tTotal number of vertices in Refined Mesh: %d' % vertices.shape[0])
        print('\tTotal number of cells in Refined Mesh: %d\n' % cells.shape[0])
        
        mesh = mesh_gen.create(vertices, cells)
        mesh_file = mesh_dir+"mesh_l%d.h5" % (k+k0+1)
        print(mesh.hmax(), mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)
            # meshio.writeMesh_h5(vertices, cells, meshFile)


def genNestedMeshes_cube_graded(bool_quasiuniform, bool_mesh_plot, bool_mesh_write):
    
    if bool_quasiuniform:
        mesh_dir = base_mesh_dir+"cube_gr/weight5E-1/"
    else:
        mesh_dir = base_mesh_dir+"cube_gr/weight3E-1/"
    print("\nMesh directory: ", mesh_dir, "\n")

    polyhedron = gm.CubeDomain()
    sin_corners_refine_weights = [0.5, 0.5]
    sin_edges_refine_weights = [0.3]
    polyhedron.set_refine_weights(sin_corners_refine_weights, sin_edges_refine_weights)
    
    l0 = 0
    num_refinements = 3
    mesh_gen = mgen.GradedMeshGenerator3d(polyhedron)
    gen_refined_mesh(mesh_dir, mesh_gen, l0, num_refinements, bool_quasiuniform, bool_mesh_plot, bool_mesh_write)
    

def genNestedMeshes_prism_graded(bool_quasiuniform, bool_mesh_plot, bool_mesh_write):
    
    if bool_quasiuniform:
        mesh_dir = base_mesh_dir+"prism_gr/weight5E-1/"
    else:
        mesh_dir = base_mesh_dir+"prism_gr/weight3E-1/"
    print("\nMesh directory: ", mesh_dir, "\n")
    
    polyhedron = gm.PrismDomain()
    sin_corners_refine_weights = [0.5, 0.5]
    sin_edges_refine_weights = [0.3]
    polyhedron.set_refine_weights(sin_corners_refine_weights, sin_edges_refine_weights)
    
    l0 = 0
    num_refinements = 6
    mesh_gen = mgen.GradedMeshGenerator3d(polyhedron)
    gen_refined_mesh(mesh_dir, mesh_gen, l0, num_refinements, bool_quasiuniform, bool_mesh_plot, bool_mesh_write)


def genNestedMeshes_fichera_cube_graded(bool_quasiuniform, bool_mesh_plot, bool_mesh_write):
    
    if bool_quasiuniform:
        mesh_dir = base_mesh_dir+"ficheraCube_gr/weight5E-1/"
    else:
        mesh_dir = base_mesh_dir+"ficheraCube_gr/weight3E-1/"
    print("\nMesh directory: ", mesh_dir, "\n")
    
    polyhedron = gm.FicheraCubeDomain()
    corners_refine_weights = [0.3]
    sing_edges_refine_weights = [0.3, 0.3, 0.3]
    polyhedron.set_refine_weights(corners_refine_weights, sing_edges_refine_weights)
    
    l0 = 0
    num_refinements = 6
    mesh_gen = mgen.GradedMeshGenerator3d(polyhedron)
    gen_refined_mesh(mesh_dir, mesh_gen, l0, num_refinements, bool_quasiuniform, bool_mesh_plot, bool_mesh_write)


if __name__ == "__main__":
    
    # gen_initial_mesh()
    
    bool_quasiuniform_ = False
    bool_mesh_plot_ = False
    bool_mesh_write_ = True
    
    startTime = time.process_time()
    # genNestedMeshes_cube_graded(bool_quasiuniform_, bool_mesh_plot_, bool_mesh_write_)
    # genNestedMeshes_prism_graded(bool_quasiuniform_, bool_mesh_plot_, bool_mesh_write_)
    genNestedMeshes_fichera_cube_graded(bool_quasiuniform_, bool_mesh_plot_, bool_mesh_write_)
    print("Elapsed time: ", time.process_time()-startTime)
    
    
# End of file
