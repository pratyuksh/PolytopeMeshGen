#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import h5py
import numpy as np
import time

import generators.graded_meshes_2d as mgen
import generators.geometry2d as gm
import generators.io_mesh as meshio


# Base mesh directory
base_mesh_dir = "./meshes/"

# types of graded meshes
gr_type = ["weight3E-1", "weight1E-1", "weight4E-2"]


def gen_initial_mesh():
    # mesh_dir = base_mesh_dir + "gammaShaped_gr/"
    mesh_dir = base_mesh_dir + "lShaped_gr/"

    dat_mesh_file = mesh_dir + "mesh_l0.dat"
    vertices, cells = meshio.readMesh2d_dat(dat_mesh_file)

    h5_mesh_file = mesh_dir + "mesh_l0.h5"
    meshio.writeMesh_h5(vertices, cells, h5_mesh_file)


def genNestedMeshes_graded_gammaShapedDomain(bool_mesh_plot, bool_mesh_write):
    deg = 3
    angle = 3 * np.pi / 2

    singular_corners = np.array([2, 3, 4])

    mesh_dir = base_mesh_dir + "gammaShaped_gr/"+gr_type[deg-1]+"/"

    refine_weights = None
    refine_weights_uniform = [0.5, 0.5, 0.5]
    if deg == 1:
        refine_weights = [0.5, 0.3, 0.5]
    elif deg == 2:
        refine_weights = [0.5, 0.1, 0.5]
    elif deg == 3:
        refine_weights = [0.5, 0.04, 0.5]

    polygon = gm.AngularDomain(angle)
    polygon.set_singular_corners(singular_corners)
    # polygon.set_refine_weights(refine_weights)

    l0 = 0
    init_mesh_file = base_mesh_dir + "gammaShaped_gr/mesh_l%d.h5" % l0
    init_mesh = meshio.read_mesh_h5(init_mesh_file)
    if bool_mesh_plot:
        meshio.plot_mesh(init_mesh)

    mesh_gen = mgen.GradedMeshGenerator2d(polygon)

    mesh = init_mesh
    num_refinements = 8
    for k in range(0, num_refinements):

        if k + l0 + 1 <= 1:
            polygon.set_refine_weights(refine_weights_uniform)
        else:
            polygon.set_refine_weights(refine_weights)

        print('Running mesh refinement ...')

        print('\tTotal number of vertices in Initial Mesh: %d' % mesh.num_vertices())
        print('\tTotal number of cells in Initial Mesh: %d' % mesh.num_cells())

        vertices, cells = mesh_gen.refine(mesh)

        print('\n\tTotal number of vertices in Refined Mesh: %d' % vertices.shape[0])
        print('\tTotal number of cells in Refined Mesh: %d\n' % cells.shape[0])

        mesh = mesh_gen.create(vertices, cells)
        mesh_file = mesh_dir + "mesh_l%d.h5" % (k + l0 + 1)
        print(mesh.hmax(), mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


def genNestedMeshes_graded_lShapedDomain(bool_mesh_plot, bool_mesh_write):
    deg = 3
    angle = 3 * np.pi / 2

    singular_corners = np.array([3, 4, 5])

    mesh_dir = base_mesh_dir + "lShaped_gr/"+gr_type[deg-1]+"/"

    refine_weights = None
    refine_weights_uniform = [0.5, 0.5, 0.5]
    if deg == 1:
        refine_weights = [0.5, 0.3, 0.5]
    elif deg == 2:
        refine_weights = [0.5, 0.1, 0.5]
    elif deg == 3:
        refine_weights = [0.5, 0.04, 0.5]

    polygon = gm.LShapedDomain()
    polygon.set_singular_corners(singular_corners)
    # polygon.set_refine_weights(refine_weights)

    l0 = 0
    init_mesh_file = base_mesh_dir + "lShaped_gr/mesh_l%d.h5" % l0
    init_mesh = meshio.read_mesh_h5(init_mesh_file)
    if bool_mesh_plot:
        meshio.plot_mesh(init_mesh)

    mesh_gen = mgen.GradedMeshGenerator2d(polygon)

    mesh = init_mesh
    # num_refinements = 8
    num_refinements = 7
    for k in range(0, num_refinements):

        if k + l0 + 1 <= 1:
            polygon.set_refine_weights(refine_weights_uniform)
        else:
            polygon.set_refine_weights(refine_weights)

        print('Running mesh refinement ...')

        print('\tTotal number of vertices in Initial Mesh: %d' % mesh.num_vertices())
        print('\tTotal number of cells in Initial Mesh: %d' % mesh.num_cells())

        vertices, cells = mesh_gen.refine(mesh)

        print('\n\tTotal number of vertices in Refined Mesh: %d' % vertices.shape[0])
        print('\tTotal number of cells in Refined Mesh: %d\n' % cells.shape[0])

        mesh = mesh_gen.create(vertices, cells)
        mesh_file = mesh_dir + "mesh_l%d.h5" % (k + l0 + 1)
        print(mesh.hmax(), mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


if __name__ == "__main__":
    gen_initial_mesh()

    bool_mesh_plot_ = False
    bool_mesh_write_ = True

    start_time = time.clock()
    # genNestedMeshes_graded_gammaShapedDomain(bool_mesh_plot_, bool_mesh_write_)
    genNestedMeshes_graded_lShapedDomain(bool_mesh_plot_, bool_mesh_write_)
    print("Elapsed time: ", time.clock() - start_time)

# End of file
