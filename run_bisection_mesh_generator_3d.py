#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import h5py
import numpy as np
import dolfin as df
import mshr
import time

import generators.bisection_meshes_3d as mgen
import generators.geometry3d as gm
import generators.io_mesh as meshio


# Base mesh directory
base_mesh_dir = "./meshes/"


def gen_initial_mesh_from_dat_file():
    mesh_dir = base_mesh_dir + "prism_quasiUniform/"

    dat_mesh_file = mesh_dir + "mesh_l0.dat"
    vertices, cells = meshio.readMesh3d_dat(dat_mesh_file)

    h5_mesh_file = mesh_dir + "mesh_l0.h5"
    meshio.writeMesh_h5(vertices, cells, h5_mesh_file)


def gen_initial_mesh_from_off_file():
    mesh_dir = base_mesh_dir + "prism_quasiUniform/"

    off_mesh_file = mesh_dir + "mesh_l0.off"
    domain = mshr.Surface3D(off_mesh_file)
    mesh = mshr.generate_mesh(domain, 5)

    polyhedron = gm.PrismDomain()
    mesh_gen = mgen.BisectionMeshGenerator3d(polyhedron)
    mesh = mesh_gen.bisection_uniform(mesh, 0.5)
    print('Max meshwidth: ', mesh.hmax())
    print('Min meshwidth: ', mesh.hmin())

    # h5_mesh_file = mesh_dir+"mesh_l0.h5"
    # meshio.write_mesh_h5(mesh, h5_mesh_file)


def genNestedMeshes_quasiUniform_unitCube(bool_mesh_plot, bool_mesh_write):
    mesh_dir = base_mesh_dir + "cube_quasiUniform/"

    levels = np.array([1, 2, 3, 4, 5, 6])
    M = 2 ** levels
    print(M)

    for i in range(0, levels.size):

        mesh_file = mesh_dir + "mesh_l" + str(int(levels[i])) + ".h5"
        mesh = df.UnitCubeMesh(M[i], M[i], M[i])

        print(mesh.hmax(), mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


def genNestedMeshes_quasiUniform_prism(bool_mesh_plot, bool_mesh_write):
    # deg = 1
    # h0_init = 0.5
    # angle = 1.5 * np.pi

    polyhedron = gm.PrismDomain()
    mesh_dir = base_mesh_dir + "prism_quasiUniform/"

    l0 = 0
    init_mesh_file = base_mesh_dir + "prism_quasiUniform/mesh_l%d.h5" % l0
    init_mesh = meshio.read_mesh_h5(init_mesh_file)
    if bool_mesh_plot:
        meshio.plot_mesh(init_mesh)

    levels = np.array([1, 2, 3, 4], dtype=np.int32)
    # levels = np.array([1,2,3,4,5,6,7,8], dtype=np.int32)
    h = 1 / (2 ** levels)
    print(h)

    mesh_gen = mgen.BisectionMeshGenerator3d(polyhedron)

    mesh = init_mesh
    for i in range(0, levels.size):

        mesh_file = mesh_dir + "mesh_l" + str(levels[i]) + ".h5"
        mesh = mesh_gen.bisection_uniform(mesh, h[i])

        print(mesh.hmax(), mesh.num_vertices(), '\n', mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


def genNestedMeshes_bisecRefine_prism(bool_mesh_plot, bool_mesh_write):
    deg = 2
    # h0_init = 0.5
    angle = 1.5 * np.pi

    zeta = np.pi / angle
    delta = 1 - zeta

    polyhedron = gm.PrismDomain()

    singular_edges = np.array([13])
    refine_weights = [delta]
    refine_flags = [1] * len(refine_weights)

    mesh_dir = None
    if deg == 1:
        mesh_dir = base_mesh_dir + "prism_bisecRefine/linear/"
    elif deg == 2:
        mesh_dir = base_mesh_dir + "prism_bisecRefine/quadratic/"
    elif deg == 3:
        mesh_dir = base_mesh_dir + "prism_bisecRefine/cubic/"

    polyhedron.set_singular_edges(singular_edges)
    polyhedron.set_refine_weights([], refine_weights)
    polyhedron.set_refine_flags(refine_flags)

    l0 = 0
    init_mesh_file = base_mesh_dir + "prism_quasiUniform/mesh_l%d.h5" % l0
    init_mesh = meshio.read_mesh_h5(init_mesh_file)
    if bool_mesh_plot:
        meshio.plot_mesh(init_mesh)

    levels = np.array([1, 2, 3], dtype=np.int32)
    # levels = np.array([1,2,3,4,5,6,7,8], dtype=np.int32)
    h = 1 / (2 ** levels)
    print(h)

    mesh_gen = mgen.BisectionMeshGenerator3d(polyhedron)

    mesh = init_mesh
    for i in range(0, levels.size):

        mesh_file = mesh_dir + "mesh_l" + str(levels[i]) + ".h5"
        mesh = mesh_gen.bisection(deg, mesh, h[i])

        print(mesh.hmax(), mesh.num_vertices(), '\n', mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


if __name__ == "__main__":
    # gen_initial_mesh_from_dat_file()
    # gen_initial_mesh_from_off_file()

    bool_mesh_plot_ = True
    bool_mesh_write_ = False

    start_time = time.clock()
    # genNestedMeshes_quasiUniform_unitCube(bool_mesh_plot_, bool_mesh_write_)
    # genNestedMeshes_quasiUniform_prism(bool_mesh_plot_, bool_mesh_write_)
    genNestedMeshes_bisecRefine_prism(bool_mesh_plot_, bool_mesh_write_)
    print("Elapsed time: ", time.clock() - start_time)

    # h5_mesh_file = base_mesh_dir+"prism_quasiUniform/mesh_l0.h5"
    # pvd_mesh_file = base_mesh_dir+"prism_quasiUniform/mesh_l0.pvd"
    # h5_mesh_file = base_mesh_dir+"prism_bisecRefine/linear/mesh_l2.h5"
    # pvd_mesh_file = base_mesh_dir+"prism_bisecRefine/linear/mesh_l2.pvd"
    # meshio.convert_mesh_h5_to_pvd(h5_mesh_file, pvd_mesh_file)

# END OF FILE
