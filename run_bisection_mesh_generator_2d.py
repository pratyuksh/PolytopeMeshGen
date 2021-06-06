#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import h5py
import fenics
import numpy as np

import generators.geometry2d as gm
import generators.bisection_meshes_2d as mgen
import generators.io_mesh as meshio


# Base mesh directory
base_mesh_dir = "./meshes/"

# types of bisection-tree refined meshes
br_type = ["deg0", "deg1", "deg2", "deg3"]


def gen_initial_mesh():
    # mesh_dir = base_mesh_dir + "square_twoPiecewise/"
    mesh_dir = base_mesh_dir + "square_fourPiecewise/"

    dat_mesh_file = mesh_dir + "mesh_l0.dat"
    vertices, cells = meshio.readMesh2d_dat(dat_mesh_file)

    h5_mesh_file = mesh_dir + "mesh_l0.h5"
    meshio.writeMesh_h5(vertices, cells, h5_mesh_file)


def genNestedMeshes_quasiUniform_unitSquareDomain(bool_mesh_plot, bool_mesh_write):
    mesh_dir = base_mesh_dir + "unitSquare_qu/"

    levels = np.array([1, 2, 3, 4, 5, 6, 7, 8])
    M = 2 ** levels
    print(M)

    for i in range(0, levels.size):

        mesh_file = mesh_dir + "mesh_l" + str(int(levels[i])) + ".h5"
        mesh = fenics.UnitSquareMesh(M[i], M[i])

        print(mesh.hmax(), mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


def genNestedMeshes_quasiUniform_gammaShapedDomain(bool_mesh_plot, bool_mesh_write):
    # deg = 1
    h0_init = 0.5
    angle = 1.5 * np.pi

    polygon = gm.AngularDomain(angle)
    mesh_dir = base_mesh_dir + "gammaShaped_qu/"

    levels = np.array([1, 2, 3, 4, 5, 6, 7, 8], dtype=np.int32)
    h = 1 / (2 ** levels)
    print(h)

    mesh_gen = mgen.BisectionMeshGenerator2d(polygon)
    for i in range(0, levels.size):

        mesh_file = mesh_dir + "mesh_l" + str(levels[i]) + ".h5"
        mesh = mesh_gen.bisection_uniform(h0_init, h[i])

        print(mesh.hmax(), mesh.num_vertices(), '\n', mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


def genNestedMeshes_bisecRefined_gammaShapedDomain(bool_mesh_plot, bool_mesh_write):
    deg = 1
    h0_init = 0.5
    angle = 1.5 * np.pi

    zeta = np.pi / angle
    delta = 1 - zeta

    polygon = gm.AngularDomain(angle)

    singular_corners = np.array([2, 3, 4])
    refine_weights = [0, delta, 0]
    # refine_flags = [1]*len(refine_weights)
    refine_flags = [0, 1, 0]

    mesh_dir = base_mesh_dir + "gammaShaped_br/"+br_type[deg]+"/"

    polygon.set_singular_corners(singular_corners)
    polygon.set_refine_weights(refine_weights)
    polygon.set_refine_flags(refine_flags)

    levels = np.array([1, 2, 3, 4, 5, 6, 7], dtype=np.int32)
    h = 1 / (2 ** levels)
    print(deg)
    print(h)

    mesh_gen = mgen.BisectionMeshGenerator2d(polygon)
    mesh_init = mesh_gen.uniform(h0_init)
    mesh_file = mesh_dir + "mesh_l" + str(0) + ".h5"
    if bool_mesh_plot:
        meshio.plot_mesh(mesh_init)
    if bool_mesh_write:
        meshio.write_mesh_h5(mesh_init, mesh_file)

    for i in range(0, levels.size):

        mesh_file = mesh_dir + "mesh_l" + str(levels[i]) + ".h5"
        mesh = mesh_gen.bisection(deg, h0_init, h[i])

        print(mesh.hmax(), mesh.num_vertices(), '\n', mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


def genNestedMeshes_quasiUniform_lShapedDomain(bool_mesh_plot, bool_mesh_write):
    h0_init = 0.5

    polygon = gm.LShapedDomain()
    mesh_dir = base_mesh_dir + "lShaped_qu/"

    levels = np.array([1, 2, 3, 4, 5, 6, 7, 8], dtype=np.int32)
    h = 1 / (2 ** levels)
    print(h)

    mesh_gen = mgen.BisectionMeshGenerator2d(polygon)
    for i in range(0, levels.size):

        mesh_file = mesh_dir + "mesh_l" + str(levels[i]) + ".h5"
        mesh = mesh_gen.bisection_uniform(h0_init, h[i])

        print(mesh.hmax(), mesh.num_vertices(), '\n', mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


def genNestedMeshes_bisecRefined_lShapedDomain(bool_mesh_plot, bool_mesh_write):
    deg = 3
    h0_init = 0.5
    angle = 1.5 * np.pi

    zeta = np.pi / angle
    delta = 1 - zeta

    polygon = gm.LShapedDomain()

    singular_corners = np.array([3, 4, 5])
    refine_weights = [0, delta, 0]
    # refine_flags = [1]*len(refine_weights)
    refine_flags = [0, 1, 0]

    mesh_dir = base_mesh_dir + "lShaped_br/"+br_type[deg]+"/"

    polygon.set_singular_corners(singular_corners)
    polygon.set_refine_weights(refine_weights)
    polygon.set_refine_flags(refine_flags)

    # levels = np.array([1, 2, 3, 4, 5, 6, 7], dtype=np.int32)
    levels = np.array([1, 2, 3, 4, 5, 6], dtype=np.int32)
    h = 1 / (2 ** levels)
    print(deg)
    print(h)

    mesh_gen = mgen.BisectionMeshGenerator2d(polygon)
    mesh_init = mesh_gen.uniform(h0_init)
    mesh_file = mesh_dir + "mesh_l" + str(0) + ".h5"
    if bool_mesh_plot:
        meshio.plot_mesh(mesh_init)
    if bool_mesh_write:
        meshio.write_mesh_h5(mesh_init, mesh_file)

    for i in range(0, levels.size):

        mesh_file = mesh_dir + "mesh_l" + str(levels[i]) + ".h5"
        mesh = mesh_gen.bisection(deg, h0_init, h[i])

        print(mesh.hmax(), mesh.num_vertices(), '\n', mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


def genNestedMeshes_quasiUniform_squareTwoPiecewiseDomains(bool_mesh_plot, bool_mesh_write):
    mesh_dir = base_mesh_dir + "square_twoPiecewise/qu/"

    levels = np.array([1, 2, 3, 4, 5, 6])
    h = 1 / (2 ** levels)
    print(h)

    init_mesh_file = base_mesh_dir + "square_twoPiecewise/mesh_l0.h5"
    init_mesh = meshio.read_mesh_h5(init_mesh_file)
    if bool_mesh_plot:
        meshio.plot_mesh(init_mesh)

    polygon = []
    bis_refine = mgen.BisectionRefinement(polygon)

    mesh = init_mesh
    print(mesh.hmax(), mesh.hmin(), mesh.num_vertices())
    for i in range(0, levels.size):

        mesh_file = mesh_dir + "mesh_l" + str(levels[i]) + ".h5"
        mesh = bis_refine.uniform(h[i], mesh)

        print(mesh.hmax(), mesh.hmin(), mesh.num_vertices(), '\n', mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


def genNestedMeshes_bisecRefined_squareTwoPiecewiseDomains(bool_mesh_plot, bool_mesh_write):
    mesh_dir = base_mesh_dir + "square_twoPiecewise/br/"

    deg = 2

    levels = np.array([0, 1, 2, 3, 4, 5, 6])
    h = 1 / (2 ** levels)
    print(h)

    zeta = 0.6  # singular exponent
    delta = 1 - zeta

    singular_corners = np.array([5])
    refine_weights = [delta]
    refine_flags = [1]

    init_mesh_file = base_mesh_dir + "square_twoPiecewise/mesh_l0.h5"
    init_mesh = meshio.read_mesh_h5(init_mesh_file)
    if bool_mesh_plot:
        meshio.plot_mesh(init_mesh)

    mesh_dir = base_mesh_dir + "square_twoPiecewise/br/"+br_type[deg]+"/"

    polygon = gm.SquareTwoPiecewiseDomain()
    polygon.set_singular_corners(singular_corners)
    polygon.set_refine_weights(refine_weights)
    polygon.set_refine_flags(refine_flags)

    bis_refine = mgen.BisectionRefinement(polygon)
    mesh = init_mesh
    for i in range(0, levels.size):

        mesh_file = mesh_dir + "mesh_l" + str(levels[i]) + ".h5"
        mesh = bis_refine.uniform(h[i], mesh)
        mesh = bis_refine.local(deg, h[i], mesh)

        print(mesh.hmax(), mesh.num_vertices(), '\n', mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


def genNestedMeshes_quasiUniform_squareFourPiecewiseDomains(bool_mesh_plot, bool_mesh_write):
    mesh_dir = base_mesh_dir + "square_fourPiecewise/qu/"

    levels = np.array([1, 2, 3, 4, 5, 6])
    h = 1 / (2 ** levels)
    print(h)

    init_mesh_file = base_mesh_dir + "square_fourPiecewise/mesh_l0.h5"
    init_mesh = meshio.read_mesh_h5(init_mesh_file)
    if bool_mesh_plot:
        meshio.plot_mesh(init_mesh)

    polygon = []
    bis_refine = mgen.BisectionRefinement(polygon)

    mesh = init_mesh
    print(mesh.hmax(), mesh.hmin(), mesh.num_vertices())
    for i in range(0, levels.size):

        mesh_file = mesh_dir + "mesh_l" + str(levels[i]) + ".h5"
        mesh = bis_refine.uniform(h[i], mesh)

        print(mesh.hmax(), mesh.hmin(), mesh.num_vertices(), '\n', mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


def genNestedMeshes_bisecRefined_squareFourPiecewiseDomains(bool_mesh_plot, bool_mesh_write):
    mesh_dir = base_mesh_dir + "square_fourPiecewise/br/"

    deg = 2

    levels = np.array([1, 2, 3, 4, 5, 6])
    h = 1 / (2 ** levels)
    print(h)

    zeta = 0.6  # singular exponent
    delta = 1 - zeta

    singular_corners = np.array([5])
    refine_weights = [delta]
    refine_flags = [1]

    init_mesh_file = base_mesh_dir + "square_fourPiecewise/mesh_l0.h5"
    init_mesh = meshio.read_mesh_h5(init_mesh_file)
    if bool_mesh_plot:
        meshio.plot_mesh(init_mesh)

    mesh_dir = base_mesh_dir + "square_fourPiecewise/br/"+br_type[deg]+"/"

    polygon = gm.SquareTwoPiecewiseDomain()
    polygon.set_singular_corners(singular_corners)
    polygon.set_refine_weights(refine_weights)
    polygon.set_refine_flags(refine_flags)

    bis_refine = mgen.BisectionRefinement(polygon)
    mesh = init_mesh
    for i in range(0, levels.size):

        mesh_file = mesh_dir + "mesh_l" + str(levels[i]) + ".h5"
        mesh = bis_refine.uniform(h[i], mesh)
        mesh = bis_refine.local(deg, h[i], mesh)

        print(mesh.hmax(), mesh.num_vertices(), '\n', mesh_file)
        if bool_mesh_plot:
            meshio.plot_mesh(mesh)
        if bool_mesh_write:
            meshio.write_mesh_h5(mesh, mesh_file)


if __name__ == '__main__':
    bool_mesh_plot_ = False
    bool_mesh_write_ = True

    # genNestedMeshes_quasiUniform_unitSquareDomain(bool_mesh_plot_, bool_mesh_write_)
    
    # genNestedMeshes_quasiUniform_gammaShapedDomain(bool_mesh_plot_, bool_mesh_write_)
    # genNestedMeshes_bisecRefined_gammaShapedDomain(bool_mesh_plot_, bool_mesh_write_)
    
    # genNestedMeshes_quasiUniform_lShapedDomain(bool_mesh_plot_, bool_mesh_write_)
    genNestedMeshes_bisecRefined_lShapedDomain(bool_mesh_plot_, bool_mesh_write_)

    # gen_initial_mesh()
    # genNestedMeshes_quasiUniform_squareTwoPiecewiseDomains(bool_mesh_plot_, bool_mesh_write_)
    # genNestedMeshes_bisecRefined_squareTwoPiecewiseDomains(bool_mesh_plot_, bool_mesh_write_)
    
    # gen_initial_mesh()
    # genNestedMeshes_quasiUniform_squareFourPiecewiseDomains(bool_mesh_plot_, bool_mesh_write_)
    # genNestedMeshes_bisecRefined_squareFourPiecewiseDomains(bool_mesh_plot_, bool_mesh_write_)

# END OF FILE
