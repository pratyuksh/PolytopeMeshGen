#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib import rcParams
font = {'family': 'Dejavu Sans',
        'weight': 'normal',
        'size': 30}
rc('font', **font)
rcParams['lines.linewidth'] = 4
rcParams['lines.markersize'] = 18
rcParams['markers.fillstyle'] = 'none'
rcParams.update({'figure.autolayout': True})

import generators.io_mesh as meshio


# Mesh directory
meshes_dir = "./meshes/"

# Output directory
output_dir = "./figures/"


def save_plot_mesh2d(mesh, filename):
    x, y = mesh.coordinates()[:, 0], mesh.coordinates()[:, 1]
    plt.triplot(x, y, triangles=mesh.cells())

    fig = plt.gcf()
    fig.set_size_inches(15, 15)
    fig.savefig(filename,format='png', dpi=100)
    plt.close()


# Unit square domain
def plot_unitSquare_meshes():
    lx = [3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "unitSquare_qu/mesh_l%d.h5" % lx[k]
        filename = output_dir + "unitSquare_qu_mesh_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)


# Gamma-shaped domain
def plot_gammaShaped_meshes():
    lx = [2, 3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "gammaShaped_qu/mesh_l%d.h5" % lx[k]
        filename = output_dir + "gammaShaped_qu_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)

    # deg 1
    lx = [0, 2, 3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "gammaShaped_br/deg1/mesh_l%d.h5" % lx[k]
        filename = output_dir + "gammaShaped_br_deg1_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)

    # deg 2
    lx = [2, 3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "gammaShaped_br/deg2/mesh_l%d.h5" % lx[k]
        filename = output_dir + "gammaShaped_br_deg2_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)

    # deg 3
    lx = [2, 3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "gammaShaped_br/deg3/mesh_l%d.h5" % lx[k]
        filename = output_dir + "gammaShaped_br_deg3_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)


# L-shaped domain
def plot_lShaped_meshes():
    lx = [2, 3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "lShaped_qu/mesh_l%d.h5" % lx[k]
        filename = output_dir + "lShaped_qu_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)

    # deg 1
    lx = [0, 2, 3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "lShaped_br/deg1/mesh_l%d.h5" % lx[k]
        filename = output_dir + "lShaped_br_deg1_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)

    # deg 2
    lx = [2, 3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "lShaped_br/deg2/mesh_l%d.h5" % lx[k]
        filename = output_dir + "lShaped_br_deg2_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)

    # deg 3
    lx = [2, 3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "lShaped_br/deg3/mesh_l%d.h5" % lx[k]
        filename = output_dir + "lShaped_br_deg3_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)


# Square domain, two rectangular sectors
def plot_squareTwoPiecewise_meshes():
    lx = [3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "square_twoPiecewise/qu/mesh_l%d.h5" % lx[k]
        filename = output_dir + "square_twoPiecewise_qu_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)

    lx = [3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "square_twoPiecewise/br/deg2/mesh_l%d.h5" % lx[k]
        filename = output_dir + "square_twoPiecewise_br_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)


# Square domain, four rectangular sectors
def plot_squareFourPiecewise_meshes():
    lx = [3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "square_fourPiecewise/qu/mesh_l%d.h5" % lx[k]
        filename = output_dir + "square_fourPiecewise_qu_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)

    lx = [3]
    for k in range(0, len(lx)):
        meshFile = meshes_dir + "square_fourPiecewise/br/deg2/mesh_l%d.h5" % lx[k]
        filename = output_dir + "square_fourPiecewise_br_l%d.png" % lx[k]
        print(meshFile, "\n", filename)
        mesh = meshio.read_mesh_h5(meshFile)
        save_plot_mesh2d(mesh, filename)


if __name__ == "__main__":
    # plot_unitSquare_meshes()
    # plot_gammaShaped_meshes()
    plot_lShaped_meshes()
    # plot_squareTwoPiecewise_meshes()
    # plot_squareFourPiecewise_meshes()

# End of file
