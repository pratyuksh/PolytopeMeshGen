#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import numpy.linalg as la
import fenics


# CLASS: Bisection Mesh Refinement
class BisectionRefinement:

    def __init__(self, polyhedron):

        self.polyhedron = polyhedron

        self.factor = 0.49  # must be < 0.5
        self.R0 = self.calc_R0()

    # converts the geometry corner c_k to fenics.Point(c_k)
    def convert_to_fenics_format(self, z):

        x1 = z[0]
        x2 = z[1]
        x3 = z[2]
        return fenics.Point(x1, x2, x3)

    # calculates R0 for corners
    def calc_R0(self):

        num_corners = len(self.polyhedron.corners)
        R0 = np.zeros(num_corners)
        for j in range(0, num_corners):
            d_min = 1e+10
            for i in range(0, num_corners):
                if i != j:
                    d = la.norm(self.polyhedron.corners[i] - self.polyhedron.corners[j])
                    d_min = min(d_min, d)
            R0[j] = d_min * self.factor

        return R0

    # area of triangle; Heron's formula
    def area_triangle(self, a, b, c):
        p = (a + b + c) / 2
        area = np.sqrt(p * (p - a) * (p - b) * (p - c))
        return area

    # distance of given cell from edge
    def distance_from_edge(self, cell, p1, p2, e):

        vertices = cell.get_vertex_coordinates()
        num_vertices = vertices.shape[0] // 3

        d2e = np.zeros(num_vertices)
        for k in range(num_vertices):
            v = self.convert_to_fenics_format(vertices[k * 3:(k + 1) * 3])
            r1 = fenics.Point.distance(v, p1)
            r2 = fenics.Point.distance(v, p2)
            triArea = self.area_triangle(e, r1, r2)
            d2e[k] = 2 * triArea / e
            # print(k, d2e[k])

        return np.amin(d2e)

    # uniform refinement
    def uniform(self, h0, mesh):

        while mesh.hmax() > h0:
            cell_markers = fenics.MeshFunction("bool", mesh, mesh.topology().dim())
            cell_markers.set_all(True)
            mesh = fenics.refine(mesh, cell_markers)

        return mesh

    # local refinement
    def local(self, deg, h0, mesh):

        dim = 3  # space dimensions

        TOL = fenics.DOLFIN_EPS_LARGE
        singular_edges = self.polyhedron.edges[self.polyhedron.singular_edges - 1, :]
        num_sin_edges = len(singular_edges)

        # hmax_init = mesh.hmax()

        for k in range(num_sin_edges):
            if self.polyhedron.refine_flags[k] == 1:

                v1 = singular_edges[k][0]
                v2 = singular_edges[k][1]

                p1 = self.polyhedron.corners[v1]
                p2 = self.polyhedron.corners[v2]

                R0 = (self.R0[v1] + self.R0[v2]) / 2

                p1f = self.convert_to_fenics_format(p1)
                p2f = self.convert_to_fenics_format(p2)
                e = fenics.Point.distance(p1f, p2f)

                gamma = 1. - self.polyhedron.refine_weights[k]
                K0 = -np.log2(h0) * (2 * deg + dim) / (2 * gamma + dim - 2) - 1
                K = np.ceil(K0)
                NLocRef = int(dim * (K + 1) - 1)
                print(gamma, K, NLocRef)

                for l in range(NLocRef):

                    weight = 1. - gamma
                    expo = -2 * l * (deg + weight) / (dim * (2 * deg + dim))
                    h_min = h0 * np.power(2., expo)
                    R_max = R0 * np.power(2., -l / dim)
                    cell_markers = fenics.MeshFunction("bool", mesh, mesh.topology().dim())
                    cell_markers.set_all(False)

                    count = 0
                    for cell in fenics.cells(mesh):
                        h = fenics.Cell.h(cell)
                        d = self.distance_from_edge(cell, p1f, p2f, e)

                        if d < R_max:
                            if h > h_min:
                                cell_markers[cell] = True
                                count += 1
                                # print("\t", h, d)
                    print('\nRefinement loop: ', l, R_max, h_min, count)
                    mesh = fenics.refine(mesh, cell_markers)

        return mesh


# CLASS: Bisection Mesh Generator 3d
class BisectionMeshGenerator3d:

    def __init__(self, polyhedron):
        self.polyhedron = polyhedron
        self.bis_refine = BisectionRefinement(polyhedron)

        self.corners = None
        self.convert_to_fenics_format()

    # converts the geometry corner c_k to fenics.Point(c_k)
    def convert_to_fenics_format(self):
        num_corners = self.polyhedron.corners.shape[0]
        self.corners = []
        for k in range(0, num_corners):
            x1 = self.polyhedron.corners[k, 0]
            x2 = self.polyhedron.corners[k, 1]
            x3 = self.polyhedron.corners[k, 2]
            self.corners.append(fenics.Point(x1, x2, x3))

    # generates bisection refined mesh
    def bisection(self, deg, mesh_init, h0):
        mesh = self.bisection_uniform(mesh_init, h0)  # Step 1: uniform bisection refinement
        mesh = self.bisection_local(deg, mesh, h0)  # Step 2: local bisection refinement
        return mesh

    # generates _uniform_ bisection refined mesh
    def bisection_uniform(self, mesh_init, h0):
        mesh = self.bis_refine.uniform(h0, mesh_init)
        return mesh

    # generates _locally_ bisection refined mesh
    def bisection_local(self, deg, mesh_init, h0):
        mesh = self.bis_refine.local(deg, h0, mesh_init)
        return mesh

    # writes the Mesh to an output file
    def write(mesh, mesh_file):
        fenics.File(mesh_file) << mesh

# END OF FILE
