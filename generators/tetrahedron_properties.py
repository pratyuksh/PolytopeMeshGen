#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import numpy.linalg as la


class TetrahedronProperties:

    def __init__(self):
        self.ndim = 3
        self.nv = 4

    ##
    def eval_alpha(self, tet):

        mat = np.ones((self.nv, self.nv))
        for i in range(0, self.nv):
            for j in range(0, self.ndim):
                mat[i, j] = tet[i, j]

        return la.det(mat)

    # volume of a tetrahedron
    def volume(self, tet):
        return np.abs(self.eval_alpha(tet)) / 6

    # surface area of a tetrahedron
    # noinspection PyMethodMayBeStatic
    def surface_area(self, tet):

        areas = np.zeros(4)
        areas[0] = 0.5 * la.norm(np.cross(tet[2] - tet[1], tet[3] - tet[1]))
        areas[1] = 0.5 * la.norm(np.cross(tet[2] - tet[0], tet[3] - tet[0]))
        areas[2] = 0.5 * la.norm(np.cross(tet[1] - tet[0], tet[3] - tet[0]))
        areas[3] = 0.5 * la.norm(np.cross(tet[1] - tet[0], tet[2] - tet[0]))

        return np.sum(areas)

    # radius of the sphere inscribed in a tetrahedron
    # noinspection PyMethodMayBeStatic
    def inradius(self, tet):
        return 3 * self.volume(tet) / self.surface_area(tet)

    # radius of the sphere circumscribing a tetrahedron
    # noinspection PyMethodMayBeStatic
    def circumradius(self, tet):

        alpha = self.eval_alpha(tet)

        mat = np.ones((self.nv, self.nv))
        for i in range(0, self.nv):
            for j in range(0, self.ndim):
                mat[i, j + 1] = tet[i, j]
            mat[i, 0] = np.sum(mat[i, 1:4] ** 2)
        gamma = la.det(mat)

        matx = np.ones((self.nv, self.nv))
        matx[:, 0] = mat[:, 0]
        matx[:, 1] = mat[:, 2]
        matx[:, 2] = mat[:, 3]
        Dx = +la.det(matx)

        maty = np.ones((self.nv, self.nv))
        maty[:, 0] = mat[:, 0]
        maty[:, 1] = mat[:, 1]
        maty[:, 2] = mat[:, 3]
        Dy = -la.det(maty)

        matz = np.ones((self.nv, self.nv))
        matz[:, 0] = mat[:, 0]
        matz[:, 1] = mat[:, 1]
        matz[:, 2] = mat[:, 2]
        Dz = +la.det(matz)

        return np.sqrt(np.abs(Dx ** 2 + Dy ** 2 + Dz ** 2 - 4 * alpha * gamma)) / (2 * np.abs(alpha))

    # aspect ratio of a tetrahedron
    def aspect_ratio(self, tet):
        return 3 * self.inradius(tet) / self.circumradius(tet)

    def skewness(self, tet):

        # given tet
        cell_size = self.volume(tet)
        R = self.circumradius(tet)

        # equilateral tet
        a = R * np.sqrt(8. / 3.)
        optimal_cell_size = a ** 3 / (6 * np.sqrt(2))

        return (optimal_cell_size - cell_size) / optimal_cell_size

# End of file
