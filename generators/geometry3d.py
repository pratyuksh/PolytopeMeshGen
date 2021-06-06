#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import numpy.linalg as la


# CLASS POLYHEDRON #

class Polyhedron:

    def __init__(self):

        self.corners = None
        self.edges = None

        self.singular_corners = None
        self.singular_edges = None
        self.refine_flags = None
        self.refine_weights = None
        self.sin_corners_refine_weights = None
        self.sin_edges_refine_weights = None
        self.sc2se = None

        self.tol = None

    # determines if a given point is on the edge (e1,e2)
    def is_on_edge(self, e1, e2, x):

        d = la.norm(e1 - e2)
        d1 = la.norm(e1 - x)
        d2 = la.norm(e2 - x)
        # print('\t', d, d1, d2)

        val = x[1] * (e2[0] - e1[0]) - x[0] * (e2[1] - e1[1]) - (e1[1] * e2[0] - e2[1] * e1[0])
        if abs(val) <= self.tol and abs(d - (d1 + d2)) <= self.tol * d:
            return True

        return False

    # determines if a given point lies on a specified singular edge
    def point_lies_on_singular_edge(self, sin_edge_id, x):

        edge_id = self.singular_edges[sin_edge_id - 1] - 1

        e1 = self.corners[self.edges[edge_id, 0]]
        e2 = self.corners[self.edges[edge_id, 1]]
        if self.is_on_edge(e1, e2, x):
            return True

        return False

    # determines if a given point is a singular corner
    def is_a_singular_corner(self, x):

        num_singular_corners = self.singular_corners.shape[0]
        if num_singular_corners == 0:
            return False, 0

        for i in range(0, num_singular_corners):
            j = self.singular_corners[i]
            v = self.corners[j - 1]
            if la.norm(v - x) <= self.tol * la.norm(v):
                return True, i + 1

        return False, 0

    # determines if the facet (sin_corner[sin_corner_id], x) lies on a singular edge
    # connected to singular_corner[sin_corner_id]
    def on_singular_edge_connected_to_singular_corner(self, sin_corner_id, x):

        n = len(self.sc2se[sin_corner_id - 1])
        if n == 1:  # no singular edges connected
            if self.sc2se[sin_corner_id - 1][0] == 0:
                return False, 0

        for i in range(0, n):
            edge_id = self.singular_edges[self.sc2se[sin_corner_id - 1][i] - 1]
            e1 = self.corners[self.edges[edge_id - 1, 0]]
            e2 = self.corners[self.edges[edge_id - 1, 1]]

            if self.is_on_edge(e1, e2, x):
                return True, self.sc2se[sin_corner_id - 1][i]

        return False, 0

    # determines if a given point lies on of the singular edges
    def point_lies_on_singular_edges(self, x):

        num_singular_edges = self.singular_edges.shape[0]
        if num_singular_edges == 0:
            return False, 0

        for i in range(0, num_singular_edges):
            sin_edge_id = i + 1
            if self.point_lies_on_singular_edge(sin_edge_id, x):
                return True, sin_edge_id

        return False, 0

    # determines if a given edge lies on one of the singular edges
    def edge_lies_on_singular_edges(self, f):

        num_singular_edges = self.singular_edges.shape[0]
        if num_singular_edges == 0:
            return False, 0

        for i in range(0, num_singular_edges):
            if self.point_lies_on_singular_edge(i, f[0]) \
                    and self.point_lies_on_singular_edge(i, f[1]):
                return True, i + 1

        return False, 0

    # sets singular corners
    def set_singular_corners(self, singular_corners):
        self.singular_corners = singular_corners

    # sets singular edges
    def set_singular_edges(self, singular_edges):
        self.singular_edges = singular_edges

    # sets singular corners to singular edges connectivity
    def set_singular_corner_to_singular_edge_connectivity(self, sc2se):
        self.sc2se = sc2se

    # set refinement weights for corners and singular edges
    def set_refine_weights(self, sin_corners_refine_weights, sing_edges_refine_weights):
        self.sin_corners_refine_weights = np.array(sin_corners_refine_weights)
        self.sin_edges_refine_weights = np.array(sing_edges_refine_weights)
        self.refine_weights = self.sin_corners_refine_weights.tolist() + self.sin_edges_refine_weights.tolist()

    # sets local refinement flags for singular edges
    # used in bisection mesh refinement
    def set_refine_flags(self, refine_flags):
        self.refine_flags = refine_flags


# Cube domain: (0,1)^3
class CubeDomain(Polyhedron):

    def __init__(self):
        self.tol = 1e-15
        self.add_corners()
        self.add_edges()

        self.num_corners = self.corners.shape[0]
        self.num_edges = self.edges.shape[0]

        # default singular corners and edges
        self.add_singular_corners()
        self.add_singular_edges()
        self.add_singular_corner_to_singular_edge_connectivity()

        # default refinement weights
        self.sin_corners_refine_weights = np.array([0.5] * self.singular_corners.shape[0])
        self.sin_edges_refine_weights = np.array([0.5] * self.singular_edges.shape[0])
        self.refine_weights = self.sin_corners_refine_weights.tolist() + self.sin_edges_refine_weights.tolist()

    def add_corners(self):
        self.corners = np.array([[0, 0, 0],
                                 [1, 0, 0],
                                 [1, 1, 0],
                                 [0, 1, 0],
                                 [0, 0, 1],
                                 [1, 0, 1],
                                 [1, 1, 1],
                                 [0, 1, 1]], dtype=np.float64)

    def add_edges(self):
        self.edges = np.array([[0, 1],
                               [1, 2],
                               [2, 3],
                               [0, 3],
                               [4, 5],
                               [5, 6],
                               [6, 7],
                               [4, 7],
                               [0, 4],
                               [1, 5],
                               [2, 6],
                               [3, 7]], dtype=np.int32)

    # indexing starts from 1
    def add_singular_corners(self):
        self.singular_corners = np.array([1, 2])

    # edge indexing starts from 1, index 0 implies no singular edges are connected
    def add_singular_edges(self):
        self.singular_edges = np.array([1], dtype=np.int32)

    # edge indexing starts from 1, index 0 implies no singular edges are connected
    def add_singular_corner_to_singular_edge_connectivity(self):
        self.sc2se = [[1], [1]]


# Prism domain: ((0,1)^2 \ T) \times (0,1); T := ((0,0), (1,0), (0.5,0.5))
class PrismDomain(Polyhedron):

    def __init__(self):
        self.tol = 1e-15
        self.add_corners()
        self.add_edges()

        self.num_corners = self.corners.shape[0]
        self.num_edges = self.edges.shape[0]

        # default singular edges
        self.add_singular_corners()
        self.add_singular_edges()
        self.add_singular_corner_to_singular_edge_connectivity()

        # default refinement weights
        self.sin_corners_refine_weights = np.array([0.5] * self.singular_corners.shape[0])
        self.sin_edges_refine_weights = np.array([0.5] * self.singular_edges.shape[0])
        self.refine_weights = self.sin_corners_refine_weights.tolist() + self.sin_edges_refine_weights.tolist()

    def add_corners(self):
        self.corners = np.array([[0, 0, 0],
                                 [1, 0, 0],
                                 [0.5, 0.5, 0],
                                 [0, 1, 0],
                                 [1, 1, 0],
                                 [0, 0, 1],
                                 [1, 0, 1],
                                 [0.5, 0.5, 1],
                                 [0, 1, 1],
                                 [1, 1, 1]], dtype=np.float64)

    def add_edges(self):
        self.edges = np.array([[0, 2],
                               [1, 2],
                               [1, 3],
                               [3, 4],
                               [0, 4],
                               [5, 7],
                               [6, 7],
                               [6, 8],
                               [8, 9],
                               [5, 9],
                               [0, 5],
                               [1, 6],
                               [2, 7],
                               [3, 8],
                               [4, 9]], dtype=np.int32)

    # indexing starts from 1
    def add_singular_corners(self):
        self.singular_corners = np.array([3, 8])

    # edge indexing starts from 1
    def add_singular_edges(self):
        self.singular_edges = np.array([13], dtype=np.int32)

    # edge indexing starts from 1, index 0 implies no singular edges are connected
    def add_singular_corner_to_singular_edge_connectivity(self):
        self.sc2se = [[1], [1]]


# Fichera Cube domain: [-1,1]^3 \ [0,1]^3
class FicheraCubeDomain(Polyhedron):

    def __init__(self):
        self.tol = 1e-15
        self.add_corners()
        self.add_edges()

        self.num_corners = self.corners.shape[0]
        self.num_edges = self.edges.shape[0]

        # default singular edges
        self.add_singular_corners()
        self.add_singular_edges()
        self.add_singular_corner_to_singular_edge_connectivity()

        # default refinement weights
        self.sin_corners_refine_weights = np.array([0.5] * self.singular_corners.shape[0])
        self.sin_edges_refine_weights = np.array([0.5] * self.singular_edges.shape[0])
        self.refine_weights = self.sin_corners_refine_weights.tolist() + self.sin_edges_refine_weights.tolist()

    def add_corners(self):
        self.corners = np.array([[-1, -1, -1],
                                 [+1, -1, -1],
                                 [+1, +1, -1],
                                 [-1, +1, -1],
                                 [0, 0, 0],
                                 [+1, 0, 0],
                                 [+1, +1, 0],
                                 [0, +1, 0],
                                 [-1, -1, +1],
                                 [+1, -1, +1],
                                 [+1, 0, +1],
                                 [0, 0, +1],
                                 [0, +1, +1],
                                 [-1, +1, +1]], dtype=np.float64)

    def add_edges(self):
        self.edges = np.array([[0, 1],
                               [1, 2],
                               [2, 3],
                               [0, 3],
                               [4, 5],
                               [5, 6],
                               [6, 7],
                               [4, 7],
                               [8, 9],
                               [9, 10],
                               [10, 11],
                               [11, 12],
                               [12, 13],
                               [8, 13],
                               [0, 8],
                               [1, 9],
                               [4, 11],
                               [5, 10],
                               [3, 13],
                               [7, 12],
                               [2, 6]], dtype=np.int32)

    # corner indexing starts from 1
    def add_singular_corners(self):
        self.singular_corners = np.array([5])

    # edge indexing starts from 1
    def add_singular_edges(self):
        self.singular_edges = np.array([5, 8, 17], dtype=np.int32)

    # edge indexing starts from 1, index 0 implies no singular edges are connected
    def add_singular_corner_to_singular_edge_connectivity(self):
        self.sc2se = [[1, 2, 3]]

# End of file
