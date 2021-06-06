#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import dolfin as df
import numpy as np
from generators.graded_meshes_base import GradedMeshGeneratorBase


# CLASS: Graded Mesh Generator 2d, triangulation
class GradedMeshGenerator2d(GradedMeshGeneratorBase):

    def __init__(self, polygon):

        self.ndim = 2
        self.polygon = polygon

        self.nvc = 3  # number of vertices per cell
        self.nec = 3  # number of edges per cell
        self.ncc = 4  # number of children per cell

    # marks the edges; 1: interior edge, 2: edge connected to a singular corner
    # generates coordinates of child vertices on the edges
    def create_child_vertices(self, mesh):

        mesh.init(1)
        num_edges = mesh.num_entities(1)
        vertices = mesh.coordinates()

        edge_locality = None
        child_vertices = np.zeros((num_edges, self.ndim))  # array of new vertices
        for k in range(0, num_edges):

            edge_marked = False
            edge = df.Edge(mesh, k)
            edge_vertices_idG = edge.entities(0)
            edge_coords = vertices[edge_vertices_idG]

            for j in range(0, 2):
                corner_flag, corner_id = self.polygon.is_a_singular_corner(edge_coords[j, :])
                # corner_flag, corner_id = self.polygon.is_a_singular_corner(edge_vertices_idG[j])

                if corner_flag:
                    edge_locality = [2, corner_id, j]  # v-type
                    edge_marked = True
                    break

            if not edge_marked:  # interior edge
                edge_locality = [1, 0]
                edge_marked = True

            # default
            weight = 0.5
            v0 = edge_coords[0]
            v1 = edge_coords[1]

            # v-refine
            if edge_locality[0] == 2:
                weight = self.polygon.refine_weights[edge_locality[1] - 1]
                if edge_locality[2] == 1:
                    v0 = edge_coords[1]
                    v1 = edge_coords[0]

            child_vertices[k, :] = v0 * (1 - weight) + v1 * weight
            # print(k, edge_locality, '\n', edge_coords, '\n', child_vertices[k,:], "\n")

        return child_vertices

    # refines all the mesh cells
    def refine(self, mesh):

        nn0 = mesh.num_vertices()
        vertices = mesh.coordinates()
        # cells = mesh.cells()

        child_vertices = self.create_child_vertices(mesh)

        all_cells = np.zeros((self.ncc * mesh.num_cells(), self.nvc), dtype=np.uintp)
        for cell in df.cells(mesh):

            cell_idG = cell.index()
            parent_v = cell.entities(0)

            child_v = np.zeros(self.nec, dtype=np.uintp)
            idL = 0
            for edge in df.edges(cell):
                child_v[idL] = nn0 + edge.index()
                idL += 1

            all_cells[self.ncc * cell_idG: self.ncc * (cell_idG + 1), :] = self.create_child_triangles(parent_v,
                                                                                                       child_v)

        all_vertices = np.zeros((nn0 + child_vertices.shape[0], self.ndim))
        all_vertices[0: nn0, :] = vertices
        all_vertices[nn0:, :] = child_vertices

        return all_vertices, all_cells

# END OF FILE
