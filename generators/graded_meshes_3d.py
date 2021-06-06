#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import dolfin as df
import numpy as np
from generators.graded_meshes_base import GradedMeshGeneratorBase


# CLASS: Anisotropic Graded Meshes 3d, tetrahedral elements
class GradedMeshGenerator3d(GradedMeshGeneratorBase):

    def __init__(self, polyhedron):

        self.ndim = 3
        self.polyhedron = polyhedron

        self.nvc = 4  # number of vertices per cell
        self.nec = 6  # number of edges per cell
        self.ncc = 8  # number of children per cell

        self.nsc = polyhedron.singular_corners.shape[0]  # no. of singular corners
        self.refine_weights = polyhedron.refine_weights

    # get weight from edge_locality
    def get_weight(self, edge_locality):
        weight = 0.5  # default
        if edge_locality[0] > 1:
            weight = self.polyhedron.refine_weights[edge_locality[1] - 1]
        return weight

    # get child vertex from edge_locality
    def get_child_vertex(self, edge_locality, edge_coords):

        weight = self.get_weight(edge_locality)
        v0 = edge_coords[0]
        v1 = edge_coords[1]

        if edge_locality[0] == 2 \
                or edge_locality[0] == 3:
            if edge_locality[2] == 1:
                v0 = edge_coords[1]
                v1 = edge_coords[0]

        if edge_locality[0] == 4:
            if edge_locality[3] == 1:
                v0 = edge_coords[1]
                v1 = edge_coords[0]

        return v0 * (1 - weight) + v1 * weight

    # marks the edges
    # 1: interior edge or singular edge  
    # 2: edge connected to a singular corner
    # 3: edge connected to a singular edge through a vertex
    # 4: edge is singular + an edge vertex is a singular
    # generates coordinates of child vertices on the edges
    def create_child_vertices(self, mesh, bool_quasiuniform):

        mesh.init(1)
        num_edges = mesh.num_entities(1)
        vertices = mesh.coordinates()

        child_vertices = np.zeros((num_edges, self.ndim))  # array of child vertices

        # quasi-uniform refinement
        if bool_quasiuniform:
            print("\n\tQuasi-uniform mesh refinement running...")
            for k in range(0, num_edges):
                edge = df.Edge(mesh, k)
                edge_vertices_id_global = edge.entities(0)
                edge_coords = vertices[edge_vertices_id_global]
                child_vertices[k, :] = 0.5 * (edge_coords[0, :] + edge_coords[1, :])
            return child_vertices

        # anisotropic refinement
        print("\n\tAnisotropic mesh refinement running...")
        edges_locality = []
        # mark edges, compute child vertices for interior-, v_e- , e-, ev- type edges
        for k in range(0, num_edges):

            edge_marked = False
            edge = df.Edge(mesh, k)
            edge_vertices_id_global = edge.entities(0)
            edge_coords = vertices[edge_vertices_id_global]

            for j in range(0, 2):
                sin_corner_flag, sin_corner_id = self.polyhedron.is_a_singular_corner(edge_coords[j, :])

                if sin_corner_flag:

                    # sin_edge_flag = False
                    v = None
                    if j == 0:
                        v = edge_coords[1, :]
                    elif j == 1:
                        v = edge_coords[0, :]

                    sin_edge_flag, sin_edge_id = self.polyhedron.on_singular_edge_connected_to_singular_corner(
                        sin_corner_id, v)

                    if sin_edge_flag:
                        edge_locality = [4, sin_corner_id, sin_edge_id + self.nsc, j]  # ev-type
                    else:
                        edge_locality = [2, sin_corner_id, j]  # v-type

                    edge_marked = True
                    edges_locality.append(edge_locality)
                    break

            if not edge_marked:

                flag1, sin_edge_id = self.polyhedron.edge_lies_on_singular_edges(edge_coords)
                if flag1:
                    edge_locality = [1, 0]  # e-type
                    edge_marked = True
                    edges_locality.append(edge_locality)
                else:
                    for m in range(0, 2):
                        flag2, sin_edge_id = self.polyhedron.point_lies_on_singular_edges(edge_coords[m])
                        if flag2:
                            edge_locality = [3, sin_edge_id + self.nsc, m]  # v_e-type
                            edge_marked = True
                            edges_locality.append(edge_locality)
                            break

            if not edge_marked:  # interior edge
                edge_locality = [1, 0]
                edge_marked = True
                edges_locality.append(edge_locality)

            child_vertices[k, :] = self.get_child_vertex(edges_locality[k], edge_coords)

        # loop over v-type elements and choose minimum weights, if an ev-type edge is present in the connected cells.
        mesh.init(0, 1)  # Build connectivity between vertices and edges
        for k in range(0, num_edges):

            if edges_locality[k][0] == 2:

                edge = df.Edge(mesh, k)
                edge_vertices_id_global = edge.entities(0)
                edge_coords = vertices[edge_vertices_id_global]

                v0_id_local = edges_locality[k][2]
                v0_id_global = edge.entities(0)[v0_id_local]

                v0 = edge_coords[0]
                v1 = edge_coords[1]
                if v0_id_local == 1:
                    v0 = edge_coords[1]
                    v1 = edge_coords[0]

                # search if an edge connected to v0 is ev-type
                # if yes, then weight = min(kappa_v0, kappa_e)
                kappa_v0 = self.get_weight(edges_locality[k])
                weight = kappa_v0  # default

                # gather the connected edges
                vertex = df.Vertex(mesh, v0_id_global)
                connected_edges_id_global = vertex.entities(1)
                connected_edges_id_global = connected_edges_id_global.tolist()
                connected_edges_id_global.remove(k)
                # print(k, v0_id_global, connected_edges_id_global)

                # search for ev-type connected edges
                for j in range(0, len(connected_edges_id_global)):
                    connected_edge = connected_edges_id_global[j]
                    if edges_locality[connected_edge][0] == 4:
                        connected_edge_locality = edges_locality[connected_edge]
                        kappa_e = self.polyhedron.refine_weights[connected_edge_locality[2] - 1]
                        weight = min(kappa_v0, kappa_e)
                        # print("\t", connected_edge, kappa_v0, kappa_e, weight)
                        break

                child_vertices[k, :] = v0 * (1 - weight) + v1 * weight

            # edge = df.Edge(mesh, k)
            # edge_vertices_id_global = edge.entities(0)
            # edge_coords = vertices[edge_vertices_id_global]
            # print(k, edge_locality, '\n', edge_coords, '\n', child_vertices[k,:], "\n")

        return child_vertices

    # refines all the mesh cells
    # uses ufc numbering and info in cells_locality
    def refine(self, mesh, bool_quasiuniform):

        nn0 = mesh.num_vertices()
        vertices = mesh.coordinates()

        child_vertices = self.create_child_vertices(mesh, bool_quasiuniform)

        all_cells = np.zeros((self.ncc * mesh.num_cells(), self.nvc), dtype=np.uintp)
        for cell in df.cells(mesh):

            cell_id_global = cell.index()
            parent_v = cell.entities(0)
            parent_coords = vertices[parent_v]

            child_v = np.zeros(self.nec, dtype=np.uintp)
            child_coords = np.zeros((self.nec, self.ndim))
            id_local = 0
            for edge in df.edges(cell):
                child_v[id_local] = nn0 + edge.index()
                child_coords[id_local, :] = child_vertices[edge.index()]
                id_local += 1

            all_cells[self.ncc * cell_id_global: self.ncc * (cell_id_global + 1), :] \
                = self.create_child_tetrahedrons(parent_v, parent_coords, child_v, child_coords)

        all_vertices = np.zeros((nn0 + child_vertices.shape[0], self.ndim))
        all_vertices[0: nn0, :] = vertices
        all_vertices[nn0:, :] = child_vertices

        return all_vertices, all_cells

# End of file
