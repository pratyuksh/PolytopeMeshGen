#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import dolfin as df
import numpy as np
import generators.tetrahedron_properties as tetprop


# BASE CLASS: Graded Mesh Generator
class GradedMeshGeneratorBase:

    def __init__(self):

        self.ndim = None
        self.ncc = None
        self.nvc = None

    # creates new mesh given vertices and cells using FEniCs MeshEditor
    def create(self, vertices, cells):

        editor = df.MeshEditor()
        mesh = df.Mesh()

        celltype = None
        if self.ndim == 2:
            celltype = "triangle"
        elif self.ndim == 3:
            celltype = "tetrahedron"

        editor.open(mesh, celltype, self.ndim, self.ndim, 1)
        editor.init_vertices(vertices.shape[0])  # number of vertices
        editor.init_cells(cells.shape[0])  # number of cells

        for k in range(0, vertices.shape[0]):
            editor.add_vertex(k, vertices[k, :])

        for k in range(0, cells.shape[0]):
            editor.add_cell(k, cells[k, :])

        editor.close()
        return mesh

    # creates 4 child triangles from a parent triangle
    # uses ufc numbering
    def create_child_triangles(self, parent_v, child_v):

        child_triangles = np.zeros((self.ncc, self.nvc), dtype=np.uintp)

        # child 1
        child_triangles[0, 0] = parent_v[0]
        child_triangles[0, 1] = child_v[1]
        child_triangles[0, 2] = child_v[2]

        # child 2
        child_triangles[1, 0] = parent_v[1]
        child_triangles[1, 1] = child_v[2]
        child_triangles[1, 2] = child_v[0]

        # child 3
        child_triangles[2, 0] = parent_v[2]
        child_triangles[2, 1] = child_v[0]
        child_triangles[2, 2] = child_v[1]

        # child 4
        child_triangles[3, 0] = child_v[0]
        child_triangles[3, 1] = child_v[1]
        child_triangles[3, 2] = child_v[2]

        return child_triangles

    # creates 8 child tetrahedrons from the parent tetrahedrons
    # uses ufc numbering
    def create_child_tetrahedrons(self, parent_v, parent_coords,
                                  child_v, child_coords):  # , cell_vertices, new_vertices_coords, new_vertices_idG):

        child_tetrahedrons = np.zeros((self.ncc, self.nvc), dtype=np.uintp)

        # 4 child tetrahedrons connected to the vertices of the parent tetrahedron
        child_tetrahedrons[0, :] = np.array([parent_v[0], child_v[5], child_v[4], child_v[3]], dtype=np.uintp)
        child_tetrahedrons[1, :] = np.array([parent_v[1], child_v[5], child_v[2], child_v[1]], dtype=np.uintp)
        child_tetrahedrons[2, :] = np.array([parent_v[2], child_v[4], child_v[2], child_v[0]], dtype=np.uintp)
        child_tetrahedrons[3, :] = np.array([parent_v[3], child_v[3], child_v[1], child_v[0]], dtype=np.uintp)

        # partitioning the octagon to 8 new tetrahedrons
        tet_prop = tetprop.TetrahedronProperties()

        child_type1 = np.zeros((4, self.nvc), dtype=np.uintp)  # common edge: (x_{01}, x_{23})
        child_type1[0, :] = np.array([5, 4, 3, 0], dtype=np.uintp)
        child_type1[1, :] = np.array([5, 4, 2, 0], dtype=np.uintp)
        child_type1[2, :] = np.array([5, 1, 3, 0], dtype=np.uintp)
        child_type1[3, :] = np.array([5, 1, 2, 0], dtype=np.uintp)

        child_type2 = np.zeros((4, self.nvc), dtype=np.uintp)  # common edge: (x_{02}, x_{13})
        child_type2[0, :] = np.array([5, 4, 3, 1], dtype=np.uintp)
        child_type2[1, :] = np.array([5, 4, 2, 1], dtype=np.uintp)
        child_type2[2, :] = np.array([0, 4, 3, 1], dtype=np.uintp)
        child_type2[3, :] = np.array([0, 4, 2, 1], dtype=np.uintp)

        child_type3 = np.zeros((4, self.nvc), dtype=np.uintp)  # common edge: (x_{03}, x_{12})
        child_type3[0, :] = np.array([5, 4, 3, 2], dtype=np.uintp)
        child_type3[1, :] = np.array([5, 1, 3, 2], dtype=np.uintp)
        child_type3[2, :] = np.array([0, 4, 3, 2], dtype=np.uintp)
        child_type3[3, :] = np.array([0, 1, 3, 2], dtype=np.uintp)

        aspect_ratios = np.zeros((3, 4))
        for i in range(0, 4):
            aspect_ratios[0, i] = tet_prop.aspect_ratio(child_coords[child_type1[i]])
            aspect_ratios[1, i] = tet_prop.aspect_ratio(child_coords[child_type2[i]])
            aspect_ratios[2, i] = tet_prop.aspect_ratio(child_coords[child_type3[i]])

        # vec_ar = np.array([aspect_ratios[0].min(), aspect_ratios[1].min(), aspect_ratios[2].min()])
        vec_ar = np.amin(aspect_ratios, axis=1)
        # print(aspect_ratios, '\n', vec_ar, '\n\n')

        child_type = np.copy(child_type1)
        if vec_ar[1] == vec_ar.max():
            child_type = np.copy(child_type2)
        elif vec_ar[2] == vec_ar.max():
            child_type = np.copy(child_type3)

        for i in range(0, 4):
            child_tetrahedrons[i + 4, :] = child_v[child_type[i]]

        return child_tetrahedrons

# END OF FILE
