#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import numpy.linalg as la


### CLASS GEOMETRY 2D###

class Geometry2d():

    def __init__(self, geom_file):
        
        self.tol = 1e-15
        self.parse(geom_file)
        
    
    def parse(self, geom_file):
    
        # open geometry.txt file to read
        f = open(geom_file, "r")
  
        # read mesh info
        line = f.readline() # skip comment: mesh dimension
        line = f.readline()
        self.ndim = int(line[0])
        
        # read vertex coordinates
        while True: # skip newlines and comments
            line = f.readline()
            if not (line.startswith('#') or line.startswith('\n')):
                break
        
        self.num_vertices = int(line[0])
        self.vertices = np.zeros((self.num_vertices, self.ndim))
        for k in range(0,self.num_vertices):
            line = f.readline().split(',')
            for i in range(0, self.ndim):
                self.vertices[k,i] = float(line[i])
        
        # read boundary facets
        while True: # skip newlines and comments
            line = f.readline()
            if not (line.startswith('#') or line.startswith('\n')):
                break
        
        self.num_bdry_facets = int(line[0])
        self.bdry_facets = np.zeros((self.num_bdry_facets, 2), dtype=np.int32)
        for k in range(0,self.num_bdry_facets):
            line = f.readline().split(',')
            for i in range(0,2):
                self.bdry_facets[k,i] = int(line[i])
    
    
    # determines the attribute index for a given mesh boundary facet
    def get_bdry_attrib_id(self,be):
        
        for i in range(0,self.num_bdry_facets):
            bf_id = i+1
            bf1 = self.vertices[self.bdry_facets[i,0]-1, :]
            bf2 = self.vertices[self.bdry_facets[i,1]-1, :]
            
            if self.is_on_bdry_edge(bf1, bf2, be[0]) and self.is_on_bdry_edge(bf1, bf2, be[1]):
                return bf_id
        
        return 0
    
    
    # determines if a given facet is on the boundary
    def is_bdry_facet(self,e):
        
        for i in range(0,self.num_bdry_facets):
            bf_id = i+1
            bf1 = self.vertices[self.bdry_facets[i,0]-1, :]
            bf2 = self.vertices[self.bdry_facets[i,1]-1, :]
            
            if self.is_on_bdry_edge(bf1, bf2, e[0]) and self.is_on_bdry_edge(bf1, bf2, e[1]):
                return True, bf_id
        
        return False, bf_id
        
    
    # determines if a given point is on the boundary
    def is_on_bdry_edge(self, e1, e2, z):
        
        val =  z[1]*(e2[0]-e1[0]) - z[0]*(e2[1]-e1[1]) - (e1[1]*e2[0]-e2[1]*e1[0])
        if abs(val) <= self.tol and abs(la.norm(z-e1) + la.norm(z-e2) - la.norm(e1-e2)) <= self.tol:
            return True
        
        return False
    
    
    # determines if a given point is on the boundary
    def is_on_bdry(self,z):
        
        for i in range(0,self.num_bdry_facets):
            bf1 = self.vertices[self.bdry_facets[i,0], :]
            bf2 = self.vertices[self.bdry_facets[i,1], :]
            if self.is_on_bdry_edge(e1, e2, z):
                return True
        
        return False



class Geometry3d():

    def __init__(self, geom_file):
        
        self.tol = 1e-15
        self.parse(geom_file)
        
    
    def parse(self, geom_file):
    
        # open geometry.txt file to read
        f = open(geom_file, "r")
  
        # read mesh info
        line = f.readline() # skip comment: mesh dimension
        line = f.readline()
        self.ndim = int(line[0])
        
        # read vertex coordinates
        while True: # skip newlines and comments
            line = f.readline()
            if not (line.startswith('#') or line.startswith('\n')):
                break
        
        self.num_vertices = int(line)
        self.vertices = np.zeros((self.num_vertices, self.ndim))
        for k in range(0,self.num_vertices):
            line = f.readline().split(',')
            for i in range(0, self.ndim):
                self.vertices[k,i] = float(line[i])
        
        # read boundary facets
        while True: # skip newlines and comments
            line = f.readline()
            if not (line.startswith('#') or line.startswith('\n')):
                break
        
        self.num_bdry_facets = int(line)
        self.bdry_facets = np.zeros((self.num_bdry_facets, self.ndim), dtype=np.int32)
        for k in range(0,self.num_bdry_facets):
            line = f.readline().split(',')
            for i in range(0,self.ndim):
                self.bdry_facets[k,i] = int(line[i])

    # determines the attribute index for a given mesh boundary facet
    def get_bdry_attrib_id(self,be):
    
        for i in range(0,self.num_bdry_facets):
            bf_id = i+1
            el = self.vertices[self.bdry_facets[i,:]-1, :]
            #print("\n\n", i, "\n", el)
            if self.is_on_bdry_facet(el, be[0]) and self.is_on_bdry_facet(el, be[1]) and self.is_on_bdry_facet(el, be[2]):
                return bf_id
        
        return 0

    # determines if a given point is on the boundary
    def is_on_bdry_facet(self, el, P):
        
        v0 = el[1,:] - el[0,:]
        v1 = el[2,:] - el[0,:]
        v2 = P - el[0,:]
        
        n1 = np.cross(v0, v1)
        n2 = np.cross(v2, v1)
        
        #print("\n\n", P)
        #print(np.abs(np.dot(n1,n2)), la.norm(n1)*la.norm(n2), np.abs(np.abs(np.dot(n1,n2)) - la.norm(n1)*la.norm(n2)) > self.tol)
        
        # point P not in plane
        if ( np.abs(np.abs(np.dot(n1,n2)) - la.norm(n1)*la.norm(n2)) > self.tol ):
            return False
        
        # if point P is in plane, check if in triangle
        dotv00 = np.dot(v0,v0)
        dotv11 = np.dot(v1,v1)
        dotv01 = np.dot(v0,v1)
        dotv20 = np.dot(v2,v0)
        dotv21 = np.dot(v2,v1)
    
        s = (dotv11*dotv20 - dotv01*dotv21)/(dotv00*dotv11 - dotv01*dotv01)
        t = (dotv00*dotv21 - dotv01*dotv20)/(dotv00*dotv11 - dotv01*dotv01)
        #print(dotv00, dotv11, dotv01, dotv20, dotv21)
        #print(s, t)
        if ((s+self.tol >= 0) and (t+self.tol >= 0) and ((1-s-t+self.tol) >= 0)):
            return True
        else:
            return False
        
        #print( self.sameSide(P, el[2,:], el[0,:], el[1,:]) )
        #print( self.sameSide(P, el[0,:], el[1,:], el[2,:]) )
        #print( self.sameSide(P, el[1,:], el[2,:], el[0,:]) )

    """
    def sameSide (self, p1,p2, a,b):

        cp1 = np.cross(b-a, p1-a)
        cp2 = np.cross(b-a, p2-a)
        if (np.dot(cp1, cp2) >= 0):
            return True
        else:
            return False
    """
    
    
# main function calls
if __name__=='__main__':
    
    geom_file = "geometry_files/square_domain.txt"
    square_domain = Geometry2d(geom_file)
    
    
# END OF FILE
