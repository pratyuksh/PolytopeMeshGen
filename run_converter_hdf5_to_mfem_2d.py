#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

import converters.geometry as gm
import converters.converter_hdf5_to_mfem as mconv
import generators.io_mesh as meshio


# Base input mesh directory
base_in_mesh_dir = "./meshes/"

# Base output mesh directory
base_out_mesh_dir = "./meshes_MFEM/"

# types of bisection-tree refined meshes
br_type = ["deg0", "deg1", "deg2", "deg3"]

# types of graded meshes
gr_type = ["weight3E-1", "weight1E-1", "weight4E-2"]


def call_converter_hdf5_to_mfem_2d(geometry_type, mesh_type, level, deg):
    
    print(geometry_type, mesh_type, level, deg)
    
    geom_file = None
    h5_mesh_file = None
    mfem_mesh_file = None
    
    # Unit-square geometry
    if geometry_type in ['unit-square', 'unit_square', 'unitSquare']:
        
        geom_file = "geometry_files/unitSquare.txt"

        if mesh_type in ['qu', 'quasi-uniform', 'quasi_uniform', 'quasiUniform']:
        
            h5_mesh_file = base_in_mesh_dir+"unitSquare_qu/mesh_l%d.h5"%level
            mfem_mesh_file = base_out_mesh_dir+"unitSquare_qu/mesh_l%d.mesh"%level
        
        else:
            print("Unknown mesh type!")
            sys.exit()
        
    
    # Gamma-shaped geometry
    elif geometry_type in ['gamma-shaped', 'gamma_shaped', 'gammaShaped']:
        
        geom_file = "geometry_files/gammaShaped.txt"
        
        if mesh_type in ['qu', 'quasi-uniform', 'quasi_uniform', 'quasiUniform']:
        
            h5_mesh_file = base_in_mesh_dir+"gammaShaped_qu/mesh_l%d.h5"%level
            mfem_mesh_file = base_out_mesh_dir+"gammaShaped_qu/mesh_l%d.mesh"%level
        
        elif mesh_type in ['br', 'bisection-refined', 'bisection_refined', 'bisectionRefined']:
            
            h5_mesh_file = base_in_mesh_dir+"gammaShaped_br/"+br_type[deg]+"/mesh_l%d.h5"%level
            mfem_mesh_file = base_out_mesh_dir+"gammaShaped_br/"+br_type[deg]+"/mesh_l%d.mesh"%level        
        
        elif mesh_type in ['gr', 'graded']:
            
            h5_mesh_file = base_in_mesh_dir+"gammaShaped_gr/"+gr_type[deg-1]+"/mesh_l%d.h5"%level
            mfem_mesh_file = base_out_mesh_dir+"gammaShaped_gr/"+gr_type[deg-1]+"/mesh_l%d.mesh"%level            
        
        else:
            print("Unknown mesh type!")
            sys.exit()
    
     # L-shaped geometry
    elif geometry_type in ['l-shaped', 'l_shaped', 'lShaped']:
        
        geom_file = "geometry_files/lShaped.txt"
        
        if mesh_type in ['qu', 'quasi-uniform', 'quasi_uniform', 'quasiUniform']:
        
            h5_mesh_file = base_in_mesh_dir+"lShaped_qu/mesh_l%d.h5"%level
            mfem_mesh_file = base_out_mesh_dir+"lShaped_qu/mesh_l%d.mesh"%level
        
        elif mesh_type in ['br', 'bisection-refined', 'bisection_refined', 'bisectionRefined']:
            
            h5_mesh_file = base_in_mesh_dir+"lShaped_br/"+br_type[deg]+"/mesh_l%d.h5"%level
            mfem_mesh_file = base_out_mesh_dir+"lShaped_br/"+br_type[deg]+"/mesh_l%d.mesh"%level        
        
        elif mesh_type in ['gr', 'graded']:
            
            h5_mesh_file = base_in_mesh_dir+"lShaped_gr/"+gr_type[deg-1]+"/mesh_l%d.h5"%level
            mfem_mesh_file = base_out_mesh_dir+"lShaped_gr/"+gr_type[deg-1]+"/mesh_l%d.mesh"%level
        
        else:
            print("Unknown mesh type!")
            sys.exit()
    elif geometry_type in ['square-twoPiecewise', 'square_twoPiecewise', 'squareTwoPiecewise']:
        
        geom_file = "geometry_files/square_twoPiecewise.txt"

        if mesh_type in ['qu', 'quasi-uniform', 'quasi_uniform', 'quasiUniform']:
        
            h5_mesh_file = base_in_mesh_dir+"square_twoPiecewise/qu/mesh_l%d.h5"%level
            mfem_mesh_file = base_out_mesh_dir+"squareTwoPiecewise/qu/mesh_l%d.mesh"%level
        
        elif mesh_type in ['br', 'bisection-refined', 'bisection_refined', 'bisectionRefined']:
            
            h5_mesh_file = base_in_mesh_dir+"square_twoPiecewise/br/"+br_type[deg]+"/mesh_l%d.h5"%level
            mfem_mesh_file = base_out_mesh_dir+"squareTwoPiecewise/br/"+br_type[deg]+"/mesh_l%d.mesh"%level   
        
        else:
            print("Unknown mesh type!")
            sys.exit()
    
    elif geometry_type in ['square-fourPiecewise', 'square_fourPiecewise', 'squareFourPiecewise']:
        
        geom_file = "geometry_files/square_fourPiecewise.txt"

        if mesh_type in ['qu', 'quasi-uniform', 'quasi_uniform', 'quasiUniform']:
        
            h5_mesh_file = base_in_mesh_dir+"square_fourPiecewise/qu/mesh_l%d.h5"%level
            mfem_mesh_file = base_out_mesh_dir+"squareFourPiecewise/qu/mesh_l%d.mesh"%level
        
        elif mesh_type in ['br', 'bisection-refined', 'bisection_refined', 'bisectionRefined']:
            
            h5_mesh_file = base_in_mesh_dir+"square_fourPiecewise/br/"+br_type[deg]+"/mesh_l%d.h5"%level
            mfem_mesh_file = base_out_mesh_dir+"squareFourPiecewise/br/"+br_type[deg]+"/mesh_l%d.mesh"%level
        
        else:
            print("Unknown mesh type!")
            sys.exit()
            
    else:
        
        print("Unknown geometry type!")
        sys.exit()
        
    print(geom_file)
    print(h5_mesh_file)
    print(mfem_mesh_file)
    
    geom = gm.Geometry2d(geom_file)
    ppMesh = mconv.Hdf5ToMfemConverter(geom)
    ppMesh.convert_mesh_hdf5_to_mfem(h5_mesh_file, mfem_mesh_file)

   
if __name__=='__main__':
    
    n = len(sys.argv)
    if n != 5:
        print("Total arguments passed: ", n)
        print("But 5 arguments required!")
        sys.exit()
    
    geometry_type = str(sys.argv[1])
    mesh_type = str(sys.argv[2])
    level = int(sys.argv[3])
    
    deg = int(sys.argv[4])
    if mesh_type in ['qu', 'quasi-uniform', 'quasi_uniform', 'quasiUniform']:
        deg = 1
    
    # call converter
    call_converter_hdf5_to_mfem_2d(geometry_type, mesh_type, level, deg)


# End of file
