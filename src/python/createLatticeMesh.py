# -*- coding: utf-8 -*-
"""
Create lattice mesh.
Khanh Trinh 
Intelligent Systems Division
NASA Ames Research Center
"""

import sys
from  abqiface import *

def createLatticeMesh(numX, numY, numZ, elemType, offX, offY, offZ):
    
    mesh = ABAQUS_mesh()
    pitch = 2
    offX = 0.
    offY = 0.
    offZ = 0.
    numBeams = 4
    numVoxels = 0
    numSharedCorners = 0
    for i in range(numX):
        for j in range(numY):
            for k in range(numZ):
                x = offX + pitch*(i)
                y = offY + pitch*(j)
                z = offZ + pitch*(k)
                numVoxels = numVoxels + 1
                print('adding voxel ',numVoxels,x,y,z)
                voxel = Voxel_1(pitch, x, y, z, numBeams)
                # voxel.printVoxelData()
                # add voxel
                sharedNodes = mesh.addVoxel(voxel)
                numSharedCorners = numSharedCorners + sharedNodes
    print('Num shared corners: ', numSharedCorners)         
    return mesh

