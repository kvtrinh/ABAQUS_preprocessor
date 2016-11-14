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

def addHemisphereShells( cx, cy, cz, rad, refinement):
    # create hemisphere
    nodeList = []
    elemList = []    # elemList stores node index
    cornerNodeList = []
    posZfeets = []
    negZfeets = []
     # patch 1
    for i in range(refinement):
        nodeList.append([cx+0, cy - rad, cz + 0])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        nodeList.append([cx+0, cy - rad, cz + -rad/3.])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        nodeList.append([cx+- rad/3., cy - rad, cz + -rad/3.])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        nodeList.append([cx- rad/3., cy - rad, cz + 0.])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        elemList.append([0,3,2,1])
        negZfeets.append(len(elemList) - 1)
    # patch 2
    for i in range(refinement):
        nodeList.append([cx- rad/3., cy - rad, cz + rad/3.])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        nodeList.append([cx+0, cy - rad, cz + rad/3.])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        elemList.append([0,5,4,3])
        negZfeets.append(len(elemList) - 1)
    # patch 3
    for i in range(refinement):
        nodeList.append([cx- 0., cy - 0.3*rad, cz - rad])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        nodeList.append([cx-rad, cy - 0.3*rad, cz - rad])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        elemList.append([1,2,7,6])  
    # patch 4
    for i in range(refinement):
        nodeList.append([cx- rad, cy - 0.3*rad, cz - 0])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        elemList.append([3,8,7,2])  
    # patch 5
    for i in range(refinement):
        nodeList.append([cx- rad, cy - 0.3*rad, cz +rad])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        elemList.append([4,9,8,3])
    # patch 6
    for i in range(refinement):
        nodeList.append([cx- 0., cy - 0.3*rad, cz +rad])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        elemList.append([10,9,4,5])
    # patch 7
    for i in range(refinement):
        nodeList.append([cx+rad/3., cy - rad, cz -rad/3.])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        nodeList.append([cx+rad/3., cy - rad, cz +0])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        elemList.append([0,1,11,12])
        posZfeets.append(len(elemList) - 1)
    # patch 8
    for i in range(refinement):
        nodeList.append([cx+rad/3., cy - rad, cz +rad/3.])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        elemList.append([5,0,12,13])
        posZfeets.append(len(elemList) - 1)
    # patch 9
    for i in range(refinement):
        nodeList.append([cx+rad, cy- 0.3*rad , cz -rad])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        elemList.append([1,6,14,11]) 
    # patch 10
    for i in range(refinement):
        nodeList.append([cx+rad, cy- 0.3*rad , cz ])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        elemList.append([12,11,14,15]) 
    # patch 11
    for i in range(refinement):
        nodeList.append([cx+rad, cy- 0.3*rad , cz+rad ])
        cornerNodeList.append(len(nodeList)-1)
    for i in range(refinement):
        elemList.append([13,12,15,16]) 
    # patch 12    
    for i in range(refinement):
        elemList.append([10,5,13,16]) 
        
  
    return [nodeList, elemList, cornerNodeList, posZfeets, negZfeets] 
