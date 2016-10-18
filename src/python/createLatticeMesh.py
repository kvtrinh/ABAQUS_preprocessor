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


sys.path.append("../../src/python");
import abqiface

print('createRobotMesh')


if __name__ == "__main__":
    
    # print the python version
    print("\n\nUsing python version:")
    print(sys.version)
    print(" \n")

    mesh = createLatticeMesh(2,2,4,None,0.,0.,0.)
    print('final mesh geometry: ',len(mesh.nodeList), ' nodes, ',len(mesh.elemList),' elements')
    geomFile = open('geometry.inp','w')
    geomFile.write('*HEADING\n')
    geomFile.write('Robot walking model\n')
    geomFile.write('*NODE, NSET=allnodes\n')
    for i in range (len(mesh.nodeList)):
        nodeId = i+1
        geomFile.write('%d, %f, %f, %f\n'%(nodeId, mesh.nodeList[i][0],mesh.nodeList[i][1],mesh.nodeList[i][2]))

    geomFile.write('*ELEMENT, TYPE=B31, ELSET=BEAM\n')
    for i in range (len(mesh.elemList)):
        geomFile.write('%d, %d, %d\n'%(i+1, mesh.elemList[i][0]+1,mesh.elemList[i][1]+1))    
    geomFile.close()
    #mesh.printNodeList()
    #mesh.printElemList()
    print('end createRobotMesh')    
