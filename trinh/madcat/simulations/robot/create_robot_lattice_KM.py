# -*- coding: utf-8 -*-
"""
Create 2x2x4 walking robot structure lattice to get K and M matrices.
Khanh Trinh 
Intelligent Systems Division
NASA Ames Research Center
"""

import sys
sys.path.append("../../src/python");
from  abqiface import *
from  createLatticeMesh import *

print('createRobotMesh')

if __name__ == "__main__":
    
    # print the python version
    print("\n\nUsing python version:")
    print(sys.version)
    print(" \n")

    # set up file
    ############## file name INPUT ######
    geomFile = open('robot_K_M.inp','w')
    geomFile.write('*HEADING\n')
    geomFile.write('Robot walking lattice structure\n')
    
    # create and write lattice mesh
    mesh = ABAQUS_mesh()
    ############## lattice size INPUT ######
    mesh.addVoxelLatticeMesh(2,2,4,None,0.,0.,0.,False)
    print('lattice mesh geometry: ',len(mesh.nodeList), ' nodes, ',len(mesh.elemList),' elements')
    # write lattice nodes and elements
    geomFile.write('*NODE, NSET=lattice_nodes\n')
    nsetList = []
    for i in range (len(mesh.nodeList)):
        nodeId = i+1    # nodeId is ABAQUS id
        geomFile.write('%d, %f, %f, %f\n'%(nodeId, mesh.nodeList[i][0],mesh.nodeList[i][1],mesh.nodeList[i][2]))
        nsetList.append(nodeId)  # nodeId is ABAQUS id
        mesh.addNset('lattice_nodes', nsetList)
    mesh.addNset('FEET_NODES',[366, 407, 75, 448, 117, 489, 159, 33])
    mesh.writeNset(geomFile,'FEET_NODES')        
    
    geomFile.write('*ELEMENT, TYPE=B31, ELSET=BEAM\n')
    beamList = []
    for i in range (len(mesh.elemList)):
        geomFile.write('%d, %d, %d\n'%(i+1, mesh.elemList[i][0]+1,mesh.elemList[i][1]+1))
        beamList.append(i+1)
    mesh.addElset('BEAM',beamList)

    # write voxel beam section set
    ptop = [0,0,1]
    mid1 = [.5,.5,0]
    n1_s1 = [ptop[0]-mid1[0], ptop[1]-mid1[1],ptop[2]-mid1[2]]
    mid2 = [-.5,.5,0]
    n1_s2 = [ptop[0]-mid2[0], ptop[1]-mid2[1],ptop[2]-mid2[2]]
    n1_s3 = n1_s1
    n1_s4 = n1_s2
    py = [0,1,0]
    mid5 = [.5,0,.5]
    n1_s5 = [py[0]-mid5[0], py[1]-mid5[1],py[2]-mid5[2]]
    px = [1,0,0]
    mid6 = [0,.5,.5]
    n1_s6 = [px[0]-mid6[0], px[1]-mid6[1],px[2]-mid6[2]]
    mid7 = [-.5,0,.5]
    n1_s7 = [py[0]-mid7[0], py[1]-mid7[1],py[2]-mid7[2]]
    mid8 = [0,-.5,.5]
    n1_s8 = [px[0]-mid8[0], px[1]-mid8[1],px[2]-mid8[2]]
    n1_s9 = n1_s7
    n1_s10 = n1_s8
    n1_s11 = n1_s5
    n1_s12 = n1_s6
    n1 = [n1_s1,n1_s2,n1_s3,n1_s4,n1_s5,n1_s6,n1_s7,n1_s8,n1_s9,n1_s10,n1_s11,n1_s12]
    for i in range(len(mesh.elsetList['voxelSectionList'])):
        name = 'voxel_strut'+str(i+1)
        n1X = n1[i][0]
        n1Y = n1[i][1]
        n1Z = n1[i][2]
        mesh.addElset(name,mesh.elsetList['voxelSectionList'][i])  # storing ABAAQUS id
        add_one_to_ID = 1
        mesh.writeElset(geomFile,name,add_one_to_ID) 
        geomFile.write('*BEAM SECTION, SECTION=RECTANGULAR, ELSET=%s, MATERIAL=ultem_2200_polyetherimide_E\n'%(name))
        geomFile.write(' .072, .072\n')
        geomFile.write(' %f, %f, %f\n'%(n1X, n1Y, n1Z))

    matlLib = Material_Library()
    matlLib.writeAbaqusMatlLines('ultem_2200_polyetherimide_E', geomFile)
    matlLib.writeAbaqusMatlLines('Aluminum_E', geomFile)
    
    geomFile.write('*STEP\n')
    geomFile.write('*MATRIX GENERATE, STIFFNESS\n')
    geomFile.write('*MATRIX OUTPUT, STIFFNESS\n')    
    geomFile.write('*END STEP\n')
    geomFile.write('*STEP\n')
    geomFile.write('*MATRIX GENERATE, MASS\n')
    geomFile.write('*MATRIX OUTPUT, MASS\n')    
    geomFile.write('*END STEP\n')    
    geomFile.close()

    print('end createRobotMesh')    
