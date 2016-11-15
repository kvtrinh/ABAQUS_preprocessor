# -*- coding: utf-8 -*-
"""
Create 2x2x4 walking robot input deck for ABAQUS.
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
    geomFile = open('robot_simulation_4.inp','w')
    geomFile.write('*HEADING\n')
    geomFile.write('Robot walking model\n')
    
    # create and write lattice mesh
    mesh = ABAQUS_mesh()
    ############## lattice size INPUT ######
    numBeamsPerStrut = 4
    mesh.addVoxelLatticeMesh(2,2,4,numBeamsPerStrut,None,0.,0.,0.,True)
    #mesh.addVoxelLatticeMesh(2,2,4,None,0.,0.,0.)
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
        geomFile.write('%d, %d, %d, %d\n'%(i+1, mesh.elemList[i][0]+1,mesh.elemList[i][1]+1,mesh.elemList[i][2]+1))
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
        

    # create feet
    print('adding shell elements for feet')
    feetNodeList = mesh.nsetList['FEET_NODES']
    posZfrictionElemList = []
    negZfrictionElemList = []
    for latticeCornerNodeId in feetNodeList:
        nodeIndex = latticeCornerNodeId-1
        cx = mesh.nodeList[nodeIndex][0]
        cy = mesh.nodeList[nodeIndex][1]
        cz = mesh.nodeList[nodeIndex][2]
        rad = .2
        refinement = 1
        submesh = addHemisphereShells(cx, cy, cz, rad, refinement)
        newNodes = submesh[0]
        newElems = submesh[1]
        cornerNodes = submesh[2]
        nodeShift = len(mesh.nodeList)+1
        elemShift = len(mesh.elemList)+1        
        geomFile.write('*NODE, NSET=surface_feet_nodes\n')
        for i in range(len(newNodes)):
            mesh.addNode([newNodes[i][0],newNodes[i][1],newNodes[i][2]])
            mesh.writeNodeLineLastNode(geomFile)
        geomFile.write('*ELEMENT, TYPE=S4R, ELSET=SHELLS\n')
        
        for i in range(len(newElems)):
            n1 = newElems[i][0] + nodeShift
            n2 = newElems[i][1] + nodeShift
            n3 = newElems[i][2] + nodeShift
            n4 = newElems[i][3] + nodeShift
            mesh.add4Nodes_Element(n1,n2,n3,n4)
            mesh.writeElementLineLastElement(geomFile)


        # creating beams connecting connection node to shell corner nodes
        geomFile.write('*ELEMENT, TYPE=B31, ELSET=BEAM_feet_connectors\n')
        for i in range(len(cornerNodes)):
            mesh.add2Nodes_Element(latticeCornerNodeId, i+ nodeShift)
            mesh.writeElementLineLastElement(geomFile)

        posZfeets = submesh[3]
        negZfeets = submesh[4]
        for elemId in posZfeets:
            if cz <3:
                posZfrictionElemList.append(elemId + 1 + elemShift)
            else:
                negZfrictionElemList.append(elemId + 1 + elemShift)          
        for elemId in negZfeets:
            if cz <3:
                posZfrictionElemList.append(elemId + 1 + elemShift)
            else:
                negZfrictionElemList.append(elemId + 1 + elemShift)             
        
        
    geomFile.write('*SHELL SECTION, ELSET=SHELLS, MATERIAL=Aluminum_E\n')
    geomFile.write(' .1\n')
    geomFile.write('*BEAM SECTION, SECTION=CIRC, ELSET=BEAM_feet_connectors, MATERIAL=Aluminum_E\n')
    geomFile.write(' .2\n')
    mesh.addElset('negZfrictionElems',negZfrictionElemList)
    mesh.addElset('posZfrictionElems',posZfrictionElemList)
    mesh.writeElset(geomFile, 'negZfrictionElems')
    mesh.writeElset(geomFile, 'posZfrictionElems')

    
    # actuators and springs (springs are for visual only)
    ############## actuator nodes INPUT ######
    ##mesh.addNset('FEET_NODES',[478, 156, 438, 398, 358, 33, 74, 115])
    ##mesh.addNset('FEET_NODES',[486, 160, 445, 404, 363, 34, 76, 118])
    an1 = 387
    an2 = 551
    an3 = 55
    an4 = 223
    geomFile.write('*ELEMENT, TYPE=CONN3D2, ELSET=ACTUATOR1\n')
    mesh.add2Nodes_Element(an1, an2)
    mesh.writeElementLineLastElement(geomFile)
    geomFile.write('*ELEMENT, TYPE=CONN3D2, ELSET=ACTUATOR2\n')
    mesh.add2Nodes_Element(an3, an4)
    mesh.writeElementLineLastElement(geomFile)
    geomFile.write('*CONNECTOR SECTION, ELSET=ACTUATOR1\n')
    geomFile.write('axial\n')
    geomFile.write('*CONNECTOR SECTION, ELSET=ACTUATOR2\n')
    geomFile.write('axial\n')    
    geomFile.write('*ELEMENT, TYPE=SPRINGA, ELSET=SPRINGS\n')
    mesh.add2Nodes_Element(an1, an2)
    mesh.writeElementLineLastElement(geomFile)
    mesh.add2Nodes_Element(an3, an4)
    mesh.writeElementLineLastElement(geomFile)
    geomFile.write('*SPRING,ELSET=SPRINGS\n')
    geomFile.write('\n')
    geomFile.write('1.E1,\n')    
    
    # create contact surface for table
    # reference node for table
    ############## reference node INPUT ######
    table_y = -1.71
    mesh.addNode([0., table_y, 0.])
    refNodeId = len(mesh.nodeList)
    geomFile.write('*NODE, NSET=ground_ref\n')
    mesh.writeNodeLineLastNode(geomFile)
    ############## table corner nodes INPUT ######
    geomFile.write('*NODE, NSET=RIGID_NODES\n')
    
    mesh.addNode([-10,table_y,10])
    gn1 = len(mesh.nodeList)
    mesh.writeNodeLineLastNode(geomFile)
    mesh.addNode([10, table_y, 10])
    gn2 = len(mesh.nodeList)
    mesh.writeNodeLineLastNode(geomFile)
    mesh.addNode([10, table_y, -10])
    gn3 = len(mesh.nodeList)
    mesh.writeNodeLineLastNode(geomFile)
    mesh.addNode([-10, table_y, -10])
    gn4 = len(mesh.nodeList)
    mesh.writeNodeLineLastNode(geomFile)
    mesh.add1Node_Element(refNodeId)
    geomFile.write('*ELEMENT, TYPE=MASS, ELSET=REF_MASS\n')
    mesh.writeElementLineLastElement(geomFile)
    mesh.add4Nodes_Element(gn1, gn2, gn3, gn4)
    #geomFile.write('\n')
    geomFile.write('*MASS, ELSET=REF_MASS\n')
    geomFile.write('1000., \n')
    geomFile.write('*ELEMENT, TYPE=R3D4, ELSET=GROUND\n')
    mesh.writeElementLineLastElement(geomFile)
    #geomFile.write('*SURFACE,TYPE=NODE, NAME=FEET\n')
    #geomFile.write('FEET_NODES\n')
    geomFile.write('*SURFACE,TYPE=ELEMENT, NAME=FEET\n')
    geomFile.write('SHELLS\n')
    geomFile.write('*SURFACE,TYPE=ELEMENT, NAME=posZ_FEET\n')
    geomFile.write('posZfrictionElems\n')
    geomFile.write('*SURFACE,TYPE=ELEMENT, NAME=negZ_FEET\n')
    geomFile.write('negZfrictionElems\n')     
    geomFile.write('*SURFACE, TYPE=ELEMENT, NAME=GROUND\n')
    geomFile.write('GROUND,SPOS\n')
    geomFile.write('*RIGID BODY, ELSET=GROUND, REFNODE=ground_ref\n')
    geomFile.write('*SURFACE INTERACTION, NAME=friction_posZ\n')
    geomFile.write('*FRICTION\n')
    geomFile.write(' 0.00001,\n')
    geomFile.write('*SURFACE BEHAVIOR,PENALTY\n')
    geomFile.write('*CONTACT PAIR,INTERACTION=friction_posZ\n')
    geomFile.write('posZ_FEET,GROUND\n')
    geomFile.write('**\n')
    geomFile.write('*SURFACE INTERACTION, NAME=friction_negZ\n')
    geomFile.write('*FRICTION\n')
    geomFile.write(' 0.75,\n')
    geomFile.write('*SURFACE BEHAVIOR,PENALTY\n')
    geomFile.write('*CONTACT PAIR,INTERACTION=friction_negZ\n')
    geomFile.write('negZ_FEET,GROUND\n')
    geomFile.write('**\n')    
##    geomFile.write('\n')
##    geomFile.write('\n')
##    geomFile.write('\n')
##    geomFile.write('\n')

    ############## BC INPUTs ######
    geomFile.write('*BOUNDARY\n')
    geomFile.write('ground_ref,1,6,0.\n')
    geomFile.write('*AMPLITUDE,TIME=TOTAL TIME,VALUE=RELATIVE,NAME=MIN1\n')
    geomFile.write('0.,0.,.2, 0.,1.2,1.0,2.2,0.\n')
    geomFile.write('3.2,1.,4.2,0.,5.2,1.,6.2,1.\n')
    geomFile.write('*AMPLITUDE,TIME=TOTAL TIME,VALUE=RELATIVE,NAME=MIN24\n')
    geomFile.write('0.,0.,.1, 0.,1.,1.0,2.,0.\n')
    geomFile.write('3.,-1.,4.0,0.,5.,1.,6.,0.\n')
    geomFile.write('7.,-1.,8., 0.,9.,1.0,10.,0.\n')
    geomFile.write('11.,-1.,12., 0.,13.,1.0,14.,0.\n')
    geomFile.write('15.,-1.,16., 0.,17.,1.0,18.,0.\n')
    geomFile.write('19.,-1.,20., 0.,21.,1.0,22.,0.\n')
    geomFile.write('23.,-1.,24., 0.,25.,1.0,26.,0.\n')     
    geomFile.write('27.,-1.,28., 0.,29.,1.0,30.,0.\n')

    geomFile.write('31.,-1.,32., 0.,33.,1.0,34.,0.\n')
    geomFile.write('35.,-1.,36., 0.,37.,1.0,38.,0.\n')
    geomFile.write('39.,-1.,40., 0.,41.,1.0,42.,0.\n')
    geomFile.write('43.,-1.,44., 0.,45.,1.0,46.,0.\n')     
    geomFile.write('47.,-1.,48., 0.,49.,1.0,50.,0.\n')
    geomFile.write('51.,-1.,52., 0.,53.,1.0,54.,0.\n')
    geomFile.write('55.,-1.,56., 0.,57.,1.0,58.,0.\n')
    geomFile.write('59.,-1.,60., 0.,61.,1.0,62.,0.\n')



    
    matlLib = Material_Library()
    matlLib.writeAbaqusMatlLines('ultem_2200_polyetherimide_E', geomFile)
    matlLib.writeAbaqusMatlLines('Aluminum_E', geomFile)
##    geomFile.write('*MATERIAL, NAME=Matl_1\n')
##    geomFile.write('*ELASTIC\n')
##    geomFile.write(' 3.E6,.3\n')
##    geomFile.write('*DENSITY\n')
##    geomFile.write(' .0005\n')
##    geomFile.write('**\n')
    geomFile.write('*STEP,INC=40000,NLGEOM\n')
    geomFile.write('*DYNAMIC\n')
    geomFile.write('1.0E-4,60.,1.0E-8,.05\n')    # full run
    geomFile.write('*DLOAD\n')
    geomFile.write('BEAM,GRAV,99810.,0.,-1.,0.\n')
    geomFile.write('*CONNECTOR MOTION, AMPLITUDE=MIN24, TYPE=DISPLACEMENT\n')
    geomFile.write('ACTUATOR1,1,.2\n')
    geomFile.write('*CONNECTOR MOTION, AMPLITUDE=MIN24, TYPE=DISPLACEMENT\n')
    geomFile.write('ACTUATOR2,1,.2\n')
    geomFile.write('*END STEP\n')
    geomFile.close()

    #mesh.printNodeList()

    print('end createRobotMesh')    
