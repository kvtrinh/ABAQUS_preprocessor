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
    geomFile = open('robot_simulation.inp','w')
    geomFile.write('*HEADING\n')
    geomFile.write('Robot walking model\n')
    
    # create and write lattice mesh
    ############## lattice size INPUT ######
    mesh = createLatticeMesh(2,2,4,None,0.,0.,0.)
    print('lattice mesh geometry: ',len(mesh.nodeList), ' nodes, ',len(mesh.elemList),' elements')
    # write lattice nodes and elements
    geomFile.write('*NODE, NSET=lattice_nodes\n')
    nsetList = []
    for i in range (len(mesh.nodeList)):
        nodeId = i+1    # nodeId is ABAQUS id
        geomFile.write('%d, %f, %f, %f\n'%(nodeId, mesh.nodeList[i][0],mesh.nodeList[i][1],mesh.nodeList[i][2]))
        nsetList.append(nodeId)  # nodeId is ABAQUS id
        mesh.addNset('lattice_nodes', nsetList)
    mesh.addNset('FEET_NODES',[478, 156, 438, 398, 358, 33, 74, 115])
    mesh.writeNset(geomFile,'FEET_NODES')        
    
    geomFile.write('*ELEMENT, TYPE=B31, ELSET=BEAM\n')
    beamList = []
    for i in range (len(mesh.elemList)):
        geomFile.write('%d, %d, %d\n'%(i+1, mesh.elemList[i][0]+1,mesh.elemList[i][1]+1))
        beamList.append(i+1)
    mesh.addElset('BEAM',beamList)
    # write beam section
    geomFile.write('*BEAM SECTION, SECTION=RECTANGULAR, ELSET=BEAM, MATERIAL=STEEL\n')
    geomFile.write(' 1., 2.\n')

    # actuators and springs (springs are for visual only)
    ############## actuator nodes INPUT ######
    an1 = 378
    an2 = 538
    an3 = 54
    an4 = 218
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
    mesh.addNode(0., -.76, 0.)
    refNodeId = len(mesh.nodeList)
    geomFile.write('*NODE, NSET=ground_ref\n')
    mesh.writeNodeLineLastNode(geomFile)
    ############## table corner nodes INPUT ######
    geomFile.write('*NODE, NSET=RIGID_NODES\n')
    mesh.addNode(-10,-1.1,10)
    gn1 = len(mesh.nodeList)
    mesh.writeNodeLineLastNode(geomFile)
    mesh.addNode(10, -1.1, 10)
    gn2 = len(mesh.nodeList)
    mesh.writeNodeLineLastNode(geomFile)
    mesh.addNode(10, -1.1, -10)
    gn3 = len(mesh.nodeList)
    mesh.writeNodeLineLastNode(geomFile)
    mesh.addNode(-10, -1.1, -10)
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
    geomFile.write('*SURFACE, TYPE=ELEMENT, NAME=GROUND\n')
    geomFile.write('GROUND,SPOS\n')
    geomFile.write('*RIGID BODY, ELSET=GROUND, REFNODE=ground_ref\n')
    geomFile.write('*SURFACE INTERACTION, NAME=INT1\n')
    geomFile.write('*FRICTION\n')
    geomFile.write(' 0.005,\n')
    geomFile.write('*SURFACE BEHAVIOR,PENALTY\n')
    geomFile.write('*CONTACT PAIR,INTERACTION=INT1\n')
    geomFile.write('FEET,GROUND\n')
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
    geomFile.write('*AMPLITUDE,TIME=TOTAL TIME,VALUE=RELATIVE,NAME=MIN2\n')
    geomFile.write('0.,0.,.1, 0.,1.,1.0,2.,0.\n')
    geomFile.write('.,-1.,4.0,0.,5.,1.,6.,0.\n')
    geomFile.write('*MATERIAL, NAME=STEEL\n')
    geomFile.write('*ELASTIC\n')
    geomFile.write(' 3.E6,.3\n')
    geomFile.write('*DENSITY\n')
    geomFile.write(' .0005\n')
    geomFile.write('**\n')
    geomFile.write('*STEP,INC=10000,NLGEOM\n')
    geomFile.write('*DYNAMIC\n')
    geomFile.write('1.0E-4,6.,1.0E-8,.01\n')
    geomFile.write('*DLOAD\n')
    geomFile.write('BEAM,GRAV,9810.,0.,-1.,0.\n')
    geomFile.write('*CONNECTOR MOTION, AMPLITUDE=MIN2, TYPE=DISPLACEMENT\n')
    geomFile.write('ACTUATOR1,1,.5\n')
    geomFile.write('*CONNECTOR MOTION, AMPLITUDE=MIN2, TYPE=DISPLACEMENT\n')
    geomFile.write('ACTUATOR2,1,-.5\n')
    geomFile.write('*END STEP\n')
    geomFile.close()

    print('end createRobotMesh')    
