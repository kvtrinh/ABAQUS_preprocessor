# -*- coding: utf-8 -*-
"""Test the beam model.

Khanh Trinh
Intelligent Systems Division
NASA Ames Research Center

Version 0: October 17, 2016
"""

import shutil
#import numpy
#import scipy
import sys

#import abaqus
#import abaqusConstants

sys.path.append("../../src/python");
import abqiface


if __name__ == "__main__":
    
    # print the python version
    print("\n\nUsing python version:")
    print(sys.version)
    print(" \n")

    # copy required files from model directory
    shutil.copyfile('./model/geometry.inp','./geometry.inp')
    shutil.copyfile('./model/elements.inp','./elements.inp')
    shutil.copyfile('./model/bcs.inp','./bcs.inp')
    shutil.copyfile('./model/materials.inp','./materials.inp')
    shutil.copyfile('./model/steps.inp','./steps.inp')
#    shutil.copyfile('./model/inputFileOnly.cae','./inputFileOnly.cae')    

    jobName = 'ex1'
    curRun = abqiface.ABAQUS_run(jobName)
    curRun.setGeometry('geometry.inp')
    curRun.setElements('elements.inp')
    curRun.setBCs('bcs.inp')
    curRun.setMaterials('materials.inp')
    curRun.setSteps('steps.inp')

#    curRun.setMatlE('STEEL',30.E6)
    curRun.writeInpFile()

##    # assign the run name from the object (should be 'ex1')
##    run_name = curRun.getInpFileName();
##    print run_name
##
##    # open prebuilt model
##    abaqus.openMdb('inputFileOnly.cae');
##
##    mdb.JobFromInputFile(
##        name = run_name, 
##        inputFileName = run_name + '.inp',
##        type=abaqusConstants.ANALYSIS, 
##        atTime=None, 
##        waitMinutes=0, 
##        waitHours=0, 
##        queue=None, 
##        memory=90, 
##        memoryUnits=abaqusConstants.PERCENTAGE, 
##        getMemoryFromAnalysis=True, 
##        explicitPrecision=abaqusConstants.SINGLE, 
##        nodalOutputPrecision=abaqusConstants.SINGLE, 
##        userSubroutine='', 
##        scratch='', 
##        resultsFormat=abaqusConstants.ODB, 
##        parallelizationMethodExplicit=abaqusConstants.DOMAIN, 
##        numDomains=1, 
##        activateLoadBalancing=False, 
##        multiprocessingMode=abaqusConstants.DEFAULT, 
##        numCpus=1)
##
##
##    ####### 1st run, E = 30.E6
##    mdb.jobs[run_name].submit(consistencyChecking=OFF)
##
##    jobName = run_name; #inputFileNameWithoutInp
##    # read K matrix produced by curRun
##    K_fileName = jobName + '_STIF1.mtx'
##    F_fileName = jobName + '_LOAD2.mtx' 
##    
##    print "Reading load matrix 1"
##    K = abqiface.read_stiff_mtx(K_fileName, 36, ndof_per_node=6)
##    print "Done 1."
##    
##    print 'stiff ness K for E=40E6, nu=0.25'
##    for i in range(36):
##        for j in range(36):
##            print 'K[%d][%d]= %d\n'%(i,j,K[i,j])
##
##    print 'Load vector'
##    F =  abqiface.read_load_mtx(F_fileName, 36, ndof_per_node=6)
##    for i in range(36):
##        print 'F[%d]= %f\n'%(i,F[i])


##    ##### 2nd run, E = 40.E6   ######
##    curRun.setMatlE('STEEL',40.E6)
##    curRun.writeInpFile()
##    ### keep same file names and just overwrite those files
##    mdb.jobs[run_name].submit(consistencyChecking=OFF)
##    
##    # read K matrix produced by curRun
##    K_fileName = jobName + '_STIF1.mtx'
##    F_fileName = jobName + '_LOAD2.mtx' 
##    
##    print "Read load matrix 2"
##    K = abqiface.read_stiff_mtx(K_fileName, 36, 6)
##    print "Done 2"
##    
##    print 'stiffness K for E=40E6, nu=0.25'
##    for i in range(36):
##        for j in range(36):
##            print 'K[%d][%d]= %d\n'%(i,j,K[i,j])
##
##    print 'Load vector'
##    F =  abqiface.read_load_mtx(F_fileName, 36, 6)
##    for i in range(36):
##        print 'F[%d]= %f\n'%(i,F[i])






