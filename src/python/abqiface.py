# -*- coding: utf-8 -*-
"""
Abaqus interface module.
Khanh Trinh & Stefan Schuet
Intelligent Systems Division
NASA Ames Research Center
"""

import os
import re
import shutil
#import numpy
#import scipy.sparse

print('Hi Khanh')

class ABAQUS_run:  
    """ 
    Class to generate ABAQUS input file from partitioned ABAQUS input file.
    """
    
    numInpOutputs = 0
    def __init__(self, inpFileName):
        self.inpFileName = inpFileName       # string
        self.geom = None                     # text file, part of ABAQUS input deck
        self.elements = None                 # text file, part of ABAQUS input deck
        self.materials = None                # text file, part of ABAQUS input deck
        self.interactions = None
        self.bcs = None
        self.loads = None
        self.steps = None
        self.staticStep = None
        self.writeStiffnessStep = None
        self.writeLoadStep = None
        self.flag_writeSteps = 0         
        self.flag_writeStatic = 0
        self.flag_writeStiffness = 0
        self.flag_writeLoad = 0
      
    def setGeometry(self, fileName):
        self.geom = fileName
        
    def setElements(self, fileName):
        self.elements = fileName
        
    def setBCs(self, fileName):
        self.bcs = fileName
        
    def setLoads(self, fileName):
        self.loads = fileName
        
    def setMaterials(self, fileName):
        self.materials = fileName

    def setInteractions(self, fileName):
        self.interactions = fileName
        
    def setSteps(self, fileName):
        self.steps = fileName
        
    def setStaticStep(self, fileName):
        self.staticStep = fileName
        
    def setWriteStiffnessStep(self, fileName):
        self.writeStiffnessStep = fileName
        
    def setWriteLoadStep(self, fileName):
        self.writeLoadStep = fileName
        

    def getInpFileName(self):
        #curDir = os.getcwd()
        inpFileName = self.inpFileName+'.inp'
        print("inpFileName:  %s" % inpFileName)
        return self.inpFileName

    def setMatlE(self,matl,E):
        oldMatlFile =  open(self.materials)
        tempFileName = 'tempMaterials.inp'
        newMatlFile = open(tempFileName,'w')
        foundMatl = False
        replaceEline = -1
        numLine = 0
        for line in oldMatlFile:
            numLine = numLine + 1
            if numLine == replaceEline:
                values = line.split(',')
                #print values
                #print(len(values)
                replaceEline = -1
                # values[1] included a '\n'
                newLine = str(E)+', '+values[1] 
                newMatlFile.write(newLine)
            else:  
                if matl in line:
                    foundMatl = True
                    replaceEline = numLine + 2
                newMatlFile.write(line)
        if foundMatl:
            print('Found material name %s\n'%matl)
            shutil.move(tempFileName, self.materials)
        

    def setMatlNu(self,matl,nu):
        oldMatlFile =  open(self.materials)
        tempFileName = 'tempMaterials.inp'
        newMatlFile = open(tempFileName,'w')
        foundMatl = False
        replaceEline = -1
        numLine = 0
        for line in oldMatlFile:
            numLine = numLine + 1
            if numLine == replaceEline:
                values = line.split(',')
                #print(values
                #print(len(values)
                replaceEline = -1
                newLine = values[0]+', '+str(nu) +os.linesep
                newMatlFile.write(newLine)
            else:  
                if matl in line:
                    foundMatl = True
                    replaceEline = numLine + 2
                newMatlFile.write(line)
        if foundMatl:
            print('Found material name %s\n'%matl)
            shutil.move(tempFileName, self.materials)

    def modifySurfacePressure(self, loadFile, surfaceName, pressureValue):
        # assume loads are defined in a step
        # assume there is a "*Dsload" line
        oldLoadFile =  open(loadFile)
        tempFileName = 'tempLoadFile.inp'
        newLoadFile = open(tempFileName,'w')
        dsLoadlineSection = False
        dsLoadlineNumber = -1
        modifiedPressure = False
        numLine = 0
        changedFile = False
        for line in oldLoadFile:
            newLine = line
            numLine = numLine + 1
            if "Dsload" in line:
                dsLoadlineSection = True
                dsLoadlineNumber = numLine
                modifiedPressure = False

            if dsLoadlineSection:
                print(line)
                if numLine == dsLoadlineNumber:
                    pass
                else:
                    if "*" in line:
                        print('* in line')
                        dsLoadlineSection = False
                        if modifiedPressure == False:
                            extraLine = surfaceName + ', P, ' + str(pressureValue) +os.linesep
                            newLoadFile.write(extraLine)
                            print('*** printing newline cause surf not found ' + extraLine)
                            changedFile = True
                            
                    if surfaceName in line:                    
                        newLine = surfaceName + ', P, ' + str(pressureValue) +os.linesep
                        modifiedPressure = True
                        print('change line to ' + newLine)
                        changedFile = True
            newLoadFile.write(newLine)
            
        if changedFile:
            print('Set pressure on surface %s to %f\n'%(surfaceName, pressureValue))
            shutil.move(tempFileName, loadFile)        


    def setTractionSeparationKs(self,interaction,Knn_Kss_Ktt):
        Knn = Knn_Kss_Ktt[0]
        Kss = Knn_Kss_Ktt[1]
        Ktt = Knn_Kss_Ktt[2]
        oldInteractFile =  open(self.interactions)
        tempFileName = 'tempInteractions.inp'
        newInteractFile = open(tempFileName,'w')
        foundInteract = False
        replaceEline = -1
        numLine = 0
        for line in oldInteractFile:
            numLine = numLine + 1
            if numLine == replaceEline:
                values = line.split(',')
                print(values)
                print(len(values))
                replaceEline = -1
                # values[1] included a '\n'
                newLine = str(Knn)+', '+str(Kss) + ', ' + str(Ktt) + os.linesep
                print('change line to ' + newLine)
                newInteractFile.write(newLine)
                #newInteractFile.write(line)
            else:  
                if interaction in line:
                    if '*Surface Interaction' in line:
                        foundInteract = True
                        replaceEline = numLine + 3
                newInteractFile.write(line)
        if foundInteract:
            print('Found interaction name %s\n'%interaction)
            shutil.move(tempFileName, self.interactions)
        
    def setSurfaceInteractionThickness(self,interaction,thickness):
        oldInteractFile =  open(self.interactions)
        tempFileName = 'tempInteractions.inp'
        newInteractFile = open(tempFileName,'w')
        foundInteract = False
        replaceEline = -1
        numLine = 0
        for line in oldInteractFile:
            numLine = numLine + 1
            if numLine == replaceEline:
                values = line.split(',')
                print(values)
                print(len(values))
                replaceEline = -1
                newLine = str(str(thickness) + os.linesep)
                print('change line to ' + newLine)
                newInteractFile.write(newLine)
                #newInteractFile.write(line)
            else:  
                if interaction in line:
                    if '*Surface Interaction' in line:
                        foundInteract = True
                        replaceEline = numLine + 1
                newInteractFile.write(line)
        if foundInteract:
            print('Found interaction name %s\n'%interaction)
            shutil.move(tempFileName, self.interactions)

    def addSteps(self):
        print('add  steps to run input')
        self.flag_writeSteps = 1

    def addStaticStep(self):
        print('add static step to run input')
        self.flag_writeStatic = 1

    def addWriteSteps(self):
        print('add write load step to run input')
        self.flag_writeSteps = 1        

    def addWriteStiffnessStep(self):
        print('add write stiffness step to run input')
        self.flag_writeStiffness = 1

    def addWriteLoadStep(self):
        print('add write load step to run input')
        self.flag_writeLoad = 1

    def clearStepFlags(self):
        self.flag_writeSteps = 0
        self.flag_writeStatic = 0
        self.flag_writeStiffness = 0
        self.flag_writeLoad = 0
            
    def writeInpFile(self):
        print('write input file %s'% self.inpFileName)
        fileName = self.inpFileName + '.inp'
        inpFile = open(fileName,'w')

        # mesh geometry
        if self.geom != None:
            geomFile = open(self.geom)
            str = geomFile.read()
            inpFile.write(str)
            geomFile.close()

        # element definitions
        if self.elements != None:
            elemFile = open(self.elements)
            str = elemFile.read()
            inpFile.write(str)
            elemFile.close()

        # materials definitions
        if self.materials != None:
            matlFile = open(self.materials)
            str = matlFile.read()
            inpFile.write(str)
            matlFile.close()

        # bcs definitions
        if self.bcs != None:
            bcsFile = open(self.bcs)
            str = bcsFile.read()
            inpFile.write(str)
            bcsFile.close()

        # interactions definitions
        if self.interactions != None:
            interactionsFile = open(self.interactions)
            str = interactionsFile.read()
            inpFile.write(str)
            interactionsFile.close()
            
        # steps definitions
        numberOfSteps = 0
        if self.flag_writeSteps == 1:
            print('write step file name: %s' %self.steps)
            if self.steps != None:
                stepFile = open(self.steps)
                str = stepFile.read()
                inpFile.write(str)
                #print str
                numberOfSteps = numberOfSteps + 1
                stepFile.close()

        if self.flag_writeStatic== 1:
            print('write load step file name: %s' %self.staticStep)
            if self.staticStep != None:
                stepFile = open(self.staticStep)
                str = stepFile.read()
                inpFile.write(str)
                #print str
                numberOfSteps = numberOfSteps + 1
                stepFile.close()

        if self.flag_writeStiffness == 1:
            print('write load step file name: %s' %self.writeStiffnessStep)
            if self.writeStiffnessStep != None:
                stepFile = open(self.writeStiffnessStep)
                str = stepFile.read()
                inpFile.write(str)
                #print str
                numberOfSteps = numberOfSteps + 1
                stepFile.close()

        if self.flag_writeLoad == 1:
            print('write load step file name: %s' %self.writeLoadStep)
            if self.writeLoadStep != None:
                stepFile = open(self.writeLoadStep)
                str = stepFile.read()
                inpFile.write(str)
                #print str
                numberOfSteps = numberOfSteps + 1
                stepFile.close() 

        if numberOfSteps < 1:
            print('Error in ABAQUS input file: no step is defined')
                

        inpFile.close()

#
# End ABAQUS_run
#

def writeLoadVector(fileName, unitLoadFileName, p1, p2, p3):
    inpFile = open(unitLoadFileName,'r')
    numLine = 0
    newLoadFile = open(fileName, 'w')
    for line in inpFile:
        numLine = numLine + 1            
        values = line.split(',')
        if len(values) != 3:
            newLoadFile.write(line)
            continue
        if line.startswith('*'):
            newLoadFile.write(line)
            continue
        dof = int(values[1])
        newLoad = float(values[2])
        if dof == 1:
            newLoad = float(values[2])*p1
        if dof == 2:
            newLoad = float(values[2])*p2
        if dof == 3:
            newLoad = float(values[2])*p3            
                            
        newLine = values[0]+', '+values[1]+', '+str(newLoad) +os.linesep
        newLoadFile.write(newLine)
    inpFile.close()
    newLoadFile.close()

       
def read_stiff_mtx(filename, ndof, ndof_per_node=3, output_sparse=False):
    """Read abaqus matrix (mtx) file and return the matrix."""

    f = open(filename, 'r')

    # run through lines in file and build the i,j,val lists
    rind = []
    cind = []
    val  = []
    
    for line in f:
        vals = line.split(",")
        rn = int(vals[0])-1
        rdof = int(vals[1])-1
        cn = int(vals[2])-1
        cdof = int(vals[3])-1
        
        # Ignore negative row/col indices.
        if rn < 0:
            continue;
        if cn < 0:
            continue;
        
        row = rn*ndof_per_node + rdof
        col = cn*ndof_per_node + cdof
        
        #if row<0 or col<0:
            #continue
        
        val.append(float(vals[4]))
        rind.append(row)
        cind.append(col)
        
        if row != col:
            val.append(float(vals[4]))
            rind.append(col)
            cind.append(row)
    
    f.close()
    
    #k = scipy.sparse.coo_matrix((val,(rind,cind)), shape=(ndof,ndof))
    k = scipy.sparse.csc_matrix((val,(rind,cind)), shape=(ndof,ndof))
    
    if not output_sparse:
        k = numpy.array(k.todense())
        
    return k

    
    
def read_load_mtx(filename, ndof, ndof_per_node=3):
    """Read abaqus matrix (mtx) file and return the load vector."""

    f = open(filename, 'r')
    load = numpy.zeros(shape=ndof)

    numLine = 0
    for line in f:
        numLine = numLine + 1
        if numLine < 3:
            continue
        vals = line.split(",")
        rn = int(vals[0]) - 1
        rdof = int(vals[1]) - 1
        val = float(vals[2])
        rowGlobaId = rn*ndof_per_node + rdof
        load[rowGlobaId] = val
        
    f.close()

    return load



def read_displacement_vector(filename, ndof, ndof_per_node=3):
    """Return displacement vector read from file."""

    fid  = open(filename, 'r')
    u    = numpy.zeros(shape=(ndof))
    
    # Build regular expression (regex) to match solution header
    header_pattern = re.compile("\s*NODE\s*FOOT\-\s*U1\s*U2\s*U3\s*")
    footer_pattern = re.compile("\s*MAXIMUM")
    
    # Skip everything before the header pattern
    for line in fid:
        if header_pattern.match(line) is not None:
            break

    numLine = 0
    
    for line in fid:

        if footer_pattern.match(line) is not None:
            break;

        parts = line.split()
        
        # the block below should handle the file formatting
        try:
            #node = int(parts[0].strip()) - 1
            #u1   = float(parts[1].strip())
            #u2   = float(parts[2].strip())
            #u3   = float(parts[3].strip())
            
            #u[node*3 + 0] = u1
            #u[node*3 + 1] = u2
            #u[node*3 + 2] = u3
            
            # General version that should work when ndof_per_node not equal
            # to 3, but this is untested.
            node = int(parts[0].strip()) - 1
            for ind in range(ndof_per_node):
                u[node*ndof_per_node + ind] = float(parts[ind+1].strip())
            
        except ValueError:
            continue
        except IndexError:
            continue
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise
        
           
    fid.close()

    return u
