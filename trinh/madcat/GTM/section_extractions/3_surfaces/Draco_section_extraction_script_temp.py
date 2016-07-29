# draco script to extract plane section from surfaces of a part (wing, fuslage, etc.)

####### input ###########
option = 1;
# CAD file in ?
satFile = 'GTM_3wing_surfaces.sat'
baseName = 'GTM_wing_3s'
sectionName = 'GTM_757_wing_3s'
constCoordIndex = 2
h = .25         # element length.  
tol = .025
segmentCoordinates = [11., 24.5, 61.5]   # segment coordinates
segmentIntervals = [10, 10]        # number of sections for each segment

Geometry.load(satFile)
print 'sat file: %s'%satFile
faces = Geometry.faces
print 'faces length:'
print len(faces)
fName = baseName + '.m'
gFile = open(fName, 'w')
gFile.write('% data for aircraft part\n')
gFile.write('%s_data = ['%baseName);

constantSections = []
segmentId = 0
curTraverseCoord = segmentCoordinates[0]
numTotalSections = 1          # number of sections with data based on inputs.  Actual sections created may be less.
for i in range(len(segmentIntervals)):
	numTotalSections  = numTotalSections + segmentIntervals[i]
print 'estimated number of total sections: %d'%numTotalSections

for delta_step in range(numTotalSections):
	print '****** new delta step  ******'
	stepSize = 0.
	if delta_step > 0:
		for j in range(1,len(segmentIntervals)):
			#print 'j %d'%j
			if curTraverseCoord > (segmentCoordinates[j]-tol):
				#print 'curTraverseCoord %f'% curTraverseCoord
				#print 'segmentCoordinates[j] %f'%segmentCoordinates[j]
				segmentId = j
		#print '***segmentId %d'%segmentId
		stepSize = (segmentCoordinates[segmentId+1] - segmentCoordinates[segmentId])/segmentIntervals[segmentId]
		curTraverseCoord = curTraverseCoord + stepSize;

	constantSections.append(curTraverseCoord)
	#print 'step %d'%delta_step
	#print 'step size %f'%stepSize
	print 'constant station value: %f'%curTraverseCoord
	for faceId in faces:
		#print 'face %s'%faceId
		if isinstance(faceId, str):
			continue
		if faceId.name == '_body11':
			#continue
			pass
		#print 'boundingBox'
		bd = faceId.boundingBox
		if constCoordIndex == 1:
			if bd.start.x > curTraverseCoord:
				continue
			if bd.end.x < curTraverseCoord:
				continue
		if constCoordIndex == 2:
			if bd.start.y > curTraverseCoord:
				continue
			if bd.end.y < curTraverseCoord:
				continue
		if constCoordIndex == 3:
			if bd.start.z > curTraverseCoord:
				continue
			if bd.end.z < curTraverseCoord:
				continue
		#print 'Sheet.FromFace'
		sh = Sheet.FromFace(faceId)
		####### input ###########
		minNonConst1 = -200.
		minNonConst2 = -200.
		maxNonConst1 = 800.
		maxNonConst2 = 800.
		if constCoordIndex == 1:
			plane = Sheet.Box(PositionXYZ(curTraverseCoord,minNonConst1, minNonConst2),PositionXYZ(curTraverseCoord,maxNonConst1,maxNonConst2))
		if constCoordIndex == 2:
			plane = Sheet.Box(PositionXYZ(minNonConst1, curTraverseCoord, minNonConst2),PositionXYZ(maxNonConst1, curTraverseCoord, maxNonConst2))
		if constCoordIndex == 3:
			plane = Sheet.Box(PositionXYZ(minNonConst1, minNonConst2, curTraverseCoord),PositionXYZ(maxNonConst1,maxNonConst2,curTraverseCoord))

        #print 'Bool'
		wire = Bool.NonregIntersection(sh,plane)
		plane.delete()
		print 'number of edges %d' %len(wire.edges)
		for e in wire.edges:
			l = e.length
			if l < h:
				continue;
			numPoints = 99     #int(round(l/h))
			if numPoints < 1:
				continue
		#numPoints = 1
		#print 'numPoints %d'%numPoints
			du = 1./(numPoints)
			for node in range(0,numPoints):
				u = node*du
				p =  e.getPosition(u)
			#print p
				gFile.write('%f %f %f\n'%(p.x, p.y, p.z))
		wire.delete()
		sh.delete()
    
gFile.write('];\n');
gFile.write('const_stations = [\n');
for cs in constantSections:
	gFile.write('%f\n'%(cs))
gFile.write('];\n');
gFile.close()

