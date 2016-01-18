import abaqus
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import numpy
import math
# import scipy.optimize 
#import math
#Define input parameter variables
baL=76.2
p=1.0/4.0
p2=1.0/16.0
rb=1.5875
R1=355.6
R2=355.6
R3=355.6
rf1=19.05
rf2=19.05
rf3=6.35
rtw=12.7
Tb1=11.90625
Tb2=19.05
h1=Tb1+Tb2
h2=76.5125
h3=93.6625
h4=46.83125
h5=85.725
xb1=0.0
yb1=0.0
wb=76.2
H=h1+h3+h4
baH=Tb1
wtop=74.6125/2.0
wm=8.334375
RL=1066.8 #Length of the rail
lin=266.7 #horizontal distance between inclusion center to near end of the rail
rin=1.5 #Inclusion radius
lh=15 #number inclusion locations at horizontal distance in the rail head
lvlim=12 # number of inclusion locations at vertical distance in the rail.
#Arc center points
Xc1=xb1-wm-R1
Yc1=yb1+h5
Xc2=xb1-wm-R2
Yc2=yb1+h5
Xc3=xb1
Yc3=yb1+H-R3
ct1=(Xc1,Yc1)
ct2=(Xc2,Yc2)
ct3=(Xc3,Yc3)
Tt1=wtop*p
h=h1+h3+Tt1
htop=(h4-Tt1)*p2
#force application parameters
Lf=21.6 #Length of the force application surface
#wf=76.2 #width of the force application surface
#Area_f=101.6*76.2
#Force=263000.0/8*4.4482
PLoad=1360.7
a=9.2
b=7.8
Etie=500.0
Eincl=200000.0 #basically no inclusion
Erail=200000.0
#Pressure1=Force/Area_f
Pressure2=-2.0
L_of1=RL-lin-Lf/2
L_of2=RL-lin+Lf/2
#global mesh size
Elsize=6.0
minf=0.10
devf=0.20
def makeInclusion(myModel,rad,xc,yc,zc,RailPart):
#rad=radius of inclusion
#xc=center of inclusion from x axis
#yc=center of inclusion from y axis
#zc=center of inclusion of z axis
#Design of inclusion gap
	# RailPart.PartitionFaceByDatumPlane(datumPlane=
		# RailPart.datums[8], faces=RailPart.faces.findAt(((xc+1.0,yc,0.0), ), ))
	SketchInG=myModel.ConstrainedSketch(gridSpacing=1.0, name='InclusionGap', 
		sheetSize=2400.0, transform=
		RailPart.MakeSketchTransform(
		sketchPlane=RailPart.datums[9], 
		sketchPlaneSide=SIDE1, 
		sketchUpEdge=RailPart.edges.findAt((35.410084, 164.283968, 187.325), ), 
		sketchOrientation=TOP, origin=(xc, 0.0, 0.0)))
	SketchInG.sketchOptions.setValues(decimalPlaces=6)
	SketchInG.sketchOptions.setValues(gridAuto=OFF, gridSpacing=1.0)
	SketchInG.sketchOptions.setValues(numMinorGridIntervals=1)
	RailPart.projectReferencesOntoSketch(filter=
		COPLANAR_EDGES, sketch=SketchInG)
	# Construction of an arc with a center and two points
	InArc=SketchInG.ArcByCenterEnds(center=(-zc,yc), 
		direction=CLOCKWISE, point1=(-zc, yc-rad), point2=(-zc,yc+rad))
	InLine=SketchInG.Line(point1=(-zc, yc-rad), point2=(-zc, yc+rad))
	LineCo=SketchInG.geometry.findAt((-zc,yc+rad/2))
	SketchInG.VerticalConstraint(addUndoState=False, entity=LineCo)	
	##Vertical Construction Line
	CoLine1=SketchInG.ConstructionLine(point1=(
		-zc, 60.0), point2=(-zc, -60.0))
	#Select Entity for vertical Constraint	
	VertEnt=SketchInG.geometry.findAt((-zc, 0.0))
	#Apply vertical constraint
	SketchInG.VerticalConstraint(addUndoState=False, entity=VertEnt)
	##Horizontal Construction Line (not needed now)
		# CoLine2=SketchInG.ConstructionLine(point1=(
			# -850.0, yc), point2=(-750.0, yc))
		# HorEnt=SkecthInG.geometry.findAt((-800.0, yc))
		# SketchInG.VerticalConstraint(addUndoState=False, entity=HorEnt)
	
	#Turn on the construction line
	SketchInG.sketchOptions.setValues(constructionGeometry=ON)
	
	#Assign the centreline for revolving
	#VertEnt=SketchInG.geometry.findAt((-zc, 0.0))
	SketchInG.assignCenterline(line=VertEnt)
	
	#Revolve the spherical cut
	RailPart.CutRevolve(angle=360.0, flipRevolveDirection=OFF, sketch=
    SketchInG, sketchOrientation=TOP, 
    sketchPlane=RailPart.datums[9],  sketchPlaneSide=SIDE1, sketchUpEdge=
    RailPart.edges.findAt((35.410084, 164.283968, 187.325), ))
	del SketchInG
	#Design of inclusion	
	sketch2=myModel.ConstrainedSketch(name='inclusionsection', sheetSize=400.0)	
	sketch2.ArcByCenterEnds(center=(xb1, yc),
		direction=COUNTERCLOCKWISE, point1=(xb1+rad, yc), point2=(xb1-rad, yc))
	sketch2.Line(point1=(xb1-rad, yc), point2=(xb1+rad, yc))
	sketch2.ConstructionLine(point1=(xb1, yc), point2=(xb1+20.0, yc))
	sketch2.sketchOptions.setValues(constructionGeometry=ON)
	sketch2.assignCenterline(line=sketch2.geometry.findAt((6.0,yc),))
	partInclusion=myModel.Part(dimensionality=THREE_D, name='inclusion', type=
		DEFORMABLE_BODY)
	partInclusion.BaseSolidRevolve(angle=360.0, 
		flipRevolveDirection=OFF, sketch=sketch2)	
	# #Partition the cut surface
	# RailPart.PartitionFaceByDatumPlane(datumPlane=RailPart.datums[6], faces=
		# RailPart.faces.findAt(((xc,yc,zc+rad), ), ) )
	# RailPart.PartitionFaceByDatumPlane(datumPlane=RailPart.datums[9], faces=
		# RailPart.faces.findAt(((xc,yc,zc-rad), ), ))
	# RailPart.PartitionFaceByDatumPlane(datumPlane=RailPart.datums[9], faces=
		# RailPart.faces.findAt(((xc,yc,zc+rad), ), ))	
	# RailPart.PartitionFaceByDatumPlane(datumPlane=RailPart.datums[8], faces=
		# RailPart.faces.findAt(((xc,yc,zc+rad), ), ))

g_incr=open('inpnumb.txt','w')
pnumb=[]
for l in range(lh):
	if (l<((lh-1)/2+1)):
		xci1=4.0*l
		x_a=abs(xci1)
	elif (l>=((lh-1)/2+1)):
		xci1=-4.0*(l-(lh-1)/2)
		x_a=abs(xci1)
	if (x_a <(wm-rin)):
		lv=1
	else:
		lv=1
	for k in range(lv):
		#length of the rail
		dintp=10.32+k*10 #vertical distance between inclusion center and the top of the rail center
		din=yb1+H-dintp #vertical distance between inclusion center and the bottom of the rail center
		#Translation of inclusion
		Trin=RL-lin #Horizontal distance between inclusion center to the far end of the rail
		##Inclusion center: 
		yci1=din
		zci1=RL-lin
		p=lvlim*l+k
		# #Find the intersection of the rail surface arc and top rail side face
		# def f1(variables1):
			# (x,y)=variables1
			
			# #Equation of the arc
			# ArcEq=(x-xb1)**2+(y-yb1-H+R3)**2-R3**2
			
			# #Equation of the line
			# LineEq1=y*(h4-Tt1)*p2-x*(H-h1-h3-Tt1)-(yb1+h1+h3+Tt1)*(h4-Tt1)*p2-(xb1-wtop)*(H-h1-h3-Tt1)
			# return[ArcEq,LineEq1]
		# intersec1=scipy.optimize.fsolve(f1,(15,170)) 
		# #Find the intersection of the curving and inclined face on the shoulder of the rail
		# def f2(variables2):
			# (x,y)=variables2
			# Arc2Eq=(x-xb1-wm-R2)**2+(y-yb1-h5)**2-R2**2
			# Line2Eq=y*wtop-x*(-Tt1)-(yb1+h1+h3+Tt1)*wtop+(xb1-wtop)*(-Tt1)
			# return[Arc2Eq,Line2Eq]
		# intersec2=scipy.optimize.fsolve(f2,(13,90))
		# #Find the intersection of the curving and inclined face on foot of the rail
		# def f3(variables3):
			# (x,y)=variables3
			# Arc3Eq=(x-xb1-wm-R1)**2+(y-yb1-h5)**2-R1**2
			# Line3Eq=y*wb-x*h1-yb1*wb+(xb1-wtop)*h1
			# return[Arc3Eq, Line3Eq]
		# intersec3=scipy.optimize.fsolve(f3,(13,90))	
		# #Define vertices
		vtces=((xb1,yb1),(xb1,yb1+H),(xb1-wtop,yb1+h1+h3+Tt1),(xb1-wm,yb1+h5),
			(xb1-wb,yb1+Tb1),(xb1-wb,yb1))
		intersec1=(-35.07055819,169.71638149)
		intersec2=(-10.77456866,127.31239217)
		intersec3=(-13.1036804,27.6803299)
		#Start sketching
		myModel=mdb.Model(name='SimTrack'+str(p))
		sketch1=myModel.ConstrainedSketch(name='Rail section', sheetSize=400.0)
		sketch1.sketchOptions.setValues(decimalPlaces=6)
		sketch1.sketchOptions.setValues(gridAuto=OFF, gridSpacing=1.0)
		sketch1.sketchOptions.setValues(numMinorGridIntervals=1)
		sketch1.Line(point1=vtces[0],point2=vtces[1])
		sketch1.ArcByCenterEnds(center=ct3, direction=COUNTERCLOCKWISE, 
			point1=vtces[1],point2=intersec1)
		sketch1.Line(point1=intersec1, point2=vtces[2])
		sketch1.Line(point1=vtces[2],point2=intersec2)
		sketch1.ArcByCenterEnds(center=ct2, direction=COUNTERCLOCKWISE, 
			point1=vtces[3],point2=intersec2)
		sketch1.ArcByCenterEnds(center=ct1, direction=CLOCKWISE, 
			point1=vtces[3], point2=intersec3)
		sketch1.Line(point1=intersec3,point2=vtces[4])
		sketch1.Line(point1=vtces[4],point2=vtces[5])
		sketch1.Line(point1=vtces[5],point2=vtces[0])
		#Filleting
		#points to pick when chosing entities
		#Picking the top circular curve, "3rd fillet"
		xfc3=0.98*intersec1[0]
		yfc3=Yc3+(R3**2-(xfc3-Xc3)**2)**0.5
		xfl3=1.02*intersec1[0]
		yfl3=p2*(xfl3-intersec1[0])+intersec1[1]

		#Picking points for the second fillet
		xfc2=0.98*intersec2[0]
		yfc2=Yc2+(R2**2-(xfc2-Xc2)**2)**0.5
		xfl2=1.02*intersec2[0]
		yfl2=-p*(xfl2-intersec2[0])+intersec2[1]
		# #Picking points for 1st smaller fillet:
		xf2=0.99*vtces[2][0]
		yf2=p2*(xf2-intersec1[0])+intersec1[1]
		yf22=-p*(xf2-vtces[2][0])+vtces[2][1]

		#Picking points for the third fillet
		xfc1=0.9*intersec3[0]
		yfc1=Yc1-(R1**2-(xfc1-Xc1)**2)**0.5
		xfl1=1.02*intersec3[0]
		yfl1=p*(xfl1-intersec3[0])+intersec3[1]

		#Picking points for forth fillet
		xf4=0.98*(vtces[4])[0]
		yf4=p*(xf4-(vtces[4])[0])+(vtces[4])[1]
		xfl4=(vtces[4])[0]
		yfl4=0.98*(vtces[4])[1]

		#Picking points for fifth fillet
		xf5=xb1-wb
		yf5=(yb1+1)*0.02
		xfl5=0.98*(vtces[4])[0]
		yfl5=yb1

		#Start Filleting
		sketch1.FilletByRadius(curve1=
			sketch1.geometry[3], curve2=sketch1.geometry[4], nearPoint1=
			(xfc3,yfc3), nearPoint2=(xfl3,yfl3), radius=rf3)
		sketch1.FilletByRadius(curve1=
			sketch1.geometry[4], curve2=
			sketch1.geometry[5], nearPoint1=(
			xf2,yf2), nearPoint2=(xf2,yf22), radius=rb)
		sketch1.FilletByRadius(curve1=
			sketch1.geometry[6], curve2=
			sketch1.geometry[5], nearPoint1=(
			xfc2,yfc2), nearPoint2=(xfl2,yfl2), radius=rtw)	
		sketch1.FilletByRadius(curve1=
			sketch1.geometry[7], curve2=
			sketch1.geometry[8], nearPoint1=(
			xfc1,yfc1), nearPoint2=(xfl1,yfl1), radius=rf1)	
		sketch1.FilletByRadius(curve1=
			sketch1.geometry[8], curve2=
			sketch1.geometry[9], nearPoint1=(
			xf4,yf4), nearPoint2=(xfl4,yfl4), radius=rb)	
		sketch1.FilletByRadius(curve1=
			sketch1.geometry[9], curve2=
			sketch1.geometry[10], nearPoint1=(
			xf5,yf5), nearPoint2=(xfl5,yfl5), radius=rb)		
		myPart=myModel.Part(name='RailLeft',dimensionality=THREE_D,
			type=DEFORMABLE_BODY)
		myPart.BaseSolidExtrude(depth=RL, sketch=sketch1)

		RailPart=myModel.Part(compressFeatureList=ON, mirrorPlane=YZPLANE, name=
			'RailRight'+str(p), objectToCopy=myModel.parts['RailLeft'])
		RailPart.Mirror(keepOriginal=ON, mirrorPlane=RailPart.faces.findAt((
			0.0, 57.149999, 711.200033), ))
		del sketch1
		#Make Datum Planes
		x2b=5.92
		RailPart.DatumPlaneByPrincipalPlane(offset=
		   L_of1, principalPlane=XYPLANE) #Datums[3]
		RailPart.DatumPlaneByPrincipalPlane(offset= 
			L_of2, principalPlane=XYPLANE)#Datum[4]
		RailPart.DatumPlaneByPrincipalPlane(offset=
		  x2b, principalPlane=YZPLANE) #Datum[5]	
		RailPart.DatumPlaneByPrincipalPlane(offset=
		  xb1, principalPlane=YZPLANE) #Datum[6]
		RailPart.DatumPlaneByPrincipalPlane(offset=
		  Trin, principalPlane=XYPLANE) #Datum[7] 
		yp=40.5
		RailPart.DatumPlaneByPrincipalPlane(offset=
		  yp, principalPlane=XZPLANE) #Datum[8]
		RailPart.DatumPlaneByPrincipalPlane(offset=xci1,  principalPlane=YZPLANE) # Datum [9] 
		RailPart.DatumPlaneByPrincipalPlane(offset=H,  principalPlane=XZPLANE) # Datum [10] 
		##Partition
		RailPart.PartitionCellByDatumPlane(cells=
			RailPart.cells.findAt(((xb1, H, RL-lin), ), ), datumPlane=
			RailPart.datums[3])
		RailPart.PartitionCellByDatumPlane(cells=
			RailPart.cells.findAt(((xb1, H, RL-lin), ), ), datumPlane=
			RailPart.datums[4]) 
		RailPart.PartitionFaceByDatumPlane(datumPlane=
			RailPart.datums[6], faces=
				RailPart.faces.findAt(((xb1, H, RL-lin), ), ))
		RailPart.PartitionCellByDatumPlane(cells=
			RailPart.cells.findAt(((xb1, H, RL-lin), ), ), datumPlane=
				RailPart.datums[8])
		#make inclusion gap and inclusion part	
		makeInclusion(myModel,rin,xci1,yci1,zci1,RailPart)	
		partInclusion=myModel.parts['inclusion']
		#Partition the cut surface
		RailPart.PartitionFaceByDatumPlane(datumPlane=RailPart.datums[9], faces=
			RailPart.faces.findAt(((xci1,yci1,zci1+rin), ), ) )
		RailPart.PartitionFaceByDatumPlane(datumPlane=RailPart.datums[7], faces=
			RailPart.faces.findAt(((xci1-rin,yci1,zci1), ), ))
		RailPart.PartitionFaceByDatumPlane(datumPlane=RailPart.datums[7], faces=
			RailPart.faces.findAt(((xci1+rin,yci1,zci1), ), ))
		#Create partition for vertical load application:
		xce=-14.8
		zce=zci1
		SkeLo=myModel.ConstrainedSketch(gridSpacing=1.0, name='ellipticLoad', 
			sheetSize=2400.00, transform=RailPart.MakeSketchTransform( sketchPlane=RailPart.datums[10], 
			sketchPlaneSide=SIDE1,sketchUpEdge=RailPart.edges.findAt((-10.0, 0.0, 0.0), ), 
			sketchOrientation=RIGHT, origin=(0.0, 171.45, 0.0)))
		SkeLo.sketchOptions.setValues(decimalPlaces=6)
		SkeLo.sketchOptions.setValues(gridAuto=OFF, gridSpacing=1.0)
		SkeLo.sketchOptions.setValues(numMinorGridIntervals=1)
		RailPart.projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=SkeLo)
		SkeLo.EllipseByCenterPerimeter(axisPoint1=(zce, xce+b),axisPoint2=(zce+a,xce),center=(zce,xce))	
		RailPart.PartitionFaceBySketchThruAll(faces=RailPart.faces.findAt(((-3.292674,171.434755, 803.700012), )), sketch=
		SkeLo, sketchPlane= RailPart.datums[10], sketchPlaneSide=  SIDE1, sketchUpEdge= RailPart.edges.findAt((-10.0, 0.0, 
			0.0), ))
			
		## Start sketching the ties
		l1=114.30 #side tie width
		l2=228.6 #middle tie width
		z1=419.1 # translation distance for the middle tie
		z2=952.5 #Translation distance for the side tie
		px1=-228.6 #half length of the tie
		px2=228.6
		py1=0.0
		py2=-177.8 #depth of the tie
		vtce=((px1,py1),(px2,py1),(px2,py2),(px1,py2))
		Tiesketches1=myModel.ConstrainedSketch(name='Tiesection1', sheetSize=800.0)
		Tiesketches1.Line(point1=vtce[0],point2=vtce[1])	
		Tiesketches1.Line(point1=vtce[1],point2=vtce[2])
		Tiesketches1.Line(point1=vtce[2],point2=vtce[3])	
		Tiesketches1.Line(point1=vtce[3],point2=vtce[0])
		Tiesketches2=myModel.ConstrainedSketch(name='Tiesection2', sheetSize=800.0)
		Tiesketches2.Line(point1=vtce[0],point2=vtce[1])	
		Tiesketches2.Line(point1=vtce[1],point2=vtce[2])
		Tiesketches2.Line(point1=vtce[2],point2=vtce[3])	
		Tiesketches2.Line(point1=vtce[3],point2=vtce[0])
		#Create Tie part by extrusion
		mysidetiePart1=myModel.Part(name='SideTie1_'+str(p),dimensionality=THREE_D,
		type=DEFORMABLE_BODY)
		mysidetiePart1.BaseSolidExtrude(depth=l1, sketch=Tiesketches1)
		mysidetiePart2=myModel.Part(name='SideTie2_'+str(p),dimensionality=THREE_D,
		type=DEFORMABLE_BODY)
		mysidetiePart2.BaseSolidExtrude(depth=l1,sketch=Tiesketches2)
		mymiddletiePart1=myModel.Part(name='MiddleTie1_'+str(p),dimensionality=THREE_D,
		type=DEFORMABLE_BODY)
		mymiddletiePart1.BaseSolidExtrude(depth=l2, sketch=Tiesketches1)
		### Create Materials
		myModel.Material(name='steel rail')
		myModel.materials['steel rail'].Elastic(table=((Erail, 0.3), 
			))
		myModel.Material(name='concrete ties')
		myModel.materials['concrete ties'].Elastic(table=((Etie, 
			0.4), ))
		myModel.Material(name='steel inclusion')
		myModel.materials['steel inclusion'].Elastic(table=((Eincl, 
			0.3), ))
		myModel.HomogeneousSolidSection(material='steel rail', name=
			'rail', thickness=None)
		myModel.HomogeneousSolidSection(material='steel inclusion', 
			name='inclusion', thickness=None)
		myModel.HomogeneousSolidSection(material='concrete ties', name=
			'concrete', thickness=None)
		myModel.parts['inclusion'].SectionAssignment(offset=0.0, 
			offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
			cells=myModel.parts['inclusion'].cells.findAt(((xb1,yci1,0.0), ), )),
				sectionName='inclusion', thicknessAssignment=FROM_SECTION)	
			
		RailPart.SectionAssignment(offset=0.0,
			offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
						cells=RailPart.cells.findAt(((0.0, 169.0, 800.1), ), ((1.0, 20.0, 
			801.0), ), ((-9.143037, 109.693032, 994.833374), ), ((9.516726, 
			114.698953, 249.797251), ), )), sectionName='rail', thicknessAssignment=FROM_SECTION)
		mysidetiePart1.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
						cells=mysidetiePart1.cells.findAt(((px1/2,(py2/2),l1/2), 
							), )), sectionName='concrete', thicknessAssignment=FROM_SECTION)	
		mysidetiePart2.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
						cells=mysidetiePart2.cells.findAt(((px1/2,(py2/2),l1/2), 
							), )), sectionName='concrete', thicknessAssignment=FROM_SECTION)
		mymiddletiePart1.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
						cells=mymiddletiePart1.cells.findAt(((px1/2,(py2/2),l1/2), 
							), )), sectionName='concrete', thicknessAssignment=FROM_SECTION)
		
		mysidetiePart1.DatumPlaneByPrincipalPlane(offset=
			-169.8, principalPlane=XZPLANE) #Datum[3] for sidetie1
			
		mysidetiePart1.DatumPlaneByPrincipalPlane(offset=
			-8.0, principalPlane=XZPLANE) #Datum[4] for sidetie1
			
		mysidetiePart2.DatumPlaneByPrincipalPlane(offset=
			-169.8, principalPlane=XZPLANE) #Datum[3] for sidetie2
			
		mysidetiePart2.DatumPlaneByPrincipalPlane(offset=
			-8.0, principalPlane=XZPLANE) #Datum[4] for sidetie2
			
		mysidetiePart1.PartitionCellByDatumPlane(cells=
			mysidetiePart1.cells.findAt(((5.0, -5.0, 5.0), ), ), datumPlane=
				mysidetiePart1.datums[3])
				
		mysidetiePart1.PartitionCellByDatumPlane(cells=
			mysidetiePart1.cells.findAt(((5.0, -5.0, 5.0), ), ), datumPlane=
				mysidetiePart1.datums[4])
				
		mysidetiePart2.PartitionCellByDatumPlane(cells=
			mysidetiePart2.cells.findAt(((5.0, -5.0, 5.0), ), ), datumPlane=
				mysidetiePart2.datums[3])
				
		mysidetiePart2.PartitionCellByDatumPlane(cells=
			mysidetiePart2.cells.findAt(((5.0, -5.0, 5.0), ), ), datumPlane=
				mysidetiePart2.datums[4])
			
		##create Assembly instances:
		myModel.rootAssembly.DatumCsysByDefault(CARTESIAN)
		RailInst=myModel.rootAssembly.Instance(dependent=OFF, name='RailInst', 
				part=RailPart)
		InclInstance=myModel.rootAssembly.Instance(dependent=OFF, name='InclusionInst', 
				part=partInclusion)	
		SideTie1Inst=myModel.rootAssembly.Instance(dependent=OFF, name='SideTie1Inst', 
				part=mysidetiePart1)					
		SideTie2Inst=myModel.rootAssembly.Instance(dependent=OFF, name='SideTie2Inst', 
				part=mysidetiePart2)
		MiddleTieInst=myModel.rootAssembly.Instance(dependent=OFF, name='MiddleTieInst', 
				part=mymiddletiePart1)
		#Create a coordinate system for the load
		myModel.rootAssembly.DatumCsysByThreePoints(coordSysType=
			CARTESIAN, name='Datum Coord for Load', origin=(xce, H, zce), point1=(
			10.0, H, zce), point2=(10.0, 180.0, zce))	
		#Translate one side tie and the middle tie along the rail
		myModel.rootAssembly.translate(instanceList=('SideTie2Inst',), 
						vector=(0.0, 0.0, z2))
		myModel.rootAssembly.translate(instanceList=('MiddleTieInst',), 
						vector=(0.0, 0.0, z1))					
		# Translate the inclusion
		myModel.rootAssembly.translate(instanceList=('InclusionInst',), 
						vector=(xci1,0.0,Trin ))				
		# mdb.models['SimTrack'].parts['RailLeft'].SectionAssignment(offset=0.0, 
			# offsetField='', offsetType=MIDDLE_SURFACE, region=
			# mdb.models['SimTrack'].parts['RailLeft'].sets['Set-2'], sectionName='rail', 
			# thicknessAssignment=FROM_SECTION)
		# mdb.models['SimTrack'].Part(compressFeatureList=ON, mirrorPlane=YZPLANE, name=
			# 'RailLeft-Copy', objectToCopy=mdb.models['SimTrack'].parts['RailLeft'])
		myModel.StaticStep(description='load applied', initialInc=0.1,maxInc=0.5, 
			name='Step-1', previous='Initial')	
		myModel.ExpressionField(description='elliptic contact load', expression='sqrt(1.0-Z*Z/pow('+str(a)+',2)-X*X/pow('+str(b)+',2))', 
			localCsys=myModel.rootAssembly.datums[12], name='Contact load')
		myModel.SurfaceTraction(createStepName='Step-1', 
			directionVector=((0.0, 0.0, 0.0), (0.0, -1.0, 0.0)), distributionType=FIELD,
			field='Contact load', localCsys=myModel.rootAssembly.datums[12], magnitude=PLoad, name='ContactLoad', region=
			Region(side1Faces=RailInst.faces.findAt(((-7.90097, 171.362215, 800.099996), ), )), resultant=ON, traction=GENERAL)
						
			#Create Datum planes
		##Apply the load
		#Surface area bound:
		#curve equation: (x-xc)^2+(y-yc)^2=Rh^2
		#Rail top curve: y=yc+sqrt(Rh^2-(x-xc)^2)
		yc=H-R3
		xc=xb1
		Rh=R3
		x1=px1/20
		y1=yc+sqrt(Rh**2-(x1-xc)**2)
			# Apply the longitudinal load
		myLoadSurface=myModel.rootAssembly.Surface(name='Surf-2', side1Faces=
			myModel.rootAssembly.instances['RailInst'].faces.findAt(
			((xb1,yci1,0.0), ),((xb1,yci1,RL), ) ))
		myLoad=myModel.Pressure(amplitude=UNSET, createStepName='Step-1', 
			distributionType=UNIFORM, field='', magnitude=Pressure2, name='Load-2', region=
			myLoadSurface)
		#Apply vertical load
		# myModel.StaticStep(description='load applied', initialInc=0.1, 
			# name='Step-2', previous='Step-1')	
		# myLoadSurface=myModel.rootAssembly.Surface(name='Surf-1', side1Faces=
			# myModel.rootAssembly.instances['RailInst'].faces.findAt(
			# ((x1,y1,RL-lin), ), ))
		# myLoad=myModel.Pressure(amplitude=UNSET, createStepName='Step-2', 
			# distributionType=UNIFORM, field='', magnitude=Pressure1, name='Load-1', region=
			# myLoadSurface)
			
		# Create dummy instance
				##create a dummy element
		sketchd=myModel.ConstrainedSketch(name='dummy section', sheetSize=400.0)
		sketchd.Line(point1=(0,0),point2=(rin,0))
		sketchd.Line(point1=(rin,0),point2=(rin,rin))
		sketchd.Line(point1=(rin,rin),point2=(0,rin))	
		sketchd.Line(point1=(0,rin),point2=(0,0))
		mydummyPart=myModel.Part(name='dummy',dimensionality=THREE_D,
			type=DEFORMABLE_BODY)
		mydummyPart.BaseSolidExtrude(depth=rin, sketch=sketchd)
		mydummyPart.SectionAssignment(offset=0.0, 
					offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
					cells=mydummyPart.cells.findAt(((0,0,rin), 
					), )), sectionName='rail', thicknessAssignment=
						FROM_SECTION)

		#Creating a dummy Instance
		myModel.rootAssembly.DatumCsysByDefault(CARTESIAN)
		DummyInst=myModel.rootAssembly.Instance(dependent=ON, name='dmyInst', part=mydummyPart)
		#Translate the the dummy instance at a distance Lr
		myModel.rootAssembly.translate(instanceList=('dmyInst', ), 
				vector=(0.0, 0.0, RL))	
		#Applying Encastred BC
		myModel.EncastreBC(createStepName='Initial', name='BC-dummy', region=Region(faces=DummyInst.faces.findAt(((rin/2, rin/2, RL+rin),),)))
		#Mesh the dummy part
		mydummyPart.setElementType(elemTypes=(ElemType(elemCode=C3D20R, elemLibrary=STANDARD), ElemType(elemCode=C3D15, 
			elemLibrary=STANDARD), ElemType(elemCode=C3D10, elemLibrary=STANDARD)), 
			regions=(mydummyPart.cells.findAt(((rin/2, rin/2, rin/2), ), ), ))
		mydummyPart.seedPart(deviationFactor=0.2, minSizeFactor=0.2, size=rin)
		mydummyPart.generateMesh()
		#create dummy node
		myModel.rootAssembly.Set(name='Set-3', vertices=DummyInst.vertices.findAt(((0.0, 0.0, RL), ), ))
		##CREATING TIE PARTS	
		#Apply boundary conditions

			## Creating tie contraints
		#1. Create tie surfaces
		myRailBoSurf1=myModel.rootAssembly.Surface(name='RailBottomSurf1', side1Faces=
			myModel.rootAssembly.instances['RailInst'].faces.findAt(((49.741669, 0.0, 249.766663), ),) )
		myRailBoSurf2=myModel.rootAssembly.Surface(name='RailBottomSurf2', side1Faces=
			myModel.rootAssembly.instances['RailInst'].faces.findAt(((-49.741669, 0.0, 994.833374), ),) )	
		mysideT1Surf=myModel.rootAssembly.Surface(name='sideT1Surf', side1Faces=
			SideTie1Inst.faces.findAt(
			((xb1,yb1,l1/4), ), ))
		mysideT2Surf=myModel.rootAssembly.Surface(name='sideT2Surf', side1Faces=
			SideTie2Inst.faces.findAt(
			((xb1,yb1,RL-l1/2), ), ))
		mymiT1Surf=myModel.rootAssembly.Surface(name='miT1Surf', side1Faces=
			MiddleTieInst.faces.findAt(
			((xb1,yb1,RL/2), ), ))
		# Tie constraints
		myModel.Tie(adjust=ON, constraintEnforcement=NODE_TO_SURFACE, 
			master=myRailBoSurf1, 
			name='Constraint-1', positionToleranceMethod=COMPUTED, slave=
			mysideT1Surf, thickness=ON, 
			tieRotations=ON)
		myModel.Tie(adjust=ON, constraintEnforcement=NODE_TO_SURFACE, 
			master=myRailBoSurf2, 
			name='Constraint-2', positionToleranceMethod=COMPUTED, slave=
			mysideT2Surf, thickness=ON, 
			tieRotations=ON)
		myModel.Tie(adjust=ON, constraintEnforcement=NODE_TO_SURFACE, 
			master=myRailBoSurf1, 
			name='Constraint-3', positionToleranceMethod=COMPUTED, slave=
			mymiT1Surf, thickness=ON, 
			tieRotations=ON)

		#Surfaces at inclusion interface and Tie constraints
		##inclusion surface points
		xi2=xci1	
		yi2=yci1+sqrt(rin**2-(xi2-xci1)**2)
		zi2=zci1+sqrt(rin**2-(yi2-yci1)**2)
		zi2d=zci1-sqrt(rin**2-(yi2-yci1)**2)
		# ##create another datum plane
		# mdb.models['SimTrack'].parts['RailRight'].DatumPlaneByPrincipalPlane(offset=
		 # xci1, principalPlane=YZPLANE) #Datum[9]  

		#Partition rail inclusion surrounding surface
		# RailPart.PartitionFaceByDatumPlane(datumPlane=
			# RailPart.datums[6], faces=
				# RailPart.faces.findAt(((xi2,yi2,zci1), ), ))
		##surfaces 
		myIncSurf=myModel.rootAssembly.Surface(name='inclSurf', side1Faces=
			InclInstance.faces.findAt(((xi2,yci1,Trin+rin), ), ))
				
		myRailIncSurf=myModel.rootAssembly.Surface(name='RailInclSurf', side1Faces=
			RailInst.faces.findAt(((xi2,yi2,zi2), ),((xi2,yi2,zi2d), ), ))
		##Tie constraints
		myModel.Tie(adjust=ON, constraintEnforcement=NODE_TO_SURFACE, 
			master=myRailIncSurf, 
				name='Constraint-4', positionToleranceMethod=COMPUTED, slave=
					myIncSurf, thickness=ON, 
						tieRotations=ON)
		myRailBoSurf1=myModel.rootAssembly.Surface(name='RailBottomSurf1', side1Faces=
			myModel.rootAssembly.instances['RailInst'].faces.findAt(((xb1,yb1,l1), ),) )
		myRailBoSurf2=myModel.rootAssembly.Surface(name='RailBottomSurf2', side1Faces=
			myModel.rootAssembly.instances['RailInst'].faces.findAt(((xb1,yb1,RL-l1/2), ),) )	
			#Start Meshing
		#Apply encastred boundary conditions
		myModel.EncastreBC(createStepName='Initial', localCsys=None, 
			name='BC-1', region=Region(
			faces=myModel.rootAssembly.instances['SideTie1Inst'].faces.findAt(
			((76.200002, -177.8, 76.200002), ), )+\
			myModel.rootAssembly.instances['MiddleTieInst'].faces.findAt(
			((76.200002, -177.8, 571.500004), ), )+\
			myModel.rootAssembly.instances['SideTie2Inst'].faces.findAt(
			((76.200002, -177.8, 1028.700002), ), )))
		#Redefine surfaces
		myModel.rootAssembly.Surface(name='RailBottomSurf1', side1Faces=
			myModel.rootAssembly.instances['RailInst'].faces.findAt(((
			49.741669, 0.0, 249.766663), )))
		x1=xci1-0.666583
		z1=799.711771
		x2=xci1+0.098876
		z2=799.711771
		y1=yci1+sqrt(rin**2-(x1-xci1)**2-(z1-zci1)**2)
		y2=yci1+sqrt(rin**2-(x2-xci1)**2-(z2-zci1)**2)
		#     -0.666583, 159.843553, 799.711771), ), ((0.098876, 159.684489, 799.711771), 
			#), ))
		myModel.rootAssembly.Surface(name='RailInclSurf', side1Faces=
			myModel.rootAssembly.instances['RailInst'].faces.findAt(((
			x1, y1, z1), ), ((x2, y2, z2), ), ))
		myModel.rootAssembly.Surface(name='Surf-1', side1Faces=
			myModel.rootAssembly.instances['RailInst'].faces.findAt(((
			-3.292674, 171.434755, 817.033346), )))	

		# myModel.parts['RailRight'].setMeshControls(elemShape=HEX, 
			# regions= myModel.parts['RailRight'].cells.findAt(((xb1,yb1,l1), ),((xb1,yb1,Trin),),
			# ((xb1,yb1,RL-l1/2), ), ), technique=FREE)
		#
		#Mesh the instances 
		myModel.rootAssembly.setElementType(elemTypes=(ElemType(
			elemCode=C3D8R, elemLibrary=STANDARD, secondOrderAccuracy=OFF, 
				kinematicSplit=AVERAGE_STRAIN, hourglassControl=DEFAULT, 
					distortionControl=DEFAULT), ElemType(elemCode=C3D6, elemLibrary=STANDARD), 
						ElemType(elemCode=C3D4, elemLibrary=STANDARD)), regions=(
							RailInst.cells.findAt(((xb1,yb1,l1), ),((xb1,yb1,Trin),),((xb1,yb1+0.5,RL-l1), ), )+\
							SideTie2Inst.cells.findAt(((xb1,yb1,RL-l1/2),),)+\
							MiddleTieInst.cells.findAt(((xb1,yb1,RL/2), ), )+\
							SideTie1Inst.cells.findAt(((xb1,yb1,l1), ), ),))

		myModel.rootAssembly.setMeshControls(elemShape=TET, 
			regions= RailInst.cells.findAt(((xb1,H-0.5,Trin), ), ), technique=FREE)
		#Mesh the inclusion
		myModel.rootAssembly.setMeshControls(elemShape=TET, 
			regions= InclInstance.cells.findAt(((xci1,yci1,Trin), ), ), technique=FREE)	
		myModel.rootAssembly.setElementType(elemTypes=(ElemType(
			elemCode=C3D20R, elemLibrary=STANDARD), ElemType(elemCode=C3D15, 
				elemLibrary=STANDARD), ElemType(elemCode=C3D10, elemLibrary=STANDARD)), 
					regions=(RailInst.cells.findAt(((xb1,H-1.0,Trin), ), )+\
						InclInstance.cells.findAt(((xci1,yci1,Trin), ),),))	
		myModel.rootAssembly.seedPartInstance(deviationFactor=devf, 
			minSizeFactor=0.0125,regions=(RailInst,SideTie1Inst,MiddleTieInst,SideTie2Inst), size=Elsize)
		myModel.rootAssembly.seedPartInstance(deviationFactor=0.5, 
			minSizeFactor=0.5,regions=(InclInstance, ), size=0.1)	
		myModel.rootAssembly.seedEdgeBySize(constraint=FINER, 
			deviationFactor=0.5, edges=
				RailInst.edges.findAt(
					((xci1+rin,yci1,Trin), ),((xci1-rin,yci1,Trin), ),
					((xci1,yci1,Trin+rin), ),((xci1,yci1,Trin-rin), ), ),
					minSizeFactor=0.5, size=0.1)
		myModel.rootAssembly.seedEdgeBySize(constraint=FINER, 
				deviationFactor=0.5, edges=RailInst.edges.findAt(((-14.802219, 171.141788, 790.9), )),
					minSizeFactor=0.99, size=0.5)			
		myModel.rootAssembly.generateMesh(regions=RailInst.cells.findAt(((xci1,H-0.4,Trin), ), ))			
		myModel.rootAssembly.generateMesh(regions=
			RailInst.cells.findAt(((xb1,H-1.0,Trin), ),((xb1,yb1,l1), ),((xb1,yb1+0.5,Trin),),
			((xb1,yb1,RL-l1/2), ), ),seedConstraintOverride=ON)
		myModel.rootAssembly.generateMesh(regions=(
		   SideTie2Inst,MiddleTieInst,SideTie1Inst,InclInstance ))
		#Create sets for elements around inclusion	
		myModel.rootAssembly.Set(faces= RailInst.faces.findAt(((
				x1, y1, z1), ), ((x2, y2, z2), ), ), name='ElemSet-1')   
		#Creates node sets for periodic boundary conditions
		myModel.rootAssembly.Set(faces= SideTie1Inst.faces.findAt(((
		5.0, -50.0, 0.000), )), name='Set-1')
		myModel.rootAssembly.Set(faces= SideTie2Inst.faces.findAt(((
		5.0, -50.0, 1066.8), )), name='Set-2') 
		
		#Apply Periodic boundary conditions

		set1=myModel.rootAssembly.sets['Set-1'].nodes
		set2=myModel.rootAssembly.sets['Set-2'].nodes
		set3=myModel.rootAssembly.sets['Set-3'].nodes
		# node1=mdb.models['Model-1'].rootAssembly.instances['fullrail-1'].nodes
		node2=DummyInst.nodes
		set2new=[]
		set1new=[]
		x2s=[]
		y2s=[]
		z2s=[]
		z1s=[]
		# cst=152.4
		#compute the minimum distance in xy plane between one node from 
		#one face and each node from the other face and create a list of matched nodes
		x=[]
		y=[]
		z=[]
		def ArrNodeSet(set2,set2new):
			for n in range(len(set2)):
				x1=set2[n].coordinates[0]
				y1=set2[n].coordinates[1]
				z1=set2[1].coordinates[2]
				labe=set2[n].label
				x2s.append(x)
				y2s.append(y)
				z2s.append(z)
				set2new.append([labe,x1, y1, z1])
			z2cst=max(z2s)
		ArrNodeSet(set2,set2new)
		ArrNodeSet(set1,set1new)
		Tm=1 #Condition to match coordinates in set 2 to coordinates in set 1 else vice versa.
		NodeEqxy=[]
		set2n=[]
		def editingNodes(set1new,set2new,Tm):
			# set2n=[]
			for i in range(len(set1new)):
				d12xyp=[((set2new[n][1]-set1new[i][1])**2+(set2new[n][2]-set1new[i][2])**2)**0.5 for n in range(len(set2new))]
				dnmin=min(d12xyp) 
				mlab=d12xyp.index(dnmin) # Index of matching node in set2 to a node of index i in set1
				if (Tm==1):
					nameOfset='Set-2'
					#lab1 and lab2 are labels of matching nodes.
					lab1=set1new[i][0]
					lab2=set2new[mlab][0]
					no2nx=set1new[i][1]
					no2ny=set1new[i][2]
					no2nz=set1new[1][3]
					set2n.append([no2nx,no2ny,no2nz,i,mlab])
					NodeEqxy.append([set1new[i],lab1,lab2])
					myModel.rootAssembly.editNode(coordinate1=set2n[i][0], 
					coordinate2=set2n[i][1], coordinate3=set2n[i][2],
					nodes=myModel.rootAssembly.sets[nameOfset].nodes[set2n[i][4]])	
				elif (Tm==2):
					nameOfset='Set-1'
					labn=set1new[nj][0]
					no2nx=set2new[i][1]
					no2ny=set2new[i][2]
					no2nz=set2new[1][3]
					set2n.append([no2nx,no2ny,no2nz,mlab,lab2])
					NodeEqxy.append([set1new[i],lab1,lab2])
		#start editing nodes with new coordinates keeping the node number the same
					myModel.rootAssembly.editNode(coordinate1=set2n[i][0], 
					coordinate2=set2n[i][1], coordinate3=set2n[i][2],
					nodes=myModel.rootAssembly.sets[nameOfset].nodes[set2n[i][3]])		
		##find extra nodes on the other face that has more nodes
		#set2new-set2n
		editingNodes(set1new,set2new,Tm)
		if (range(len(set2))>range(len(set1))):
			set2node=[set2new[n][0] for n in range(len(set2new))]
			set2noden=[set2n[n][0] for n in range(len(set2n))]
			hangnode=list(set(set2node)-set(set2noden))
			#Edit these nodes 
			SetToEdit=[]
			def indexOfnodes(hangnode,set2,SetToEdit):
		#This function returns the set of nodes with their coordinates and their index in the nodes set)
				for k in range(len(hangnode)):
					NoToEdit=hangnode[k]
					for l in range(len(set2)):
						Nodenum=set2[l].label
						if (Nodenum==NoToEdit):
							xEdit=set2[l].coordinates[0]
							yEdit=set2[l].coordinates[1]
							zEdit=zcst
							SetToEdit.append([Nodenum,xEdit,yEdit,ZEdit,l])
							break
						elif (Nodenum!=NoToEdit):
							continue
			Tmm=2
			editingNodes(SetToEdit,Set1new,Tmm,'Set-2')
			#Print nodes numbers pairs of the same coordinates
			print 'Equal x y coordinates node pairs are', NodeEqxy
			print 'hanging nodes are', SetToEdit
		#	Start defining boundary conditions
		elif (range(len(set2))==range(len(set1))):
			for i in range(len(set1)):
		#create the names of node sets
				setname1='set1a-'+str(i)
				setname2='set2a-'+str(i)
				setname3='set3'
		##creates the node sets
				# node1=mdb.models['Model-1'].rootAssembly.sets['Set-1'].nodes
				# node2=mdb.models['Model-1'].rootAssembly.sets['Set-2'].nodes
				# node3=mdb.models['Model-1'].rootAssembly.sets['Set-3'].nodes
				st1=set1[set2n[i][3]:set2n[i][3]+1];
				st2=set2[set2n[i][4]:set2n[i][4]+1];
				# nodelabel1=nst1.label
				# nodelabel2=nst2.label
				# st1=set1[nodelabel1-1:nodelabel1]
				# st2=set1[nodelabel2-1:nodelabel2]
				#use global index sets to create node sets
		#Create a dummy node
				st3=set3;
		#Create the node sets
				myModel.rootAssembly.Set(name=setname1, nodes=st1)
				myModel.rootAssembly.Set(name=setname2, nodes=st2)	
				myModel.rootAssembly.Set(name=setname3, nodes=st3)	
		#Apply constraints equations for periodic boundary conditions	
				myModel.Equation(name='ConstraintX1-'+str(i), terms=((1.0, 
					setname1, 1), (-1.0, setname2, 1)))
				myModel.Equation(name='ConstraintY1-'+str(i), terms=((1.0, 
					setname1, 2), (-1.0, setname2, 2)))
				myModel.Equation(name='ConstraintZ1-'+str(i), terms=((1.0, 
					setname1, 3), (-1.0, setname2, 3),(1, setname3, 3)))	
					
			 #Set Field output requests
		myModel.fieldOutputRequests['F-Output-1'].setValues(variables=(
			'S', 'MISES', 'MISESMAX', 'TRIAX', 'MISESONLY', 'E', 'PE', 'PEEQ', 'PEMAG', 
			'EE', 'LE', 'U', 'RF', 'CF', 'CSTRESS', 'CDISP'))
		##Create Job output
		mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
				explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
				memory=90, memoryUnits=PERCENTAGE, model='SimTrack'+str(p), modelPrint=OFF, 
				multiprocessingMode=DEFAULT, name='SimTrack1_k_'+str(p), nodalOutputPrecision=SINGLE, 
				numCpus=1, queue=None, scratch='', type=ANALYSIS, userSubroutine='',waitHours=0, waitMinutes=0)	
		mdb.jobs['SimTrack1_k_'+str(p)].writeInput(consistencyChecking=OFF)
		g_incr.write(str(p)+'\n')
g_incr.close()
		