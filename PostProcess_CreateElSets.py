#These scripts are for analysing output database stresses in the rail
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random	
from odbAccess import *	
from  abaqus import session
import numpy as np
import odbAccess
# nt=10 # number of ties
# nl=2*nt-1 #number of load positions
for k in range(12):
# for k in [0,1]:
	maxs33=[]
	# odb1=openOdb('SimTrack1_k_'+str(k)+'.odb')
	Elsets=[36693,	34895,	35594,	35393,	34883,	35085,	33391,	35321,	35356,	31207,	35105,	36887]
	newSetElements = (Elsets[k],)
	fileName = 'SimTrack1_k_'+str(k)+'.odb' # name of odb file
	odb1= odbAccess.openOdb(path = fileName, readOnly = False) #open odb 
	# ElementSetFromElementLabels constructor to build new set (NewSet) into
	#instance (PART-1-1)
	odb1.rootAssembly.ElementSetFromElementLabels(name = 'NewElSet', elementLabels =
	(('RAILINST',newSetElements),))
	odb1.close() #close odb file 