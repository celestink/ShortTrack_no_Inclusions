#These scripts are for analysing output database stresses in the rail
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random	
from odbAccess import *	
from  abaqus import session
import numpy as np
# nt=10 # number of ties
# nl=2*nt-1 #number of load positions
MPrs=[]
MDirs=[]
EleMaxs=[]
m=[[0,12],[12,12],[24,3],[36,3],[48,3],[60,3],[72,3],[84,3],[96,12],[108,3],[120,3],[132,3],[144,3],[156,3],[168,3]]
#m=[[j+12,k] for j in range(15) for k in [12,12,3,3,3,3,3,3,12,3,3,3,3,3,3]
#nm=12
for p in range(15):
	g=open('RailMaxPrinc_'+str(p)+'.txt','w')
	for k in range(m[p][0],m[p][0]+m[p][1]):
	# for k in [0,1]:
	#Create a file that contain maximum stresses for every every standard deviation given to random stiffness in the ties
		gp=open('RailMaxPrinc-'+str(k)+'.txt','w')
		maxs33=[]
		maxs33=[]
		# for j in range(l0,lf1):
		# session.upgradeOdb("G:/EDissertProject/odbs04-19-14/JobAuto-"+str(j)+".odb", 
			# "G:/EDissertProject/odbfiles-1/JobAuto-"+str(j)+"-1.odb",)
		#Open output database model
		odb1=openOdb('SimTrack1_k_'+str(k)+'.odb')
		#Stresses:
		stresses=odb1.steps['Step-1'].frames[-1].fieldOutputs['S']
		InclElem=odb1.rootAssembly.elementSets['ELEMSET-1']
		InclStress=stresses.getSubset(region=InclElem)
		# f=open('RailbendingStress-'+str(k)+'.txt','w')
		dir=open('Princ-Direction-'+str(k)+'.txt','w')
		#
		# s33=[]
		MaxPrinc=[]
		ElemNo=[]
		DirVec=[]
		for i in range(0,len(InclStress.values)):
			str0=InclStress.values[i].data[0] 	#SGMxx
			str1=InclStress.values[i].data[1]	#SGMyy
			str2=InclStress.values[i].data[2]	#SGMzz
			str3=InclStress.values[i].data[3]	#SGMxy
			str4=InclStress.values[i].data[4]	#SGMxz
			str5=InclStress.values[i].data[5]	#SGMyz	
			ElLab=InclStress.values[i].elementLabel
			ST=[[str0,str3,str4],[str3,str1,str5],[str4,str5,str2]] #Stress tensor at integration point in an element
			#Compute Eigenvalues of the stress tensor (principal stresses) and eigenvectors of the stress tensors (directions)
			# STneg=[[x1*(-1) for x1 in X] for x in ST]
			evals,evecs=np.linalg.eig(ST)
			maxpr_pos=max(evals)
			# evals_neg=[X*(-1) for X in evals]
			# maxpr_neg=max(evals_neg)
			# if maxpr_pos>=maxpr_neg:
			maxpr=maxpr_pos
			indp=list(evals).index(maxpr)
			# else:
				# maxpr=-maxpr_neg
				# indp=list(evals_neg).index(maxpr_neg)	
			# astress=abs(InclStress.values[i].data[2])
			DirVec.append([evecs[0][indp],evecs[1][indp],evecs[2][indp]])
			# PrincStr[i]=dict({"Max.princ.":maxpr,"Direction":DirVec})
			# s33.append(astress)
			MaxPrinc.append(maxpr)
			ElemNo.append(ElLab)
			# f.write(str(astress)+'\n')
			dir.write(str(MaxPrinc[i])+'\t')
			[dir.write(str(DirVec[i][k])+'\t') for k in range(3)]
			dir.write('\n')
			#Find maximum bending stress
		# Ms33=max(s33)
		MPr_p=max(MaxPrinc) #highest maximum principle in the rail inclusion interface nodesets
		MaxPrinc_n=[X*(-1) for X in MaxPrinc]
		# MPr_n=max(MaxPrinc_n)
		# if MPr_p>=MPr_n:
		MPr=MPr_p
		Ind_MPr=MaxPrinc.index(MPr)
		# else:
				# MPr=-MPr_n
				# Ind_MPr=MaxPrinc_n.index(MPr_n)
		MDir=DirVec[Ind_MPr]
		ElMax=ElemNo[Ind_MPr]
		#Append maximum bending stress on the rail in a different files
		# maxs33.append(s33)
		MPrs.append(MPr)
		MDirs.append(MDir)
		EleMaxs.append(ElMax)
		# g.write(str(Ms33)+'\n')
		gp.write(str(MPr)+'\t')
		[gp.write(str(MDir[k])+'\t') for k in range(3)]
		gp.write('\n')
		g.write(str(MPr)+'\t')
		[g.write(str(MDir[k])+'\t') for k in range(3)]
		g.write(str(ElMax)+'\t')
		g.write('\n')
		# f.close()
		dir.close()
		# g.close()
		gp.close()
	g.close()