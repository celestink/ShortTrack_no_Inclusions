#These scripts are for analysing output database stresses in the rail with inclusion
#This script is an update version of the script PostProcess_v4nosoft.py
#It apply to models with inclusion that have the uniform materials as the rail therefore acting as if there was no inclusion
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
import random	
import odbAccess
from  abaqus import session
import numpy as np
# nt=10 # number of ties
# nl=2*nt-1 #number of load positions
m=[[0,12],[12,12],[24,3],[36,3],[48,3],[60,3],[72,3],[84,3],[96,12],[108,3],[120,3],[132,3],[144,3],[156,3],[168,3]]
#m=[[j+12,k] for j in range(15) for k in [12,12,3,3,3,3,3,3,12,3,3,3,3,3,3]
#nm=12
sigmf=936
b=-0.089
for p in range(15):
	MPrs=[]
	MDirs=[]
	EleMaxs=[]
	MaxPrinc=[]	
	EleMaxs=[]
	g=open('RailMaxPrinc_so_'+str(p)+'.txt','w')
	f=open('vonMiseMax_so_'+str(p)+'.txt','w')
	Samp=open('Ampl_Stress_so_'+str(p)+'.txt','w')
	Se=open('Mean_Princ_Stress_so_'+str(p)+'.txt','w')
	Seff=open('Eff_Stress_so_'+str(p)+'.txt','w')
	Meff=open('Mean_Stress_so_'+str(p)+'.txt','w')
	Nc=open('N_cycles_so_'+str(p)+'.txt','w')
	vm=open('vonMises_so_'+str(p)+'.txt','w')
	#!/usr/bin/python
	#call the files that has the elements of maximum principal stresses as well as von Mises stresses
	fh = open('C:\Users\celestink\Documents\TestingStation\EllptContact2\AllSoftInclusions\RailMaxPrinc_'+str(p)+'.txt');
	von = open('C:\Users\celestink\Documents\TestingStation\EllptContact2\AllSoftInclusions\MiseMax_'+str(p)+'.txt')
	x = []
	xv=[]
	for line in fh.readlines():
		y = [value for value in line.split()]
		x.append(int(y[4]))
	fh.close()
	print x
	for line in von.readlines():
		yv = [value for value in line.split()]
		xv.append(int(yv[1]))
	von.close()
	print xv
	for k in range(m[p][0],m[p][0]+m[p][1]):
	# for k in [0,1]:
	#Create a file that contain maximum stresses for every every standard deviation given to random stiffness in the ties
		# odb1=openOdb('SimTrack1_k_'+str(k)+'.odb')
		l=k-m[p][0]
		newSetElements = (x[l],)
		newSetElementsv = (xv[l],)
		fileName = 'SimTrack1_k_'+str(k)+'.odb' # name of odb file
		odb1= odbAccess.openOdb(path = fileName, readOnly = False) #open odb 
		
		#ElementSetFromElementLabels constructor to build new set (NewElSetso) for max principal stress element 
		#that was found with soft inclusion into instance (RAILINST)
		
		odb1.rootAssembly.ElementSetFromElementLabels(name = 'Newelsetsoft', elementLabels =
		(('RAILINST',newSetElements),))
		
		# ElementSetFromElementLabels constructor to build new set (NewElSetsovon) for max von Mises stress element 
		#that was found with soft inclusion into instance (RAILINST)		
		odb1.rootAssembly.ElementSetFromElementLabels(name = 'NewelsetsovonM', elementLabels =
		(('RAILINST',newSetElementsv),))
		#
		#Maximum stress:
		InclElem=odb1.rootAssembly.elementSets['Newelsetsoft']
		stresses=odb1.steps['Step-1'].frames[-1].fieldOutputs['S']
		InclStress=stresses.getSubset(region=InclElem)
		#Maximum Mises:
		InclElemvon=odb1.rootAssembly.elementSets['NewelsetsovonM']
		stresses=odb1.steps['Step-1'].frames[-1].fieldOutputs['S']
		InclStressv=stresses.getSubset(region=InclElemvon)
		#gp=open('RailMaxPrinc-'+str(k)+'.txt','w')
		# maxs33=[]
		# maxs33=[]
		# # for j in range(l0,lf1):
		# # session.upgradeOdb("G:/EDissertProject/odbs04-19-14/JobAuto-"+str(j)+".odb", 
			# # "G:/EDissertProject/odbfiles-1/JobAuto-"+str(j)+"-1.odb",)
		# #Open output database model
		# odb1=openOdb('SimTrack1_k_'+str(k)+'.odb')
		# #Stresses:
		# stresses=odb1.steps['Step-1'].frames[-1].fieldOutputs['S']
		# InclElem=odb1.rootAssembly.elementSets['NewElSet']
		# InclStress=stresses.getSubset(region=InclElem)
		# f=open('RailbendingStress-'+str(k)+'.txt','w')
		dir=open('Princ-Direction-'+str(k)+'.txt','w')
		#
		# s33=[]
		MaxPrinc=[]
		ElemNo=[]
		ElemNov=[]
		DirVec=[]
		vonMises=[]
		Princip0=[]
		Princip1=[]
		Princip2=[]
		for i in range(0,len(InclStress.values)):
			str0=InclStress.values[i].data[0] 	#SGMxx
			str1=InclStress.values[i].data[1]	#SGMyy
			str2=InclStress.values[i].data[2]	#SGMzz
			str3=InclStress.values[i].data[3]	#SGMxy
			str4=InclStress.values[i].data[4]	#SGMxz
			str5=InclStress.values[i].data[5]	#SGMyz
			Mise=InclStressv.values[i].mises
			vonMises.append(Mise)
			ElLab=InclStress.values[i].elementLabel
			ElLabv=InclStressv.values[i].elementLabel
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
			ElemNov.append(ElLabv)
			# f.write(str(astress)+'\n')
			dir.write(str(MaxPrinc[i])+'\t')
			[dir.write(str(DirVec[i][k])+'\t') for k in range(3)]
			dir.write('\n')
			
			#List of principal stresses
			Princip0.append(evals[0])
			Princip1.append(evals[1])
			Princip2.append(evals[2])			
		# Ms33=max(s33)#Find maximum bending stress
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
		
		vonmax=max(vonMises)#Find maximum VonMises stress
		Ind_vonMax=vonMises.index(vonmax)
		ElvonMax=ElemNo[Ind_vonMax]
		#Computing effective stresses and fatigue cycles to crack initiation
		effstr0max=max(Princip0)
		if effstr0max<0.0:
			effstr0max=0.0
			effstr0min=min(Princip0)
		else:
			effstr0min=0.0
		sigm0a=(effstr0max-effstr0min)/2.0
		sigm0m=(effstr0max+effstr0min)/2.0
		effstr1max=max(Princip1)
		if effstr1max<0.0:
			effstr1max=0.0
			effstr1min=min(Princip1)
		else:
			effstr1min=0.0
		sigm1a=(effstr1max-effstr1min)/2.0
		sigm1m=(effstr1max+effstr1min)/2.0
		effstr2max=max(Princip2)
		if effstr2max<0.0:
			effstr2max=0.0
			effstr2min=min(Princip2)
		else:
			effstr2min=0.0
		sigm2a=(effstr2max-effstr2min)/2.0
		sigm2m=(effstr2max+effstr2min)/2.0
		sigmae=1.0/sqrt(2.0)*(((sigm0a-sigm1a)**2+(sigm1a-sigm2a)**2+(sigm2a-sigm0a)**2)**0.5)
		sigmam=1.0/sqrt(2.0)*(((sigm0m-sigm1m)**2+(sigm1m-sigm2m)**2+(sigm2m-sigm0m)**2)**0.5)
		Nfun=(1.0/2.0*(sigmae/sigmf)**(1.0/b))*((1.0-sigmam/sigmf)**(-1.0/b))
		vm.write(str(vonmax)+'\n')
		Samp.write(str(sigm0a)+'\t')
		Samp.write(str(sigm1a)+'\t')
		Samp.write(str(sigm2a)+'\n')
		Se.write(str(sigm0m)+'\t')
		Se.write(str(sigm1m)+'\t')
		Se.write(str(sigm2m)+'\n')
		Seff.write(str(sigmae)+'\n')
		Meff.write(str(sigmam)+'\n')
		Nc.write(str(Nfun)+'\n')
		# g.write(str(Ms33)+'\n')
		#gp.write(str(MPr)+'\t')
		#[gp.write(str(MDir[k])+'\t') for k in range(3)]
		#gp.write('\n')
		g.write(str(MPr)+'\t')
		[g.write(str(MDir[k])+'\t') for k in range(3)]
		g.write(str(ElMax)+'\t')
		g.write('\n')
		f.write(str(vonmax)+'\t')
		f.write(str(ElvonMax)+'\t')
		f.write('\n')
		# f.close()
		dir.close()
		# g.close()
		#gp.close()
	g.close()
	f.close()
	vm.close()
	Seff.close()
	Meff.close()
	Nc.close()
	Samp.close()
	Se.close()	