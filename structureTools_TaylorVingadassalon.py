#!/usr/bin/env python
#-*- coding : utf3 -*-
import sys
import math

import matplotlib.pyplot as plt
import numpy as np


## Retourne le dictionnaire atome-masse moleculaire
# @a : chemin du fichier de donnees atome-masse moleculaire
def lireAtoms(a):

	atome = dict()
	##########################
	# Chargement dico atomes #
	with open(a,"r") as fichier:
		fic=fichier.readlines()
		for l in fic:
			line = l.split()
			atome[line[0]]=(float)(line[1])
	return atome
	

## Retourne le dictionnaire correspondant au fichier PDB lu
# @a : chemin du fichier PDB
# @b : dictionnaire atomes-masse moleculaire
def lirePDB(a,b):

	
	#### Dictionnaires ####
	residu = dict()
	chaine = dict()
	atome = dict()
	atmWeight = lireAtoms(b)
	#######################
	conf=''

	with open(a, "r") as fichier:
		fic = fichier.readlines()
		for l in fic:
			
			if (conf=='' and l[0:6].strip()=='ATOM'):
				conf=l[16]
				print("IN:"+'\''+conf+'\''+"\t"+'\''+l[0:6].strip()+'\'')
			else:
				if(conf!=' '): print("conf1")
				if(l[0:6].strip() == "ATOM"): print("conf2:"+'\''+l[0:6].strip()+'\'')
				print('\''+conf+'\''+"\t"+'\''+l[0:6].strip()+'\'')
				
			if (l[0:4]=="ATOM" and l[16]==conf): #Scan des chaines d'interet
				##### Recuperation des donnees ###
				res=l[17:20].strip()
				nomAtome = l[13:16].strip()
				chName=l[21].strip()
				
				info={'ID': (l[7:11].strip()),
						'x'	: (float)(l[31:38].strip()),
						'y' : (float)(l[39:46].strip()),
						'z' : (float)(l[47:54].strip()),
						'atmW': atmWeight[l[77].strip()]
					}
				######## Fin recuperation ########
				#~ print(res)
				 # Chaine non repertoriee donc on l'ajoute
				if (chName not in chaine.keys()):
					chaine[chName] = {}
				
				## Residu non encore repertorie
				if(res not in chaine[chName].keys()):
					chaine[chName][res]={}
				
				chaine[chName][res][nomAtome]=info

	return chaine


## Ajoute pour chaque residu du dictionnaire b la position (x,y,z) de son centre de masse
# @a: dictionnaire issu de la commande lirePDB
def ajouterCentreDeMasse(a):
	
	for i in a.keys():			# On parcours la chaine
		for j in a[i].keys():	# On parcours les residus
			
			cdm = {'x':0,'y':0,'z':0}
			masseTotale=0
			for k in a[i][j].keys(): # Extraction des infos pour chaque atome
				atmW = a[i][j][k]['atmW']
				cdm['x']+=atmW*a[i][j][k]['x']
				cdm['y']+=atmW*a[i][j][k]['y']
				cdm['z']+=atmW*a[i][j][k]['z']
				
				masseTotale+=atmW
				
			cdm['x']/=masseTotale
			cdm['y']/=masseTotale
			cdm['z']/=masseTotale
			
			a[i][j]['cdm']=cdm

## prend en entrée deux atomes (donc les dico d'info des deux atomes) (x,y,z) et retourne la distance entre eux
def distanceAtomes(a,b):
	x1 = a['x']
	y1 = a['y']
	z1 = a['z']

	x2 = b['x']
	y2 = b['y']
	z2 = b['z']
	
	distance= math.sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2))
	return distance

## Calcule les distances residu-residu
# @a: dictionnaire issu de la fonction ajouterCentreDeMasse(a) 
def distance(a):
	mat = dict()
	
	for i in a.keys(): #Chaines i
		for j in a[i].keys(): # Resisud	j
			
			if str(j) not in mat.keys():
				mat[str(j)]={}
			
			#~ for k in a.keys(): # Chaine k
			for l in a[i].keys(): #residu l
			
				if j!=l:
					if str(l) not in mat.keys() or str(j) not in mat[str(l)].keys():
						
						distance=distanceAtomes(a[i][j]['cdm'],a[i][l]['cdm'])
						mat[str(j)][str(l)]={'val':distance}
							
	return mat


#~ ## Test si deux residus sont en interaction
#~ def draw(mat):
	#~ x,y = np.grid[]
	

## Affiche proprement les distances residu-redisu
# @a: dictionnaire issu de la fonction distance(a)
def printDistance(a):
	for i in a.keys():
		for j in a[i].keys():
			print("["+i+"]["+j+"] = "+str(a[i][j]['val']))

def repeat(a,b):
	mot = str()
	for i in range(0,b):
		mot+=str(a)
	return mot

def formateMot(a,i):
	nbEspace=i-len(a)
	mot=str()
	if(nbEspace>=0):
		mot = str(a)+repeat(" ",nbEspace)
	
	return mot 
		

def createPDB(a):
	with open("PDBout.PDB", "w") as fout:
		for i in a.keys():
			for j in a[i].keys():
				print(formateMot("ATOM", 6)+str(i)+str(a[i][j]['ID']))
				# A compléter

## Prend en entrée deux dictionnaires de protéines et calcule la distance entre chaque atome à la même position
def rmsd(a,b):
	N = 0
	s=0.0
	for j in a.keys():# dans les chaines
		for k in a[j].keys(): # dans les résidus
			for l in a[j][k].keys(): #dans les dictionnaires d'info
				print("a:"+str(a[j][k][l]) +"\n"+ "b:"+str(b[j][k][l]))
				print(pow(distanceAtomes(a[j][k][l],b[j][k][l]),2))
				s+=pow(distanceAtomes(a[j][k][l],b[j][k][l]),2)
				N+=1
	print("n: "+str(N))
	return math.sqrt(s/N)


if __name__ == '__main__':
	monDico = dict()
	ficAtome=sys.argv[1]
	prot1 =sys.argv[2]
	prot2 =sys.argv[3]
	
	dicoProt1 = lirePDB(prot1,ficAtome)
	#~ dicoProt2 = lirePDB(prot2,ficAtome)
	
	print(dicoProt1)
	#Calculer rmsd en premier car on modifie le dictionnaire après
	#~ print(rmsd(dicoProt1,dicoProt2))
	
	#~ ajouterCentreDeMasse(dicoProt1)
	#~ ajouterCentreDeMasse(dicoProt2)
	
	

	
	
	#~ mat = distance(monDico)
	#~ printDistance(mat)

