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
	conf=""

	with open(a, "r") as fichier:
		fic = fichier.readlines()
		for l in fic:
			
			if (conf=="" and l[0:4]=="ATOM"):
				conf=l[16]
				
			if (l[0:4]=="ATOM" and l[16]==conf): #Scan des chaines d'interet
				##### Recuperation des donnees ###
				res=l[24:31].strip()
				nomAtome = l[13:16].strip()
				chName=l[20:23].strip()
				
				info={'ID': (float)(l[9:12].strip()),
						'x'	: (float)(l[32:39].strip()),
						'y' : (float)(l[40:47].strip()),
						'z' : (float)(l[48:55].strip()),
						'atmW': atmWeight[l[77].strip()]
					}
				######## Fin recuperation ########
				
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
	return

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
						x1 = a[i][j]['cdm']['x']
						y1 = a[i][j]['cdm']['y']
						z1 = a[i][j]['cdm']['z']
						
						x2 = a[i][l]['cdm']['x']
						y2 = a[i][l]['cdm']['y']
						z2 = a[i][l]['cdm']['z']
						
						distance= math.sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2))
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


def rmsd(a,b):
	N = int(); #Admettons qu'on possède déja N
	for i in range(N):
		print("test"+i)



if __name__ == '__main__':
	monDico = dict()
	monDico = lirePDB(sys.argv[1],sys.argv[2])
	
	ajouterCentreDeMasse(monDico)
	
	mat = distance(monDico)
	printDistance(mat)

