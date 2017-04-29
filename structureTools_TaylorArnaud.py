#!/usr/bin/env python
#-*- coding : utf3 -*-
import sys
import math


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
def lirePDB(a):

	#### Dictionnaires ####
	residu = dict()
	chaine = dict()
	atome = dict()
	listModeles=list()
	#######################
	
	with open(a, "r") as fichier:
		conf = False
		model = '0' #modèle par défaut = 0
		fic = fichier.readlines()

		for l in fic:
			
			## Il faut rajouter une dimension au dictionnaire afin de prendre en compte la configuration de la protéine
			#Si des modèles différents existent, on les sauvegardent
			if l[0:6].strip()=='MODEL' and l[9:14].strip() not in chaine.keys():
				model = l[9:14].strip()
				
			#Pour 1 modèle donné, on lit une seule configuration (s'il y en a plusieurs)
			if (conf==False and l[0:6].strip()=='ATOM'): #Pour la première ligne du fichier, l'ascii du char[0]=65279 parfois indique le début d'une zone de texte
				conf=l[16]
				
			if (l[0:4]=="ATOM" and l[16]==conf): #Scan des chaines d'interet
				##### Recuperation des donnees ###
				res=l[22:26].strip()
				nomAtome = l[12:16].strip()
				chName=l[22].strip()

				info={'ID': (l[6:11].strip()),
						'x'	: (float)(l[30:38].strip()),
						'y' : (float)(l[38:46].strip()),
						'z' : (float)(l[46:54].strip()),
					}
				######## Fin recuperation ########
				
				# Si modele non répertorié on l'ajoute
				if(model not in chaine.keys()):
					chaine[model] = {}
				
				# Chaine non repertoriee donc on l'ajoute
				if (chName not in chaine[model].keys()):
					chaine[model][chName] = {}
				
				## Residu non encore repertorie
				if(res not in chaine[model][chName].keys()):
					chaine[model][chName][res]={}
				
				chaine[model][chName][res][nomAtome]=info

	return chaine

# 
def selectElement(atomeName):
	atomeName=atomeName.strip()
	courant = ['C','H','O','N','P']
	if atomeName[0] in courant:
		return atomeName[0]
	else:
		print("Un atome non répertorié utilisé !")
	
	
#Ajoute au dictionnaire info (pour chaque atome) la masse de l'atome 'atmW'
def addAtomWeight(pathAtomes, dicoPDB):
	masseAtomes = lireAtoms(pathAtomes)
	for model in dicoPDB.keys():			
		for chaine in dicoPDB[model].keys():
			for residu in dicoPDB[model][chaine].keys():
				for atome in dicoPDB[model][chaine][residu].keys():
					elem = selectElement(atome)
					# Si on a réussi a récupérer le nom de l'atome, on associé une masse à partir du dico atome.txt
					if(elem == None): # Si un atome non répertorié est trouvé on l'affiche
						print(chaine,'	',residu,"	","\'",atome,"\'","	",elem)
					else: # Le nom de l'élément a été récupéré, on ajoute sa masse
						dicoPDB[model][chaine][residu][atome]['atmW']=masseAtomes[elem]
	
	
## Ajoute pour chaque residu du dictionnaire b la position (x,y,z) de son centre de masse
# @a: dictionnaire issu de la commande lirePDB
def ajouterCentreDeMasse(a):
	for mod in a.keys():
		for i in a[mod].keys():			# On parcours la chaine
			for j in a[mod][i].keys():	# On parcours les residus
				
				cdm = {'x':0,'y':0,'z':0}
				masseTotale=0
				for k in a[mod][i][j].keys(): # Extraction des infos pour chaque atome
					## Si la masse atomique est disponible, on fait un calcul plus précis
					if a[mod][i][j][k]['atmW'].exists:
						atmW = a[mod][i][j][k]['atmW']
					else: # Sinon on considère que tous les atomes ont une masse atomique de 1
						atmW = 1
						
					cdm['x']+=atmW*a[mod][i][j][k]['x']
					cdm['y']+=atmW*a[mod][i][j][k]['y']
					cdm['z']+=atmW*a[mod][i][j][k]['z']
					
					masseTotale+=atmW
					
				cdm['x']/=masseTotale
				cdm['y']/=masseTotale
				cdm['z']/=masseTotale
				
				a[mod][i][j]['cdm']=cdm

## prend en entrée deux atomes (donc les dico d'info des deux atomes) (x,y,z) et retourne la distance entre eux
def distanceAtomes(a,b):
	x1 = a['x']
	y1 = a['y']
	z1 = a['z']

	x2 = b['x']
	y2 = b['y']
	z2 = b['z']
	
	dist= math.sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2))
	return dist

## Calcule les distances residu-residu à partir du centre de masse des deux résidus.
# @a: dictionnaire issu de la fonction ajouterCentreDeMasse(a) 
def distance(a):
	mat = dict()
	for mod in a.keys():
		for i in a[mod].keys(): #Chaines i
			for j in a[mod][i].keys(): # Resisud	j
				
				if str(j) not in mat.keys():
					mat[str(j)]={}
				
				#~ for k in a.keys(): # Chaine k
				for l in a[mod][i].keys(): #residu l
				
					if j!=l:
						if str(l) not in mat.keys() or str(j) not in mat[str(l)].keys():
							
							dist=distanceAtomes(a[mod][i][j]['cdm'],a[mod][i][l]['cdm'])
							mat[str(j)][str(l)]={'val':dist}
							
	return mat


## Affiche proprement les distances residu-redisu
# @a: dictionnaire issu de la fonction distance(a)
def printDistance(a):
	for mod in a.keys():
		for i in a[mod].keys():
			for j in a[mod][i].keys():
				print("["+i+"]["+j+"] = "+str(a[mod][i][j]['val']))

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
	with open("PDB_out.PDB", "w") as fout:
		for mod in a.keys():
			fout.write("MODEL",repeat(" ",5),mod)
			for i in a[mod].keys():
				for j in a[mod][i].keys():
					fout.write(formateMot("ATOM", 6)+str(i)+str(a[mod][i][j]['ID']))
					

if __name__ == '__main__':
	monDico = dict()
	ficAtome=sys.argv[3]
	prot1 =sys.argv[1]
	prot2 =sys.argv[2]
	
	#~ dicoProt1 = lirePDB(prot1)
	dicoProt2 = lirePDB(prot2)
	
	#~ ajouterCentreDeMasse(dicoProt1)
	#~ ajouterCentreDeMasse(dicoProt2)

	#~ addAtomWeight(ficAtome,dicoProt2)
	
	#~ print(getConfigs(prot1))
	#~ print(getConfigs(prot2))
	#~ print(dicoProt1.keys())
	#~ mat = distance(monDico)
	#~ print(dicoProt1)
	#~ print(dicoProt2)
	
	
	#useless code
	#~ if nomAtome =='ZN':
		#~ print("Trouvé:",info['ID'])
