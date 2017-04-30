#!/usr/bin/env python
#-*- coding : utf3 -*-
import sys
import math
import glob #gestion dossieret fichier


## Retourne le dictionnaire atome-masse moleculaire
# @a : chemin du fichier de donnees atome-masse moleculaire
def lireAtoms(a):
	if(a==None):
		print("Pour utiliser les masses atomiques, vous devez fournir le fichier atomes.txt")
		return None
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
				model = int(l[9:14].strip())
				
			## Si la molécule est une molécule d'eau ou fait partie du milieu, on ne fait rien
			# Dans l'idéal il serait bien d'avoir un tableau des acides aminés possibles comme ça on supprime pas directement TIP, CLA et POT
			if(l[17:20]=='TIP' or l[17:20]=='CLA' or l[17:20]=='POT'):
				continue
			
			#Pour 1 modèle donné, on lit une seule configuration (s'il y en a plusieurs)
			if (conf==False and l[0:6].strip()=='ATOM'): #Pour la première ligne du fichier, l'ascii du char[0]=65279 parfois indique le début d'une zone de texte
				conf=l[16]
				
			if (l[0:4]=="ATOM" and l[16]==conf): #Scan des chaines d'interet
				##### Recuperation des donnees ###
				res=int(l[22:26].strip()) # Residue sequence number
				nomAtome = l[12:16].strip()
				chName=l[22].strip()

				info={'ID': (l[6:11].strip()),
						'x'	: (float)(l[30:38].strip()),
						'y' : (float)(l[38:46].strip()),
						'z' : (float)(l[46:54].strip()),
						'dom': l[72:76].strip(), # Ajoute complexité mémoire (+180Mo sur le gros jeu de données) Mais plus rapide
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

# Converti l'atomeName en nom de l'élément correspondant
def selectElement(atomeName):
	atomeName=atomeName.strip()
	courant = ['C','H','O','N','P','S']
	other = ['ZN', 'FE']
	if atomeName[0] in courant:
		return atomeName[0]
	elif atomeName[0:2] in other:
		return atomeName[0:2]
	else:
		print("Un atome non répertorié utilisé !")
	
	
#Ajoute au dictionnaire info (pour chaque atome) la masse de l'atome 'atmW'
def addAtomWeight(pathAtomes, dicoPDB):
	masseAtomes = lireAtoms(pathAtomes)
	if(masseAtomes==None):
		return None
	for model in dicoPDB.keys():			
		for chaine in dicoPDB[model].keys():
			for residu in dicoPDB[model][chaine].keys():
				for atome in dicoPDB[model][chaine][residu].keys():
					#Si l'utilisateur a rajouté un centre de masse on ne fait rien
					if atome == 'cdm':
						continue
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
					if('atmW' in a[mod][i][j][k]):
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

def formateMot(a,i,alignement="R"):
	nbEspace=i-len(str(a))
	mot=str()
	if(nbEspace>=0):
		if(alignement == 'L'):
			mot = str(a)+repeat(" ",nbEspace)
		elif(alignement== 'R'):
			mot = repeat(" ",nbEspace)+str(a)
		else:
			print("Problème d'alignement du mot !")
	
	return mot 
		
## Ok ça marche mais la sortie n'est pas encore formatée !
def createPDB(a):
	with open("PDB_out.PDB", "w") as fout:
		for mod in sorted(a.keys()): # Dans les modèles
			fout.write(str("MODEL"+repeat(" ",5)+mod+"\n")) #On écrit le modèle qu'on est entrain de traiter
			for i in sorted(a[mod].keys()): # clés des chaines: colonne 22
				for j in sorted(a[mod][i].keys()): # clés des résidus : colonne 22-26
					for d in sorted(a[mod][i][j].keys()): # Clés des dicoInfo (l'ID est un élément d'info) : colonne 12 à 16
						#ID: colonne 6:11
						#x : 30 à 38
						#y : 38 à 46
						#z : 46 à 54
						#dom: 70 à 75
						
						#~ print ("[",mod,"]","[",i,"]","[",j,"]","[",d,"]")
						#~ print(a[mod][i][j][d])
						
						#~ if(d!='cdm'): # S'il ne s'agit pas du centre de masse
						
						
						textLine = (formateMot("ATOM", 6, alignement='L')+formateMot(str(a[mod][i][j][d]['ID']),5)+formateMot(d,4,alignement='L')+repeat(" ",6)+i+formateMot(j,4)+repeat(" ",4)+
									formateMot(a[mod][i][j][d]['x'],8)+formateMot(a[mod][i][j][d]['y'],8)+formateMot(a[mod][i][j][d]['z'],8)+repeat(" ",16)+
									formateMot(a[mod][i][j][d]['dom'],5)+"\n")
						
						fout.write(textLine)
						
					

## Retourne la liste des fichiers comportant le motif recherché
# @a : chemin du dossier comportant les fichiers d'intérêt
# @b : motif voulu, ici le nom du domaine 
def lecture_dossier(a,b):
	path_ref=str(a+"/Refs/"+b) # création du chemin absolue pour accéder au fichier du domaine b dans le dossier ref
	path_frame=str(a+"/Frames/"+b) # création du chemin absolue pour accéder au fichier du domaine b dans le dossier frame
	ref=glob.glob(path_ref) # crée une liste des fichiers contenues dans le dossier ref contenant le motif b dans leur nom
	frame=glob.glob(path_frame)
	return([ref,frame])



if __name__ == '__main__':
	monDico = dict()
	
	if (len(sys.argv) == 4):
		prot1 =sys.argv[1]
		prot2 =sys.argv[2]
		ficAtome=sys.argv[3]
		print("Fichier ",ficAtome," fourni comme fichier d'atomes !")
	elif (len(sys.argv)== 3):
		prot1 =sys.argv[1]
		prot2 =sys.argv[2]
		print("Default mode !")
	else:
		print("Le format d'entrée attendu est : structureTools_TaylorArnaud.py fichier_1.PDB fichier_2.PDB [atomes.txt]")
		exit()
	
	dicoProt1 = lirePDB(prot1)
	dicoProt2 = lirePDB(prot2)
	print(dicoProt1)
	#~ for mod in sorted(dicoProt1.keys()):
		#~ for chaine in sorted(dicoProt1[mod].keys()):
			#~ for residu in sorted(dicoProt1[mod][chaine].keys()): 
				#~ # Normalement les résidus triés dans l'ordre des ID de résidus
				#~ # On s'arrête à ce niveau de trie : car les domaines sont normalement triés aussi
				#~ print(mod," ",chaine," ",residu," ",dicoProt1[mod][chaine][residu])
				
	
	#~ print(dicoProt1)
	#~ print(dicoProt2)
	
	#~ addAtomWeight(ficAtome,dicoProt1)
	#~ addAtomWeight(ficAtome,dicoProt2)
	
	#~ ajouterCentreDeMasse(dicoProt1)
	#~ ajouterCentreDeMasse(dicoProt2)
	
	
	#~ print(dicoProt2)
	
	createPDB(dicoProt1)
	
	
	#~ mat = distance(monDico)
	#~ printDistance(mat)
	
	


