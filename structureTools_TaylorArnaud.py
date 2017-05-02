#!/usr/bin/env python
#-*- coding : utf3 -*-
import os,sys
import math
import glob, shutil #gestion dossier et fichier


#optimisation multithreading
import itertools
from multiprocessing.dummy import Pool as ThreadPool 


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
				dom = l[72:76].strip()
				
				info={'ID': (l[6:11].strip()),
						'x'	: (float)(l[30:38].strip()),
						'y' : (float)(l[38:46].strip()),
						'z' : (float)(l[46:54].strip()),
					}
				######## Fin recuperation ########
				
				# Si modele non répertorié on l'ajoute
				if(model not in chaine.keys()):
					chaine[model] = {}
				
				# Si domaine non répertorié
				if(dom not in chaine[model].keys()):
					chaine[model][dom] = {}
				
				# Chaine non repertoriee donc on l'ajoute
				if (chName not in chaine[model][dom].keys()):
					chaine[model][dom][chName] = {}
				
				## Residu non encore repertorie
				if(res not in chaine[model][dom][chName].keys()):
					chaine[model][dom][chName][res]={}
				
				chaine[model][dom][chName][res][nomAtome]=info

	return chaine

## Converti l'atomeName en nom de l'élément correspondant
# @atomeName : nom de l'atome extrait du pdb
def selectElement(atomeName):
	atomeName=atomeName.strip()
	courant = ['C','H','O','N','P','S']
	other = ['ZN', 'FE','CA']
	
	if atomeName[0:2] in other:
		return atomeName[0:2]	
	elif atomeName[0] in courant:
		return atomeName[0]
	
	else:
		print("Un atome non répertorié utilisé ! ",atomeName)
		
##Ajoute au dictionnaire info (pour chaque atome) la masse de l'atome 'atmW'
# @pathAtomes : 
# @dicoPDB :
def addAtomWeight(pathAtomes, dicoPDB):
	masseAtomes = lireAtoms(pathAtomes)
	if(masseAtomes==None):
		return None
	for model in dicoPDB.keys():
		for dom in dicoPDB[model].keys():		
			for chaine in dicoPDB[model][dom].keys():
				for residu in dicoPDB[model][dom][chaine].keys():
					for atome in dicoPDB[model][dom][chaine][residu].keys():
						#Si l'utilisateur a rajouté un centre de masse on ne fait rien
						if atome == 'cdm':
							continue
						elem = selectElement(atome)
						# Si on a réussi a récupérer le nom de l'atome, on associé une masse à partir du dico atome.txt
						if(elem == None): # Si un atome non répertorié est trouvé on l'affiche
							print(chaine,'	',residu,"	","\'",atome,"\'","	",elem)
						else: # Le nom de l'élément a été récupéré, on ajoute sa masse
							dicoPDB[model][dom][chaine][residu][atome]['atmW']=masseAtomes[elem]
		
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
		for dom in a[mod].keys():
			for i in a[mod][dom].keys(): #Chaines i
				for j in a[mod][dom][i].keys(): # Resisud	j
					
					if str(j) not in mat.keys():
						mat[str(j)]={}
					
					#~ for k in a.keys(): # Chaine k
					for l in a[mod][dom][i].keys(): #residu l
					
						if j!=l:
							if str(l) not in mat.keys() or str(j) not in mat[str(l)].keys():
								dist=distanceAtomes(a[mod][dom][i][j]['cdm'],a[mod][dom][i][l]['cdm'])
								mat[str(j)][str(l)]={'val':dist}
								
	return mat

## Affiche proprement les distances residu-redisu
# @a: dictionnaire issu de la fonction distance(a)
def printDistance(a):
	for mod in a.keys():
		for dom in a[mod].keys():
			for i in a[mod][dom].keys():
				for j in a[mod][dom][i].keys():
					print("["+i+"]["+j+"] = "+str(a[mod][dom][i][j]['val']))

## Retourne le motif "a" repeter "b" fois
# @a : motif à repeter
# @b : nombre d'occurence du motif "a" souhaiter 
def repeat(a,b):
	mot = str()
	for i in range(0,b):
		mot+=str(a)
	return mot

## retourne le mot "a" compléter avec des espaces pour atteindre "i" caractère
# @a : mot à formater
# @i : nombre de caractere souhaiter
# @alignement : détermine l'ajout des espaces, de base "R" à droite du mot, sinon "L" pour ajouter à gauche du mot
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

## Crée un fichier PDB par modèle et par domaine à partir d'un fichier pdb de base contenant (ou pas) plusieurs modèles et domaines
# @a : dictionnaire contenant le fichier pdb de base 
# @b : nom du dossier dans lequel stocker les fichiers créés
def createPDB(a,dossier=""): #OBSOLETE ?!
	for mod in sorted(a.keys()): # Dans les modèles
		for dom in sorted(a[mod].keys()):
			#A chaque changement de domaine, le fichier d'écriture change
			with open(str(dossier)+"/"+str(dom)+"_"+str(mod)+".PDB", "w") as fout:
				fout.write(formateMot("MODEL",6)+"    "+formateMot(mod,6,alignement='L')+"\n")
				for i in sorted(a[mod][dom].keys()): # clés des chaines: colonne 22
					for j in sorted(a[mod][dom][i].keys()): # clés des résidus : colonne 22-26
						#On ne trie pas sur les noms d'atomes
						for d in sorted(a[mod][dom][i][j].keys()): # Clés des dicoInfo (noms atomes) (l'ID est un élément d'info) : colonne 12 à 16
							
							textLine = (formateMot("ATOM", 6, alignement='L')+formateMot(str(a[mod][dom][i][j][d]['ID']),5)+" "+formateMot(d,4,alignement='L')+repeat(" ",6)+i+formateMot(j,4)+repeat(" ",4)+
									 	formateMot(a[mod][dom][i][j][d]['x'],8)+formateMot(a[mod][dom][i][j][d]['y'],8)+formateMot(a[mod][dom][i][j][d]['z'],8)+
										formateMot(a[mod][dom][i][j][d],5)+repeat(" ",18)+formateMot(dom,3,alignement='L')+"\n")
							
							fout.write(textLine)
							
## Crée un fichier PDB par modèle et par domaine à partir d'un fichier pdb de base contenant (ou pas) plusieurs modèles et domaines
## Version optimisée pour utiliser 8 coeurs
# @a : dictionnaire contenant le fichier pdb de base 
# @b : nom du dossier dans lequel stocker les fichiers créés
def createPDBMultiThreads(a,dossier):
	path = dossier+"/"
	if(os.path.exists(path)):
		print('Le dossier existe déjà. Il va être réécrit !')
		shutil.rmtree(path)
		
	print("Création du dossier:",dossier)
	os.mkdir(path);
	
	pool = ThreadPool(4)  # On va utiliser 4 threads (si coeurs virtuels temps presque pareil)
	#Début de l'écriture: 1 modèle sur 1 thread
	pool.starmap(createThread,zip(itertools.repeat(a),itertools.repeat(path),sorted(a.keys())))
	pool.close()
	pool.join()
	
## Fonction helper de createPDB pour un modèle donné, crée les fichiers pdb associés à ce modèle
#@a :dictionnaire sur lequel on travaille
#@path :chemin du dossier de stockage
#@mod : modèle à traiter
def createThread(a,path,mod):
	for dom in sorted(a[mod].keys()):
			#A chaque changement de domaine, le fichier d'écriture change
			with open(str(path)+str(dom)+"_"+str(mod)+".PDB","w") as fout:
				fout.write(formateMot("MODEL",6,alignement='L')+repeat(" ",4)+formateMot(mod,6,alignement='L')+"\n")
				for i in sorted(a[mod][dom].keys()): # clés des chaines: colonne 22
					for j in sorted(a[mod][dom][i].keys()): # clés des résidus : colonne 22-26
						#On ne trie pas sur les noms d'atomes
						for d in sorted(a[mod][dom][i][j].keys()): # Clés des dicoInfo (noms atomes) (l'ID est un élément d'info) : colonne 12 à 16
							
							textLine = (formateMot("ATOM", 6, alignement='L')+formateMot(str(a[mod][dom][i][j][d]['ID']),5)+" "+formateMot(d,4,alignement='L')+repeat(" ",6)+i+formateMot(j,4)+repeat(" ",4)+
										formateMot(a[mod][dom][i][j][d]['x'],8)+formateMot(a[mod][dom][i][j][d]['y'],8)+formateMot(a[mod][dom][i][j][d]['z'],8)+repeat(" ",16)+
										formateMot(a[mod][dom][i][j][d],5)+repeat(" ",5)+dom+"\n")
							
							fout.write(textLine)
					
## Retourne la liste des fichiers comportant le motif recherché
# @a : chemin du dossier comportant les fichiers d'intérêt
# @b : motif voulu, ici le nom du domaine 
def lecture_dossier(a,b):
	path_ref=str(a+"Refs/"+b) # création du chemin absolue pour accéder au fichier du domaine b dans le dossier ref
	path_frame=str(a+"Frames/"+b) # création du chemin absolue pour accéder au fichier du domaine b dans le dossier frame
	ref=glob.glob(path_ref) # crée une liste des fichiers contenues dans le dossier ref contenant le motif b dans leur nom
	frame=glob.glob(path_frame)
	
	return(ref+frame)

## Retourne un dictionnaire "a" avec le centre de masse des residus
# @a : dicionnaire d'un fichier pdb
def cdm(a):

	for i in a.keys():   # parcourt les models
		for j in a[i].keys(): # parcourt les domaines
			for k in a[i][j].keys(): # parcourt les chaines
				for l in a[i][j][k].keys(): # parcourt les residus
					cdm={'x':0 , 'y':0, 'z':0, 'r':0}
					div=0
					for p in a[i][j][k][l].keys(): # parcourt les infos	
						## Si la masse atomique est disponible, on fait un calcul plus précis
						if('atmW' in a[i][j][k][l][p]):
							atmW = a[i][j][k][l][p]['atmW']
						else: # Sinon on considère que tous les atomes ont une masse atomique de 1
							atmW = 1
											
						cdm['x']+=atmW*a[i][j][k][l][p]['x']
						cdm['y']+=atmW*a[i][j][k][l][p]['y']
						cdm['z']+=atmW*a[i][j][k][l][p]['z']
						div=div+1
					
					cdm['x']/=div	
					cdm["y"]/=div
					cdm["z"]/=div
					a[i][j][k][l]["cdm"]=cdm
					# sous l'hyp le cdm est à équi-distant de tout les atomes des residus
					# le rayon du cbm 
					cdm['r']=distanceAtomes(a[i][j][k][l][p],cdm)
					
			
	return a

## Renvoie 1 si il y a contact entre les deux residus sinon 0
# @a: dictionnaire d'un residu contenant un cdm
# @b : dictionnaire d'un residu contenant un cdm
def contact_residu(a,b):
	dist=distanceAtomes(a['cdm'],b['cdm']) # distance entre les deux cdm

	sumR=a['cdm']['r']+b['cdm']['r'] 
				
	if(dist <= sumR):
		contact=1
	else:
		contact=0					
			
	return contact
	
## Renvoie un dictionnaire contenant les résidus en contacts 
# @a: dictionnaire d'une proteirn avec les cdm
# @b : dictionnaire d'une proteine avec les cdm
def contact(a,b):
	contact={}
	# PARCOURT DICT a
	for mod in a.keys(): # model
		for dom in a[mod].keys(): # domaine
			for chain in a[mod][dom].keys(): # chain
				for res in a[mod][dom][chain].keys(): # residu 
					# PARCOURT DICT	b 
					for mod2 in b.keys(): # model
						for dom2 in b[mod2].keys(): # domaine
							for chain2 in b[mod2][dom2].keys(): # cahin
								for res2 in b[mod2][dom2][chain2].keys(): # residu 
									path1= "1"+str(dom)+str(chain)+str(res)
									path2="2"+str(dom2)+str(chain2)+str(res2)
									c=contact_residu(b[mod2][dom2][chain2][res2], a[mod][dom][chain][res])
									if(c==1 and path1 not in contact.keys() ):
										contact[path1]={}	
									if(c==1 and path2 not in contact[path1].keys()):
										contact[path1][path2]=0
									if(c==1 and path2 in contact[path1].keys()):
										contact[path1][path2]+=1
	return(contact)


if __name__ == '__main__':
	#~ monDico = dict()
	
	#~ if (len(sys.argv) == 4):
		#~ prot1 =sys.argv[1]
		#~ prot2 =sys.argv[2]
		#~ ficAtome=sys.argv[3]
		#~ print("Fichier ",ficAtome," fourni comme fichier d'atomes !")
	#~ elif (len(sys.argv)== 3):
		#~ prot1 =sys.argv[1]
		#~ prot2 =sys.argv[2]
		#~ print("Default mode !")
	#~ else:
		#~ print("Le format d'entrée attendu est : structureTools_TaylorArnaud.py fichier_1.PDB fichier_2.PDB [atomes.txt]")
		#~ exit()
	
	#~ os.nice(10) # Le processus reçoit le niveau maximum de priorité !
	
	#~ dicoProt1 = lirePDB(prot1)
	#~ dicoProt2 = lirePDB(prot2)				

	#~ addAtomWeight(ficAtome,dicoProt1)
	#~ addAtomWeight(ficAtome,dicoProt2)
	
	#~ cdm(dicoProt1)
	#~ cdm(dicoProt2)
	
	#~ createPDBMultiThreads(dicoProt1,"Refs")
	#~ createPDBMultiThreads(dicoProt2,"Frames")
	
	#~ mat = distance(monDico)
	#~ printDistance(mat)
	
########################################	
	prot=sys.argv[1]
	print(prot)
	prot2=sys.argv[2]	
	print(2)
	dicoProt2 = lirePDB(prot2)
	print(3)
	dicoProt2=cdm(dicoProt2)
	print(4)
	#~ l=lecture_dossier(prot,"A*")
	#~ print(l)
	#~ for i in l :
	#~ prot1=i				
	#~ print (i)
	dicoProt= lirePDB(prot)
	dicoProt=cdm(dicoProt)
		
	ct=contact(dicoProt,dicoProt2)
	print(ct)	
	print(len(ct))
		
