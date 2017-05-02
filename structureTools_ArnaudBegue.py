#!/usr/bin/env python
#-*- coding : utf8 -*-

import sys
import math
import glob

## Retourne la liste des fichiers comportant le motif recherché
# @a : chemin du dossier comportant les fichiers d'intérêt
# @b : motif voulu, ici le nom du domaine 
def lecture_dossier(a,b):
	path_ref=str(a+"/Refs/"+b) # création du chemin absolue pour accéder au fichier du domaine b dans le dossier ref
	path_frame=str(a+"/Frames/"+b) # création du chemin absolue pour accéder au fichier du domaine b dans le dossier frame
	ref=glob.glob(path_ref) # crée une liste des fichiers contenues dans le dossier ref contenant le motif b dans leur nom
	frame=glob.glob(path_frame)
	return([ref,frame])

def lire_pdb(fichier):

	f= open(fichier, "r")

	data=f.readlines()

	# creation des dictionnaires
	chain={} 

	# creation flag pour identifier la conformation

	flag=0

	for myline in data:
		

		if (myline[0:4] == "ATOM" and flag == 0):
			
			conf=myline[16]
			
			flag=1
			
			
		if (myline[0:4] == "ATOM" and myline[16] == conf):
			
			#si la chaine n'est pas connue on l'a cree
			
			mychain=myline[21].strip() #nom de la chaine			
			if( mychain not in chain.keys()): # si la chaine n'existe pas 			
				chain[mychain]={}
				
			
			myresidu=myline[17:20].strip() # le numero de residu
			if( myresidu not in chain[mychain].keys()): # si le residu n'existe pas 					  			
				chain[mychain][myresidu] = {}
			

			nomatome=myline[13:16].strip() #nom de l'atome	
			info={
				"X":(float)(myline[31:38].strip()),
				"Y":(float)(myline[39:46].strip()),
				"Z":(float)(myline[47:54].strip()),
				"id":myline[7:11].strip()
			}

			chain[mychain][myresidu][nomatome] = info	
			
			f.close()
	
		
	return chain ;

def affichage(chain):
	## Affichage

	affiche= str(chain).replace("},","},\n")
	affiche= str(affiche).replace("},","},\n")		
	affiche=str(affiche).replace("},","},\n")		

	print( affiche) 	

def cbm(chain):

	for i in chain.keys():
		for j in chain[i].keys():
			cbm={'X':0 , 'Y':0, 'Z':0}
			div=0
			for k in chain[i][j].keys():
				
				cbm['X']+=chain[i][j][k]['X']
				cbm['Y']+=chain[i][j][k]['Y']
				cbm['Z']+=chain[i][j][k]['Z']
				div=div+1
			
			cbm['X']/=div	
			cbm["Y"]/=div
			cbm["Z"]/=div
			chain[i][j]["cbm"]=cbm
			
	return chain

def ecrire_PDB(fichier,chain):
	fic=open(fichier, "w")
	for i in chain.keys():
		for j in chain [i].keys():
			for k in chain [i][j].keys():
				
				
				a=str("ATOM"+"\t"+i+"\t"+j+"\t"+k+"\t"+str(chain[i][j][k]["X"])+"\t"+str(chain[i][j][k]["Y"])+"\t"+str(chain[i][j][k]["Z"])+"\n")
				print( a )
				fic.write(a )
	
	fic.close()
	
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

def RMSD_atom(a,b):
	som=0
	N=0
	for i in a.keys():
		for j in a[i].keys():
			for k in a[i][j].keys():				
				# realise la somme des distances 
				som=som+pow(distance(a[i][j][k],b[i][j][k]),2)
				# compte le nombre d'atome 
				N=N+1
	# calcul le RMSD
	RMSD=math.sqrt(som/N)
	print("Nombre d'atomes="+str(N))
	return RMSD		
	
## Retourne la liste des fichiers comportant le motif recherché
# @a : chemin du dossier comportant les fichiers d'intérêt
# @b : motif voulu, ici le nom du domaine 
def lecture_dossier(a,b):
	path_ref=str(a+"/Refs/"+b) # création du chemin absolue pour accéder au fichier du domaine b dans le dossier ref
	path_frame=str(a+"/Frames/"+b) # création du chemin absolue pour accéder au fichier du domaine b dans le dossier frame
	ref=glob.glob(path_ref) # crée une liste des fichiers contenues dans le dossier ref contenant le motif b dans leur nom
	frame=glob.glob(path_frame)
	return([ref,frame])
		
					
def distance(a,b):
	# correspond aux coordonnee de l'atome de la chaine 1 
	x1=a['X']
	y1=a['Y']
	z1=a['Z']
	# correspond aux coordonnee de l'atome de la chaine 2
	x2=b['X']
	y2=b['Y']
	z2=b['Z']

	dist=math.sqrt(pow((x1-x2),2)+pow((y1-y2),2)+pow((z1-z2),2))
	return dist
	
def RMSD_residu(a,b,atome):
	som=0
	N=0
	for i in a.keys():
		for j in a[i].keys():						
			# realise la somme des distances 
			som=som+pow(distance(a[i][j][atome],b[i][j][atome]),2)
			# compte le nombre d'atome 
			N=N+1

	# calcul le RMSD
	RMSD=math.sqrt(som/N)
	print("Nombre de residus ="+str(N))
	return RMSD	

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
							for cahin2 in b[mod2][dom2].keys(): # cahin
								for res2 in b[mod2][dom2][cahin2].keys(): # residu 
									path1= "1"+mod+dom+chain+str(res)
									path2="2"+mod2+dom2+cahin2+str(res2)
									c=contact_residu(b[mod2][dom2][cahin2][res2], a[mod][dom][chain][res])
									if(c==1 and path1 not in contact.keys() ):
										contact[path1]={}
										
									if(c==1 and path2 not in contact[path1].keys()):
										contact[path1][path2]=c
	return(contact)


#~ def graphe(chain):

	
if __name__=='__main__':
	a=lire_pdb(sys.argv[1])
	b=lire_pdb(sys.argv[2])
	a=cbm(a)
	b=cbm(b)
	c=RMSD_residu(a,b,"CA")
	print(c)
	
	
	
