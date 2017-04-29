#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  rmsd.py

import sys
import math
import matplotlib.pyplot as plt
import structureTools_TaylorArnaud as st



## Prend en entrée deux dictionnaires de protéines et calcule le rmsd entre ces deux protéines
def rmsd(a,b,liste=['CA']): #Atome par défaut CA
	N = 0
	res = 0
	s=0.0
	for mod in a.keys():
		for j in a[mod].keys():# dans les chaines
			for k in a[mod][j].keys(): # dans les résidus
				res+=1
				for l in a[mod][j][k].keys(): #dans les dictionnaires d'info
					if(l in liste): #On ignore le centre de masse s'il y en a un
						s+=pow(st.distanceAtomes(a[mod][j][k][l],b[mod][j][k][l]),2)
						N+=1
	print("Paires d'atomes alignées : ",N)
	print("Nombre résidus: ",res)
	if(N>0):
		return(math.sqrt(s/N))
	#Sinon aucune paire n'a été alignée : retourne None (car division par 0 impossible)

## Prend en entrée deux dictionnaires de protéines et calcule le rmsd entre chaque résidus des protéines 
## (sortie : colonne 1 position du résidus; colonne 2: RMDS entre les deux conf)
def res_rmsd(a,b,liste=['CA']): #Atome par défaut CA

	nbRes = 0
	dicoRmsd = dict()
	for mod in a.keys():
		for j in a[mod].keys():# dans les chaines
			for k in a[mod][j].keys(): # dans les résidus
				
				## pour chaque résidu k, on calcul le rmsd des atomes de ce résidu et on le stoque dans un dictionnaire
				nbRes+=1; N = 0; s=0.0;
				
				for l in a[mod][j][k].keys(): #dans les dictionnaires d'info de chaque atome
					if(l in liste): #Si l'atome fait partie de la liste des atomes sur lesquels on calcule le rmsd:
						s+=pow(st.distanceAtomes(a[mod][j][k][l],b[mod][j][k][l]),2)
						N+=1
						
				#  Si N>0 : Au moins 1 paire d'atome alignée pour le résidu : calcul et stockage du rmsd 
				if(N>0):
					dicoRmsd[float(k)]=float(math.sqrt(s/N)) #Sauvegarde du RMSD du résidu k
			
	print("Paires d'atomes alignées : ",N)
	print("Nombre résidus: ",nbRes)
	return(dicoRmsd)	

## Fonction qui permettra d'afficher le graphique RMDS vs AA positions (pas de dictionnaire en entrée car temporalité importante) (liste probablement)
def drawRMDS(a):
	listCle =list()
	listVal =list()
	for a in dicTest.keys():
		listCle.append(a)
		listVal.append(dicTest[a])
	
	print("Je dessine")
	plt.plot(listCle,listVal)
	plt.ylabel('RMSD(A)')
	plt.xlabel('aa positions')
	print("Dessin terminé !")
	plt.show()
	#A compléter cf séance 8 partie 4

def usage():
	print("Problème lors de l'appel de la fonction rmsd: ")

if __name__ == '__main__':
	monDico = dict()
	
	#faire un check sur les paramètres (si le bon nombre n'est pas entré ou que les fichiers n'existent pas, on retourne usage().
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
		
	
	# Récupérer la fonction lirePDB d'arnaud pour cette partie car on a pas besoin de l'info sur les atomes.
	# Ou ajouter à la fonction lirePDB une valeur par défaut genre : useAtomWeight = FALSE par défaut. Si True, on utilise le fichier atome.txt
	dicoProt1 = st.lirePDB(prot1)
	dicoProt2 = st.lirePDB(prot2)
	
	#Calculer rmsd 
	dicTest = res_rmsd(dicoProt1,dicoProt2)
		
	print(dicTest.keys())
	print(dicTest.values())
	print(len(dicTest))
	drawRMDS(dicTest)
