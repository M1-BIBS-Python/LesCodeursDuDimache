#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  rmsd.py

import sys
import math
import matplotlib.pyplot as plt
import structureTools_TaylorArnaud as st



## renvoie le RMSD entre deux proteines pour un atome donnée
# @a : dictionnaire de la premire proteine 
# @b : dictionnaire de la deuxieme proteine
# @liste : l'atome presente dans les proteines avec le quel on calcul le RMSD, par defaut la fonction prend le carbone alpha
def simple_rmsd(a,b,liste=['CA']): #Atome par défaut CA
	N = 0
	res = 0
	s=0.0
	for mod in a.keys():
		for dom in a[mod].keys():
			for j in a[mod][dom].keys():# dans les chaines
				for k in a[mod][dom][j].keys(): # dans les résidus
					res+=1
					for l in a[mod][dom][j][k].keys(): #dans les dictionnaires d'info
						l2=st.selectElement(l)
						if(l2 in liste): #On ignore le centre de masse s'il y en a un
							s+=pow(st.distanceAtomes(a[mod][dom][j][k][l],b[mod][dom][j][k][l]),2)
							N+=1
	print("Paires d'atomes alignées : ",N)
	#~ print("Nombre résidus: ",res)
	if(N>0):
		return(math.sqrt(s/N))
	#Sinon aucune paire n'a été alignée : retourne None (car division par 0 impossible)


## renvoie le RMSD entre deux proteines pour un atome donnée
# @a : dictionnaire de la premire proteine 
# @b : dictionnaire de la deuxieme proteine 
# @modA : correspond à la conformation souhaiter de la proteine 1 
# @modB : correspond à la conformation souhaiter de la proteine 2
# @liste : l'atome presente dans les proteines avec le quel on calcul le RMSD, par defaut la fonction prend le carbone alpha
def rmsd(a,b,modA,modB,liste=['CA']): #Atome par défaut CA
	N = 0
	res = 0
	s=0.0
	for dom in a[modA].keys(): # dans les domaines
		for j in a[modA][dom].keys():	# dans les chaines
			for k in a[modA][dom][j].keys(): # dans les résidus
				res+=1
				for l in a[modA][dom][j][k].keys(): # dans l'info
					l2=st.selectElement(l)
					if(l2 in liste): #On ignore le centre de masse s'il y en a un
						s+=pow(st.distanceAtomes(a[modA][dom][j][k][l],b[modB][dom][j][k][l]),2) # On peut écrire ça car même protéine, même chaine/résidus/atomes seul le modèle change
						N+=1
	#~ print("Paires d'atomes alignées : ",N)
	#~ print("Nombre résidus: ",res)
	if(N>0):
		return(math.sqrt(s/N))
	#Sinon aucune paire n'a été alignée : retourne None (car division par 0 impossible)

## calcule le rmsd entre chaque résidus des protéines pour un atome donnée (sortie : colonne 1 position du résidus; colonne 2: RMDS entre les deux conf)
# @a : dictionnaire de la premire proteine 
# @b : dictionnaire de la deuxieme proteine
# @liste : l'atome presente dans les proteines avec le quel on calcul le RMSD, par defaut la fonction prend le carbone alpha
def res_rmsd(a,b,liste=['CA']): #Atome par défaut CA
	nbRes = 0
	dicoRmsd = dict()
	for mod in a.keys(): # dans les modèles
		for dom in a[mod].keys(): # dans les domaines			
			for j in a[mod][dom].keys():# dans les chaines
				for k in a[mod][dom][j].keys(): # dans les résidus
					
					## pour chaque résidu k, on calcul le rmsd des atomes de ce résidu et on le stoque dans un dictionnaire
					nbRes+=1; N = 0; s=0.0;
					
					for l in a[mod][dom][j][k].keys(): #dans les dictionnaires d'info de chaque atome
						l2=st.selectElement(l)
						if(l2 in liste): #Si l'atome fait partie de la liste des atomes sur lesquels on calcule le rmsd:
							s+=pow(st.distanceAtomes(a[mod][dom][j][k][l],b[mod][dom][j][k][l]),2)
							N+=1
							
					#  Si N>0 : Au moins 1 paire d'atome alignée pour le résidu : calcul et stockage du rmsd 
					if(N>0):
						dicoRmsd[float(k)]=float(math.sqrt(s/N)) #Sauvegarde du RMSD du résidu k
				
	print("Nombre résidus: ",nbRes)
	return(dicoRmsd)	

## Fonction qui permettra d'afficher le graphique RMDS vs AA positions (pas de dictionnaire en entrée car temporalité importante) (liste probablement)
# @dic : dictionnaire de RMSD
# @lType : type de tracer, par defaut en ligne 
def drawRMDS(dic,lType='-',title="A Graph Has No Name."):
	listCle =list()
	listVal =list()
	for a in sorted(dic.keys()):
		listCle.append(a)
		listVal.append(dic[a])
	
	print("Je dessine")
	plt.plot(listCle,listVal,lType)
	plt.title(title)
	plt.ylabel('RMSD(A)')
	plt.xlabel('aa positions')
	print("Dessin terminé !")
	plt.show()




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

	#~ #Calculer rmsd
	
	dicTest = res_rmsd(dicoProt1,dicoProt2)
	print(len(dicTest))
	drawRMDS(dicTest)
