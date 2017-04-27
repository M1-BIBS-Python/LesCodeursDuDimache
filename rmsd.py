#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  rmsd.py
#  
#  Copyright 2017 Taylor <taylor@taylor-PC>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
import sys
import math
import structureTools_TaylorVingadassalon as st


## Prend en entrée deux dictionnaires de protéines et calcule la distance entre chaque atome à la même position
def rmsd(a,b,liste=['CA']): #Atome par défaut CA
	N = 0
	res =0
	s=0.0
	for j in a.keys():# dans les chaines
		for k in a[j].keys(): # dans les résidus
			res+=1
			for l in a[j][k].keys(): #dans les dictionnaires d'info
				if(l in liste): #On ignore le centre de masse s'il y en a un
					s+=pow(st.distanceAtomes(a[j][k][l],b[j][k][l]),2)
					N+=1
	print("Paires d'atomes alignées : ",N)
	print("Nombre résidus: ",res)
	print("RMSD: ",math.sqrt(s/N))
		

def usage():
	print("Problème lors de l'appel de la fonction rmsd: ")

if __name__ == '__main__':
	monDico = dict()
	
	#faire un check sur les paramètres (si le bon nombre n'est pas entré ou que les fichiers n'existent pas, on retourne usage().
	ficAtome=sys.argv[1]
	prot1 =sys.argv[2]
	prot2 =sys.argv[3]
	
	# Récupérer la fonction lirePDB d'arnaud pour cette partie car on a pas besoin de l'info sur les atomes.
	# Ou ajouter à la fonction lirePDB une valeur par défaut genre : useAtomWeight = FALSE par défaut. Si True, on utilise le fichier atome.txt
	dicoProt1 = st.lirePDB(prot1,ficAtome)
	dicoProt2 = st.lirePDB(prot2,ficAtome)
	
	#Calculer rmsd 
	rmsd(dicoProt1,dicoProt2,['CA','N','O'])
	
