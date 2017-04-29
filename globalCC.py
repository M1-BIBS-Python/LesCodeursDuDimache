#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  globalCC.py
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
import rmsd as r
import structureTools_TaylorArnaud as st

# Représenter cette évolution du RMDS (entre la prot de référence et chaque conformation) au cours du temps.


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
	
	# Lecture des PDBs	
	dico_structure_Ref = st.lirePDB(prot1)
	dico_dynamique = st.lirePDB(prot2)
	
	## Etape 1) Calculer RMSD pour chaque conformation de la protéine 
	# Création des listes des objets à traiter	
	model_ref = list(dico_structure_Ref.keys()).pop()
	res=dict()
	for model_dyn in dico_dynamique.keys():
		res[model_dyn] = r.rmsd(dico_structure_Ref,dico_dynamique,model_ref,model_dyn)
	
	#visualisation avec un diagramme en batons
	#~ r.drawBarChart(res)
	#Visualisation avec un nuage de points = Mieux
	r.drawRMDS(res,'ro')

	
