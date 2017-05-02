#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  localCC.py
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
#  along with this program if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
import structureTools_TaylorArnaud as st

if __name__ == '__main__':
 
 
 ########################################	
	#~ prot=sys.argv[1]
	ct=dict()
	
	## Pour chaque conformation données : de 0 à 5000 (automatiser la lecture des confs)
	for i in range(0,100,10) :
		l=st.lecture_dossier("./","A?_"+str(i)+".PDB") #Lecture dossier ref
		l2=st.lecture_dossier("./","B_"+str(i)+".PDB") #Lecture dossier dyn
		
		# Fusion des listes pour faire 
		l2=l2[0]+l2[1] 
		l=l[0]+l[1]	
		
		#~ print(l2)
		#~ print(l)
		## OK
		
		
		#### On a la liste des fichiers (A1,A2,A3,A4) à comparer au fichier ref(B)
		
		p2=st.lirePDB(l2[0]) # Lecture du PDB de référence
		st.cdm(p2)			 # ok 
		
		for p1 in l: # On parcours la liste des fichiers à lire
			#~ print(p1) #ok
			
			p1= st.lirePDB(p1)
			st.cdm(p1)
			st.contact(p1,p2,ct)
			
	#~ print(ct)
		
	print(len(ct))
	
	
	
	###########
	for i in ct.keys():
		#~ print(i,"\n")
		for j in ct[i].keys():
			
			ct[i][j]= {'freq':0, 'tps':0, 'val':ct[i][j]}
			ct[i][j]['freq'] = (ct[i][j]['val'])/500
			ct[i][j]["tps"]=(ct[i][j]['val'])*(10/500)
		
		
	for l in ct.keys():
		for j in ct[l].keys():
			if int(j[2:4]) in list(range(23,38)):
				print(l," : ",j)
	
	#~ print(ct)
	
	
			
