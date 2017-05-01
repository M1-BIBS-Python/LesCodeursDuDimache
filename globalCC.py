import os, sys
import rmsd as r
import structureTools_TaylorArnaud as st

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
	
	#########################################
	# Etape 1 : RMSD de chaque conformation # 
	#########################################
	
	#~ # Lecture des PDBs	
	dico_structure_Ref = st.lirePDB(prot1)
	dico_dynamique = st.lirePDB(prot2)
	
	# Création des listes des objets à traiter	
	model_ref = list(dico_structure_Ref.keys()).pop()
	res=dict()
	for model_dyn in dico_dynamique.keys():
		res[model_dyn] = r.rmsd(dico_structure_Ref,dico_dynamique,model_ref,model_dyn)
	
	#~ #Visualisation avec un nuage de points = Mieux
	r.drawRMDS(res,'ro')
	
	#~ st.createPDBMultiThreads(dicoProt1,"Refs")
	#~ st.createPDBMultiThreads(dicoProt2,"Frames")
	
	#########################################
	##### Etape 2 : RMSD des 5 domaines #####
	#########################################
	refId =list()
	refNames = os.listdir("./Refs")
	#Récupération des noms de domaines possibles
	for a in refNames:
		if(a[0:2] not in refId):
			refId.append(a[0:2])

	dicConfs = dict()
	# Parcours des domaines possibles
	for dom in refId: 
		#Extraction des conformations possibles pour le domaine d'intérêt
		currentDom=st.lecture_dossier(".",dom+'*') 	# Récupération du chemin de tous les fichiers associé au domaine dom
		confRef = st.lirePDB(currentDom[0][0]) 		# Chemin de la conformation de référence
		y = list(confRef.keys())
		y = int(y.pop())
			
		for dyn in currentDom[1]: 			# Parcours des fichiers de conformation
			fic = st.lirePDB(dyn)
			
			z = int(list(fic.keys()).pop()) # Récupération du numéro de modèle
			#~ z = int(z.pop())
			if z not in dicConfs.keys():
				dicConfs[z] ={} 			# Sauvegarde nouvelle configuration
				
			w = r.rmsd(confRef,fic,y,z)		# Calcul du RMSD pour chaque configuration
			if(w != None):					# Si calcul possible (CA alignés) 
				dicConfs[z][dom]= w			# Sauvegarde
	
	############################################
	####### Ecriture du fichier de sortie ######
	############################################
	toDraw = dict()
	with open("ficOut","w") as fic: #
		for conf in dicConfs.keys():
			fic.write("MODEL      "+str(conf)+"\n")	
			fic.write("G-RMSD     "+str(res[conf])+"\n") # Ecriture du RMSD global
			
			for dom in dicConfs[conf].keys():  
				fic.write(str(dom)+"\t"+str(dicConfs[conf][dom])+"\n")   # Ecriture du RMSD pour les 5 domaine de la configuration

				if dom not in toDraw.keys():
					toDraw[dom]={}#Dico pour la visualisation
				
				toDraw[dom][conf]=dicConfs[conf][dom] #Ajoute du RMSD pour chaque couple (domaine,config)
	
	
	####################################
	##### Représentation graphique #####
	####################################
	for dom in toDraw.keys():
		r.drawRMDS(toDraw[dom],'ro',str(dom))
		
		
