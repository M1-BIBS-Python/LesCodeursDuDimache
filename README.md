# LesCodeursDuDimanche
Projet sRNP H/ACA de BEGUE ARNAUD et VINGADASSALON TAYLOR

Python 3

******************

Cette API permet d'identier les résidus en interaction avec une interface

- structureTools_TaylorArnaud.py contient plusieur outils necessaire à globalCC et LocalCC
- rmsd.py contient les outils permettant de calculer le rmsd
- globalCC.py permet de réaliser une recherche des domaines flexibles dans la proteine
- localCC.py permet à partir des domaines determiner avec globalCC, les résidus en contact avec l'interface  

****************
Usage:

	globalCC nom_proteine_de_reference nom_proteine_dynamique [fichier_de_masse_atomique] 
	
	localCC chemin_dossier_parent_Ref&Frames
	 
