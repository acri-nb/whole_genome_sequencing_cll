import os
import shutil

# Chemin du répertoire "jmazurie"
repertoire_jmazurie = "/mnt/project3/jmazurie"

# Chemin du sous-répertoire "bam_05"
sous_repertoire_bam_05 = "/mnt/project3/jmazurie/bam_05"
sous_repertoire_bam_1 = "/mnt/project3/jmazurie/bam_1"

# Parcourir les fichiers dans le répertoire jmazurie
for nom_fichier in os.listdir(repertoire_jmazurie):
    chemin_fichier = os.path.join(repertoire_jmazurie, nom_fichier)
    
    # Vérifier si le fichier porte l'expression "_0.5"
    if "_0.5" in nom_fichier:
        # Déplacer le fichier dans le sous-répertoire bam_05
        nouveau_chemin = os.path.join(sous_repertoire_bam_05, nom_fichier)
        shutil.move(chemin_fichier, nouveau_chemin)
        print(f"Le fichier {nom_fichier} a été déplacé vers {sous_repertoire_bam_05}")
        
    if "_1" in nom_fichier:
        # Déplacer le fichier dans le sous-répertoire bam_05
        nouveau_chemin = os.path.join(sous_repertoire_bam_1, nom_fichier)
        shutil.move(chemin_fichier, nouveau_chemin)



