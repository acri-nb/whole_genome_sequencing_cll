import re


def is_chrom(chr):
    return re.match(r"^chr\d+", chr)



def report_falco_val(expression,fichier):
    valeur = None
    with open(fichier, 'r') as f:
        lignes = f.readlines()
        for ligne in lignes:
            if expression in ligne:
                valeur = ligne.split()[0]
                break
    return valeur