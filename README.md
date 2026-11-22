# bioapp

Application de bureau (Python + PySide6) pour faire de la bio-informatique plus facilement.

## Fonctions
- Chargement de FASTA, téléchargement NCBI (cache local) et recherche NCBI par mot-clé
- Dotplot (fenêtre réglable, option complément inverse)
- Alignement de deux séquences : global, semi-global, local ; ADN (barème modifiable) et protéines (BLOSUM62)
- Comparaison de trois séquences (A/B/C) avec résumé et dotplots par paire
- Significativité par mélange (type PRSS)
- Motifs Prosite (recherche sur les deux brins)
- Traduction, recherche d'ORF, hydropathie de Kyte-Doolittle
- Alignement multiple (Clustal Omega ou MAFFT installé en local)
- Exports : dotplot PNG/SVG, alignement texte et FASTA aligné
- Rich Presence Discord (facultative)

## Installation (Ubuntu / WSL)
    python3 -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt
    sudo apt install clustalo    # pour l'alignement multiple

## Lancer
    ./lancer.sh        # ou : python main.py

## Tests
    python -m pytest

## Exemples
Le dossier `exemples/` contient des FASTA pour tester : alignement manuel du cours (`m1`, `m2`),
thioredoxines (`trx_hp`, `trx_sa`, `trx2`), ORF (`orf`), région hydrophobe (`tm`), trio (`t1`, `t2`, `t3`).
