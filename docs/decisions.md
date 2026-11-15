# Choix de conception

- **Python + PySide6** : développement rapide, accès direct à Biopython, NumPy et matplotlib.
  Rust (PyO3) seulement si un goulot d'étranglement est mesuré.
- **Séparation `core/` et `gui/`** : la logique (parsing, alignement, dotplot) ne dépend
  pas de l'interface, elle est testable en ligne de commande avec pytest.
- **MVP** : charger deux FASTA, afficher un dotplot, aligner (global/local).
- **CI** : les tests de `core/` tournent à chaque push. PySide6 n'est pas installé
  en CI, donc `core/` ne doit jamais importer `gui/` ni PySide6.
- **À prévoir** : cache local des requêtes NCBI.
