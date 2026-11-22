"""Scénarios de test manuel : chaque entrée charge A/B, règle la fenêtre et lance une action.
Sans dépendance à Qt : les valeurs attendues sont aussi vérifiées par pytest."""

BAREME_DEFAUT = {"match": 1, "mismatch": -1, "ouverture": -2, "extension": -1}
BAREME_COURS = {"match": 5, "mismatch": -4, "ouverture": -10, "extension": -0.5}

TRX_HP = "MSHYIELTEENFESTIKKGVALVDFWAPWCGPCKMLSPVIDELASEYEGKAKICKVNTDEQEELSAKFGIRSIPTLLFTKDGEVVHQLVGVQTKVALKEQLNKLLG"
TRX_SA = "MAIVKVTDADFDSKVESGVQLVDFWATWCGPCKMIAPVLEELAADYEGKADILKLDVDENPSTAAKYEVMSIPTLIVFKDGQPVDKVVGFQPKENLAEVLDKHL"

SCENARIOS = [
    {"titre": "1. Alignement manuel du cours (semi-global)",
     "A": ("m1", "AAGTCATTCGCACATCG"), "B": ("m2", "AACACATCG"),
     "params": BAREME_COURS, "mode": "semi-global", "action": "aligner", "score": 31.5,
     "attendu": "Score 31.5. Sous AAGTCATTCGCACATCG, l'autre ligne doit être AA--------CACATCG "
                "(un grand gap, les 7 dernières bases alignées)."},
    {"titre": "2. Local : seule la région commune compte",
     "A": ("a2", "TTTTTTTTTTATGCGTACGTTAGCTTTTTTTTTT"),
     "B": ("b2", "GGGGGGGGGGATGCGTACGTTAGCCCCCCCCCC"),
     "params": BAREME_DEFAUT, "mode": "local", "action": "aligner", "score": 14.0,
     "attendu": "Score 14.0 et « Région alignée : A 11–24 | B 11–24 » (14 colonnes identiques)."},
    {"titre": "3. Global : les bords coûtent cher",
     "A": ("a2", "TTTTTTTTTTATGCGTACGTTAGCTTTTTTTTTT"),
     "B": ("b2", "GGGGGGGGGGATGCGTACGTTAGCCCCCCCCCC"),
     "params": BAREME_DEFAUT, "mode": "global", "action": "aligner", "score": -7.0,
     "attendu": "Score -7.0 : mêmes séquences que le scénario 2, mais les bords différents "
                "sont pénalisés."},
    {"titre": "4. Semi-global : une séquence incluse dans une autre",
     "A": ("long", "GGGACGTACGTGGG"), "B": ("court", "ACGTACGT"),
     "params": BAREME_DEFAUT, "mode": "semi-global", "action": "aligner", "score": 8.0,
     "attendu": "Score 8.0 ; la séquence courte est alignée avec des gaps sans coût aux deux "
                "bouts : ---ACGTACGT---. En global, le même couple donne 0.0."},
    {"titre": "5. Protéines : thioredoxines (BLOSUM62)",
     "A": ("trx_hp", TRX_HP), "B": ("trx_sa", TRX_SA),
     "params": BAREME_DEFAUT, "mode": "local", "action": "aligner", "score": 290.0,
     "attendu": "Ligne « protéines : BLOSUM62, gap -10/-1 », score 290.0, identité 52/101."},
    {"titre": "6. Dotplot de deux séquences presque identiques",
     "A": ("t1", "ACGTACGTACGTAAC"), "B": ("t2", "ACGTACGTACGTAAT"),
     "params": BAREME_DEFAUT, "mode": "global", "action": "dotplot", "fenetre": 3,
     "attendu": "Une diagonale principale de 12 cases (la dernière fenêtre contient la "
                "différence), plus des diagonales décalées de 4 en 4 : la séquence est répétée."},
    {"titre": "7. Dotplot contre le complément inverse",
     "A": ("t1", "ACGTACGTACGTAAC"), "B": ("t1_rc", "GTTACGTACGTACGT"),
     "params": BAREME_DEFAUT, "mode": "global", "action": "dotplot", "fenetre": 3, "rc": True,
     "attendu": "La case « B en complément inverse » est cochée : la diagonale principale est "
                "continue (13 cases). Sans la case, elle disparaît."},
]
