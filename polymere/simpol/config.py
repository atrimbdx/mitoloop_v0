#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Locus, ancres, et grille de parametres balayee par les deux voies."""
import itertools

import numpy as np

LOCI = {
    "tyrp1": dict(
        chrom="chr9", debut=10_300_000, fin=13_250_000, taille_bin=5_000,
        marge=500_000,          # flanc simule de part et d'autre
        ancres=[("PTPRD", 10_365_008, 10_380_869),
                ("element_distal", 10_854_755, 10_880_027),
                ("hub_distal", 11_314_973, 11_325_194),
                ("region_enhancer", 12_423_213, 12_459_838),
                ("promoteur_TYRP1", 12_635_660, 12_710_989),
                ("LURAP1L", 12_770_000, 12_780_000),
                ("MPDZ", 13_047_638, 13_125_453)],
        # contacts dont on teste la simultaneite : la chaine gagnee en
        # melanocyte, et le contact perdu
        hub=["element_distal", "hub_distal", "region_enhancer",
             "promoteur_TYRP1"],
        hub_fibro=["promoteur_TYRP1", "MPDZ"],
        encode={"MEL": "foreskin melanocyte", "FIBRO": "fibroblast of dermis"}),
    "oca2": dict(
        chrom="chr15", debut=27_700_000, fin=28_300_000, taille_bin=5_000,
        marge=300_000,
        ancres=[("enhancer_intragenique", 28_019_798, 28_035_447),
                ("promoteur_OCA2", 28_097_411, 28_102_625),
                ("enhancer_HERC2", 28_117_900, 28_133_816)],
        hub=["enhancer_intragenique", "promoteur_OCA2", "enhancer_HERC2"],
        hub_fibro=[],
        encode={"MEL": "foreskin melanocyte", "FIBRO": "fibroblast of dermis"}),
}

# --- voie A : extrusion de boucles ---------------------------------------
GRILLE_A = dict(
    separation_kb=[100, 200, 400],      # une cohesine toutes les N kb
    processivite_kb=[150, 300, 600],    # taille de boucle sans obstacle
    force_barriere=[0.5, 0.9],          # probabilite d'arret sur un site CTCF
)

# --- voie B : liants -----------------------------------------------------
GRILLE_B = dict(
    n_types=[4, 8, 12],
    eps_liant=[0.6, 1.2, 2.0],
)


def jeux(voie):
    """Liste ordonnee des jeux de parametres d'une voie."""
    g = GRILLE_A if voie == "A" else GRILLE_B
    cles = sorted(g)
    return [dict(zip(cles, v)) for v in itertools.product(*(g[c] for c in cles))]


def geometrie(locus):
    """Nombre de billes et decalage, flanc compris."""
    c = LOCI[locus]
    d0 = c["debut"] - c["marge"]
    f0 = c["fin"] + c["marge"]
    n = int(np.ceil((f0 - d0) / c["taille_bin"]))
    return n, d0, f0


def indices_ancres(locus):
    """Indice de bille de chaque ancre, dans le repere simule."""
    c = LOCI[locus]
    n, d0, _ = geometrie(locus)
    return {nom: int(np.clip(((s + e) / 2 - d0) // c["taille_bin"], 0, n - 1))
            for nom, s, e in c["ancres"]}


def surcharger(locus, taille_bin=None, marge=None):
    """Change la resolution ou le flanc d'un locus. Appele depuis les
    scripts par --bin et --marge, ce qui evite de modifier ce fichier pour
    un essai rapide."""
    if taille_bin:
        LOCI[locus]["taille_bin"] = int(taille_bin)
    if marge is not None:
        LOCI[locus]["marge"] = int(marge)
    return LOCI[locus]


# --- loci decrits hors de ce fichier (mitoloop_v0/polymere/preparer_locus.py)
# SIMPOL_LOCI pointe vers un JSON {nom: {chrom, debut, fin, taille_bin, marge,
# ancres: [[nom, debut, fin], ...], hub: [...], hub_fibro: [...]}} ; les loci
# qu'il decrit s'ajoutent a ceux ci-dessus.
import json as _json
import os as _os

if _os.environ.get("SIMPOL_LOCI"):
    with open(_os.environ["SIMPOL_LOCI"]) as _f:
        for _nom, _c in _json.load(_f).items():
            _c["ancres"] = [tuple(x) for x in _c["ancres"]]
            _c.setdefault("hub_fibro", [])
            _c.setdefault("encode", {})
            LOCI[_nom] = _c
