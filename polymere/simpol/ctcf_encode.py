#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Voie A, entree — sites CTCF orientes, depuis ENCODE.

Faute de ChIP-seq CTCF dans les cellules de l'etude, les jeux ENCODE des
types cellulaires les plus proches servent d'approximation : melanocyte de
prepuce pour le melanocyte, fibroblaste dermique pour le fibroblaste. Le
texte d'un article doit le dire, c'est un substitut et non une mesure.

L'orientation, qui est ce qui compte pour l'extrusion, n'est pas dans les
pics : elle vient du balayage du motif CTCF sous chaque pic. Le motif
utilise est MA0139.1 de JASPAR, dix-neuf positions, recopie ici pour
eviter une dependance reseau de plus.

SORTIE  un BED : chrom, debut, fin, nom, score du motif, brin, signal du pic
"""
import argparse
import gzip
import io
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np

# JASPAR MA0139.1 — comptes A, C, G, T par position
MA0139 = np.array([
    [87, 291, 76, 459], [167, 145, 414, 187], [281, 49, 449, 134],
    [56, 800, 21, 36], [8, 903, 0, 2], [744, 13, 65, 91],
    [40, 528, 334, 11], [107, 433, 48, 324], [851, 11, 32, 18],
    [5, 0, 903, 5], [333, 3, 566, 11], [54, 12, 504, 343],
    [12, 0, 890, 11], [56, 8, 775, 74], [104, 733, 5, 71],
    [372, 13, 507, 21], [82, 482, 307, 42], [117, 322, 73, 401],
    [402, 181, 266, 64],
], dtype=float)

BASES = {"A": 0, "C": 1, "G": 2, "T": 3}


def pwm_log_odds(comptes, pseudo=0.8):
    f = (comptes + pseudo) / (comptes.sum(1, keepdims=True) + 4 * pseudo)
    return np.log2(f / 0.25)


def encoder(seq):
    x = np.full(len(seq), -1, dtype=np.int8)
    for b, k in BASES.items():
        x[np.frombuffer(seq.upper().encode(), dtype="S1") == b.encode()] = k
    return x


def balayer(seq, pwm):
    """Meilleur score et brin sur une sequence. Renvoie (score, brin, pos)."""
    L = len(pwm)
    x = encoder(seq)
    if len(x) < L:
        return -np.inf, ".", 0
    rc = pwm[::-1, ::-1]                    # complement inverse du motif
    best = (-np.inf, ".", 0)
    idx = np.arange(L)
    for i in range(len(x) - L + 1):
        f = x[i:i + L]
        if (f < 0).any():
            continue
        for mat, brin in ((pwm, "+"), (rc, "-")):
            s = float(mat[idx, f].sum())
            if s > best[0]:
                best = (s, brin, i)
    return best


# ======================================================================
BASE = "https://www.encodeproject.org"


def _api(chemin, **params):
    """Appel JSON. Le /search/ d'ENCODE repond 404 quand rien ne
    correspond : ce n'est pas une panne, c'est un resultat vide, et le
    confondre avec une erreur faisait echouer tout le script sur une
    simple absence de donnees."""
    url = BASE + chemin + "?" + urllib.parse.urlencode(params, doseq=True)
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            return json.load(r), url
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return {"@graph": [], "total": 0}, url
        raise


def _fichier_pics(acc, assemblage):
    """Meilleur BED de pics d'une experience, ou None."""
    # format=json : sans lui, le portail renvoie la page HTML de l'experience
    det, _ = _api(f"/experiments/{acc}/", format="json")
    cands = [f for f in det.get("files", [])
             if f.get("file_format") == "bed"
             and f.get("assembly") == assemblage
             and f.get("status") == "released"
             and "peaks" in (f.get("output_type") or "")]
    # on prefere les pics IDR conservatifs, les plus propres
    cands.sort(key=lambda f: ("conservative" not in f["output_type"],
                              "IDR" not in f["output_type"]))
    return BASE + cands[0]["href"] if cands else None


def biosamples_disponibles(cible="CTCF", n=400):
    """Types cellulaires pour lesquels ENCODE a un ChIP-seq de la cible."""
    d, _ = _api("/search/", type="Experiment", **{
        "target.label": cible, "status": "released",
        "format": "json", "limit": n,
        "field": ["accession", "biosample_ontology.term_name",
                  "assay_title"]})
    noms = {}
    for e in d.get("@graph", []):
        t = (e.get("biosample_ontology") or {}).get("term_name")
        if t:
            noms[t] = noms.get(t, 0) + 1
    return dict(sorted(noms.items(), key=lambda kv: -kv[1]))


def chercher_encode(biosample, assemblage="GRCh38", cible="CTCF"):
    """URL du meilleur fichier de pics, en desserrant les filtres.

    Une requete unique est trop fragile : le libelle exact d'un type
    cellulaire change d'un jeu a l'autre, et ENCODE ne distingue pas
    « filtre trop strict » de « donnee absente ». On essaie donc du plus
    precis au plus large, et si rien ne sort on affiche ce qui existe
    plutot que de laisser l'utilisateur deviner.
    """
    commun = dict(type="Experiment", status="released", format="json",
                  limit=25, **{"target.label": cible})
    essais = [
        ("terme exact, TF ChIP-seq",
         dict(commun, **{"assay_title": "TF ChIP-seq",
                         "biosample_ontology.term_name": biosample})),
        ("terme exact, tout essai",
         dict(commun, **{"biosample_ontology.term_name": biosample})),
        ("recherche libre",
         dict(commun, searchTerm=biosample)),
    ]
    vus = []
    for libelle, params in essais:
        d, url = _api("/search/", **params)
        g = d.get("@graph", [])
        if not g:
            print(f"    {libelle} : aucun résultat")
            continue
        print(f"    {libelle} : {len(g)} expérience(s)")
        for e in g:
            acc = e["accession"]
            terme = (e.get("biosample_ontology") or {}).get("term_name", "?")
            vus.append((acc, terme))
            h = _fichier_pics(acc, assemblage)
            if h:
                print(f"    retenu {acc} — « {terme} »")
                return acc, h
        print(f"    aucune de ces expériences n'a de pics {assemblage}")

    print(f"\n  Rien pour « {biosample} ».")
    if vus:
        print("  Expériences vues, mais sans fichier de pics "
              f"{assemblage} :")
        for acc, t in vus[:10]:
            print(f"    {acc}  « {t} »")
    print(f"\n  Types cellulaires pour lesquels ENCODE a un ChIP-seq "
          f"{cible} :")
    try:
        for t, k in list(biosamples_disponibles(cible).items())[:30]:
            print(f"    {k:3d} expérience(s)  « {t} »")
    except Exception as e:
        print(f"    (liste indisponible : {e})")
    raise SystemExit(
        f"\naucun fichier de pics pour « {biosample} ».\n"
        "Reprendre un libellé de la liste ci-dessus dans config.py, champ\n"
        "`encode`, ou fournir un BED local avec --pics.")


def lire_fenetre(chemin, chrom, debut, fin):
    """Sequence de [debut, fin) dans un FASTA, en majuscules.

    Remplace pyfaidx, qui n'est pas sur conda-forge et obligeait a ajouter
    un canal pour une seule fonction. Le fichier est parcouru une fois et
    seule la fenetre demandee est gardee en memoire, soit quelques Mo : il
    n'y a donc ni index a construire ni chromosome entier a charger.
    """
    debut, fin = max(0, int(debut)), int(fin)
    morceaux, pos, dedans = [], 0, False
    with open(chemin, "rb") as f:
        for ligne in f:
            if ligne[:1] == b">":
                if dedans:
                    break
                # l'en-tete peut porter une description apres le nom
                dedans = ligne[1:].split()[0].decode() == chrom
                pos = 0
                continue
            if not dedans:
                continue
            s = ligne.strip()
            if not s:
                continue
            d, g = pos, pos + len(s)
            pos = g
            if g <= debut:
                continue
            if d >= fin:
                break
            morceaux.append(s[max(0, debut - d):max(0, fin - d)])
    if not dedans and not morceaux:
        raise KeyError(f"{chrom} absent de {chemin}")
    return b"".join(morceaux).upper().decode()


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--biosample",
                   help="nom ENCODE du type cellulaire, par exemple "
                        "'foreskin melanocyte' ou 'fibroblast of dermis'")
    p.add_argument("--lister", action="store_true",
                   help="affiche les types cellulaires pour lesquels "
                        "ENCODE a un ChIP-seq CTCF, puis s'arrête")
    p.add_argument("--genome", help="FASTA du génome")
    p.add_argument("--chrom")
    p.add_argument("--debut", type=int)
    p.add_argument("--fin", type=int)
    p.add_argument("--pics", help="BED local, au lieu d'interroger ENCODE")
    p.add_argument("--marge", type=int, default=150,
                   help="demi-largeur balayée autour du sommet du pic")
    p.add_argument("--seuil", type=float, default=8.0,
                   help="score minimal du motif, en log2 de rapport de "
                        "vraisemblance")
    p.add_argument("--sortie")
    a = p.parse_args()

    if a.lister:
        print("Types cellulaires avec un ChIP-seq CTCF dans ENCODE :")
        for t, k in biosamples_disponibles().items():
            print(f"  {k:3d} expérience(s)  « {t} »")
        return

    manquant = [n for n, v in (("--genome", a.genome), ("--chrom", a.chrom),
                               ("--debut", a.debut), ("--fin", a.fin),
                               ("--sortie", a.sortie)) if v is None]
    if manquant:
        p.error("arguments requis : " + ", ".join(manquant))
    if not a.biosample and not a.pics:
        p.error("il faut --biosample ou --pics")

    import pandas as pd
    if a.pics:
        src = a.pics
        brut = Path(a.pics).read_bytes()
    else:
        acc, url = chercher_encode(a.biosample)
        print(f"  ENCODE {acc} : {url}")
        src = url
        with urllib.request.urlopen(url, timeout=300) as r:
            brut = r.read()
    if brut[:2] == b"\x1f\x8b":
        brut = gzip.decompress(brut)
    d = pd.read_csv(io.BytesIO(brut), sep="\t", header=None, comment="#")
    d = d[(d[0] == a.chrom) & (d[2] > a.debut) & (d[1] < a.fin)]
    print(f"  {len(d)} pics dans {a.chrom}:{a.debut:,}-{a.fin:,}")
    if d.empty:
        sys.exit("aucun pic dans la fenêtre")

    # la fenetre entiere est lue une fois, les pics y sont decoupes ensuite
    f0, f1 = max(0, a.debut - a.marge - 50), a.fin + a.marge + 50
    fenetre = lire_fenetre(a.genome, a.chrom, f0, f1)
    print(f"  séquence lue : {len(fenetre):,} pb depuis {f0:,}")
    if set(fenetre) <= {"N"}:
        sys.exit(f"la fenêtre ne contient que des N dans {a.genome} : "
                 "assemblage ou nom de chromosome incohérent")
    pwm = pwm_log_odds(MA0139)

    lignes = []
    for r in d.itertuples(index=False):
        s, e = int(r[1]), int(r[2])
        # colonne 9 du narrowPeak : sommet relatif ; sinon le milieu
        sommet = s + int(r[9]) if len(r) > 9 and r[9] != -1 else (s + e) // 2
        g, dr = max(sommet - a.marge, 0), sommet + a.marge
        seq = fenetre[g - f0:dr - f0]
        if len(seq) < len(MA0139[0]):
            continue
        score, brin, pos = balayer(seq, pwm)
        if score < a.seuil or brin == ".":
            continue
        c = g + pos + len(pwm) // 2
        signal = float(r[6]) if len(r) > 6 else float(r[4])
        lignes.append((a.chrom, c - 10, c + 10, "CTCF",
                       round(score, 2), brin, signal))

    out = pd.DataFrame(lignes, columns=["chrom", "debut", "fin", "nom",
                                        "score", "brin", "signal"])
    out.to_csv(a.sortie, sep="\t", header=False, index=False)
    n_p = (out.brin == "+").sum()
    print(f"  {len(out)} sites orientés retenus  ({n_p} en +, "
          f"{len(out) - n_p} en −)  →  {a.sortie}")
    print(f"  source : {src}")


if __name__ == "__main__":
    main()
