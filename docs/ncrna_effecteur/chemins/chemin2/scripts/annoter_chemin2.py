#!/usr/bin/env python3
"""Chemin 2 - annotation des variants sur la VM (hors ligne).

Hypothese testee : l'ARN non codant nucleaire est present ou importe dans la mitochondrie
et le variant touche un element necessaire (import, liaison, repliement).

Entree  : TSV avec colonnes chr, pos, ref, alt (en-tete obligatoire, insensible a la casse)
          ou VCF (ligne #CHROM ; un ALT par ligne recommande, les ALT multiples sont separes).
          Coordonnees hg38, positions 1-based, alleles sur le brin + (convention VCF).
Sortie  : l'entree + colonne NC_chemin2 = score|niveau|justification
          niveau dans {aucun, hypothese, faible} : JAMAIS plus haut, JAMAIS de points ACMG.

Donnees (dossier ../donnees par defaut, preparees sur une machine connectee) :
  localisation_mito_ncrna.tsv (dans le dossier parent), elements_mito.tsv,
  arn_coordonnees_hg38.tsv, ref_regions_hg38.fa
Dependance : ViennaRNA (pip install ViennaRNA), testee avec 2.8.0.

Aucun acces reseau. Aucun seuil appris : toutes les constantes sont fixees a priori ci-dessous.
"""
import argparse
import math
import os
import re
import sys

try:
    import RNA  # ViennaRNA
except ImportError:  # le script reste utilisable sans pliage
    RNA = None

# ---------------------------------------------------------------------------
# Constantes fixees A PRIORI (non ajustees sur des variants, ni sur AlphaGenome,
# SpliceAI, SpliceAI2 ou AlphaFold 3).
# ---------------------------------------------------------------------------
# Points par niveau de preuve de localisation (rang 1 = APEX matrice ... 4 = fractionnement
# seul ; 5 = import in vitro, aptitude seulement). Echelle ordinale de chemins_mito.md.
POINTS_RANG = {1: 3.0, 2: 2.0, 3: 1.0, 4: 0.5, 5: 0.25}
# Poids d'un element touche : 2 = plage publiee, 1 = domaine publie mais coordonnees derivees.
# (lus dans elements_mito.tsv, colonne poids)
# Pliage : |ddG| >= 2.0 kcal/mol ~ perte d'un empilement de paires GC/CG (Turner 2004 :
# -2,4 a -3,4 kcal/mol) ; distance de paires >= max(5, 10 % de la longueur repliee).
SEUIL_DDG = 2.0
SEUIL_BPDIST_FRAC = 0.10
SEUIL_BPDIST_MIN = 5
# Repliement : ARN entier si <= 400 nt, sinon fenetre de +/- 75 nt autour du variant
# (sequence genomique du brin du transcrit ; l'epissage n'est pas reconstruit).
LONG_MAX_ENTIER = 400
DEMI_FENETRE = 75
# Familles multicopies : une copie variante parmi >= 17 -> effet dilue, appel peu fiable.
MULTICOPIE = re.compile(r"^RNA5S\d+$")
# Niveau "faible" (plafond du chemin 2) seulement si : preuve de rang <= 2 chez l'homme dans
# >= 2 publications independantes ET variant dans un element a coordonnees publiees ET non
# multicopie. Sinon "hypothese".
MIN_ETUDES_FAIBLE = 2
RANG_MAX_FAIBLE = 2

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s):
    return s.translate(COMP)[::-1]


def ici(*p):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), *p)


def lire_tsv(chemin):
    with open(chemin) as fh:
        tete = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(tete, l.rstrip("\n").split("\t"))) for l in fh if l.strip()]


def developper_symboles(txt):
    """'RMRP ; RNA5S1-RNA5S17' -> {'RMRP', 'RNA5S1', ..., 'RNA5S17'}"""
    out = set()
    for tok in re.split(r"[;,]", txt):
        tok = tok.strip()
        m = re.fullmatch(r"([A-Z0-9]+?)(\d+)-([A-Z0-9]+?)(\d+)", tok)
        if m and m.group(1) == m.group(3):
            for i in range(int(m.group(2)), int(m.group(4)) + 1):
                out.add(f"{m.group(1)}{i}")
        elif re.fullmatch(r"[A-Z0-9][A-Z0-9-]*", tok or "-") and tok not in ("-", "—"):
            out.add(tok)
    return out


class Donnees:
    def __init__(self, ddir, table):
        self.preuves = {}
        for r in lire_tsv(table):
            for g in developper_symboles(r["genes_hgnc"]):
                self.preuves.setdefault(g, []).append(r)
        self.apex_negatif = [r for r in lire_tsv(table) if r["type_preuve"] == "apex_matrice_negatif"]
        self.elements = {}
        for r in lire_tsv(os.path.join(ddir, "elements_mito.tsv")):
            self.elements.setdefault(r["gene"], []).append(r)
        # un enregistrement par gene : le plus long (evite les doublons GENCODE ribozyme/misc_RNA)
        self.genes = {}
        for r in lire_tsv(os.path.join(ddir, "arn_coordonnees_hg38.tsv")):
            g = r["gene"]
            if g not in self.preuves:
                continue
            L = int(r["fin"]) - int(r["debut"])
            if g not in self.genes or L > int(self.genes[g]["fin"]) - int(self.genes[g]["debut"]):
                self.genes[g] = r
        self.regions = []  # (chr, debut, fin, seq)
        nom, buf = None, []
        with open(os.path.join(ddir, "ref_regions_hg38.fa")) as fh:
            for l in fh:
                if l.startswith(";"):
                    continue
                if l.startswith(">"):
                    if nom:
                        self.regions.append(nom + ("".join(buf),))
                    c, se = l[1:].split()[0].split(":")
                    s, e = se.split("-")
                    nom, buf = (c, int(s), int(e)), []
                else:
                    buf.append(l.strip().upper())
            if nom:
                self.regions.append(nom + ("".join(buf),))

    def seq(self, c, s, e):
        for rc_, rs, re_, sq in self.regions:
            if rc_ == c and rs <= s and e <= re_:
                return sq[s - rs:e - rs + 1]
        return None


def exons(r):
    if r.get("exons_union"):
        return [tuple(map(int, x.split("-"))) for x in r["exons_union"].split(",")]
    return [(int(r["debut"]), int(r["fin"]))]


def chevauche(a, b, c, d):
    return a <= d and c <= b


def resume_preuves(rows):
    pos = [r for r in rows if r["type_preuve"] not in ("contre_preuve", "apex_matrice_negatif")]
    hum = [r for r in pos if "humain" in r["espece"]]
    rang_h = min([int(r["rang_preuve"]) for r in hum if r["rang_preuve"]] or [9])
    rang_all = min([int(r["rang_preuve"]) for r in pos if r["rang_preuve"]] or [9])
    dois = {r["doi"] for r in hum if r["rang_preuve"] and int(r["rang_preuve"]) <= RANG_MAX_FAIBLE}
    contre = [r["reference"].split(",")[0] for r in rows if r["type_preuve"] == "contre_preuve"]
    return rang_h, rang_all, len(dois), contre


def plier(seq):
    s = seq.replace("T", "U")
    st, e = RNA.fold(s)
    return st, e


def analyse_pliage(d, g, r, c, pos, ref, alt):
    """Retourne (ddg, bpdist, seuil_bp, apparie, note) ; None si non calculable."""
    if RNA is None:
        return None, None, None, None, "ViennaRNA absent : pliage non calcule"
    gs, ge, brin = int(r["debut"]), int(r["fin"]), r["brin"]
    vend = pos + len(ref) - 1
    if ge - gs + 1 <= LONG_MAX_ENTIER:
        ws, we, mode = gs, ge, "ARN entier"
    else:
        ws, we, mode = max(gs, pos - DEMI_FENETRE), min(ge, vend + DEMI_FENETRE), f"fenetre +/-{DEMI_FENETRE} nt (epissage ignore)"
    if pos < ws or vend > we:
        return None, None, None, None, "variant en bordure de l'ARN : pliage non calcule"
    ref_g = d.seq(c, ws, we)
    if ref_g is None:
        return None, None, None, None, "sequence de reference absente du paquet"
    off = pos - ws
    if ref_g[off:off + len(ref)] != ref.upper():
        return None, None, None, None, f"REF discordant avec la reference ({ref_g[off:off+len(ref)]} attendu)"
    alt_g = ref_g[:off] + alt.upper() + ref_g[off + len(ref):]
    if brin == "-":
        ref_t, alt_t = rc(ref_g), rc(alt_g)
        p_t = (we - vend) + 1  # position (1-based) du variant dans la fenetre, brin du transcrit
    else:
        ref_t, alt_t = ref_g, alt_g
        p_t = off + 1
    st_r, e_r = plier(ref_t)
    st_a, e_a = plier(alt_t)
    ddg = round(e_a - e_r, 2)
    apparie = st_r[p_t - 1] != "." if 0 < p_t <= len(st_r) else None
    if len(st_r) == len(st_a):
        bp = RNA.bp_distance(st_r, st_a)
    else:
        bp = None
    seuil_bp = max(SEUIL_BPDIST_MIN, math.ceil(SEUIL_BPDIST_FRAC * len(st_r)))
    return ddg, bp, seuil_bp, apparie, f"{mode}, {len(st_r)} nt"


def position_transcrit(r, pos):
    gs, ge = int(r["debut"]), int(r["fin"])
    return (ge - pos + 1) if r["brin"] == "-" else (pos - gs + 1)


def annoter(d, c, pos, ref, alt):
    vend = pos + max(len(ref), 1) - 1
    res = []
    for g, r in d.genes.items():
        if r["chr"] != c or not any(chevauche(pos, vend, a, b) for a, b in exons(r)):
            continue
        rang_h, rang_all, n_etudes, contre = resume_preuves(d.preuves[g])
        pts_p = POINTS_RANG.get(rang_all, 0.0)
        # elements
        hit, pts_e, el_txt = None, 0.0, []
        for el in d.elements.get(g, []):
            if el["type_element"] == "non_cartographie":
                el_txt.append(f"element non cartographie ({el['element']})")
                continue
            if chevauche(pos, vend, int(el["debut"]), int(el["fin"])):
                w = float(el["poids"])
                if w > pts_e:
                    pts_e, hit = w, el
        # pliage
        ddg, bp, seuil_bp, apparie, note_pl = analyse_pliage(d, g, r, c, pos, ref, alt)
        pts_f = 0.0
        if ddg is not None and abs(ddg) >= SEUIL_DDG:
            pts_f += 1
        if bp is not None and bp >= seuil_bp:
            pts_f += 1
        score = round(pts_p + pts_e + pts_f, 2)
        multi = bool(MULTICOPIE.match(g))
        ref_ok = not note_pl.startswith("REF discordant")
        faible = (rang_h <= RANG_MAX_FAIBLE and n_etudes >= MIN_ETUDES_FAIBLE and ref_ok
                  and hit is not None and hit["type_element"] == "publie_plage" and not multi)
        niveau = "faible" if faible else "hypothese"
        j = [f"{g} n.{position_transcrit(r, pos)} (GENCODE v27)",
             f"preuve mito: meilleur rang {rang_all} (homme {rang_h if rang_h < 9 else 'aucun'}), {n_etudes} etude(s) humaine(s) rang<=2"]
        if contre:
            j.append("conteste: " + ", ".join(sorted(set(contre))))
        if int(r["fin"]) - int(r["debut"]) + 1 > LONG_MAX_ENTIER:
            j.append("APEX matrice negatif pour les ARN nucleaires longs")
        if hit is not None:
            j.append(f"dans element {hit['type_element']}: {hit['element']}")
        else:
            j.append("hors element d'import cartographie" + (f" ({'; '.join(el_txt)})" if el_txt else ""))
        if ddg is not None:
            j.append(f"ddG={ddg:+.2f} kcal/mol (seuil {SEUIL_DDG})"
                     + (f", dist. paires={bp} (seuil {seuil_bp})" if bp is not None else ", indel: dist. paires non calculee")
                     + (", nt apparie dans la MFE REF" if apparie else ", nt non apparie dans la MFE REF" if apparie is not None else "")
                     + f" [{note_pl}]")
        else:
            j.append(note_pl)
        if multi:
            j.append("famille multicopie (>=17 copies): plafond hypothese")
        j.append("chemin d'hypothese, aucun critere ACMG")
        res.append((score, niveau, "; ".join(j).replace("|", "/")))
    if not res:
        return "0|aucun|hors table de localisation mitochondriale (chemin 2)"
    res.sort(key=lambda x: (x[1] == "faible", x[0]), reverse=True)
    best = res[0]
    just = " // ".join(x[2] for x in res)
    return f"{best[0]}|{best[1]}|{just}"


def lire_variants(chemin):
    with open(chemin) as fh:
        lignes = [l.rstrip("\n") for l in fh]
    tete, corps, vcf = None, [], False
    for l in lignes:
        if l.startswith("##"):
            continue
        if tete is None:
            tete = l.lstrip("#").split("\t")
            vcf = l.startswith("#CHROM")
            continue
        if l.strip():
            corps.append(l.split("\t"))
    if vcf:
        ic, ip, ir, ia = 0, 1, 3, 4
    else:
        low = [h.lower() for h in tete]
        try:
            ic, ip, ir, ia = (low.index(x) for x in ("chr", "pos", "ref", "alt"))
        except ValueError:
            sys.exit("en-tete attendu : chr, pos, ref, alt (TSV) ou VCF")
    return tete, corps, (ic, ip, ir, ia), vcf


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--variants", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--donnees", default=ici("..", "donnees"))
    p.add_argument("--table", default=ici("..", "localisation_mito_ncrna.tsv"))
    a = p.parse_args()
    d = Donnees(a.donnees, a.table)
    tete, corps, (ic, ip, ir, ia), vcf = lire_variants(a.variants)
    with open(a.out, "w") as out:
        out.write("\t".join(tete) + "\tNC_chemin2\n")
        for f in corps:
            c = f[ic] if f[ic].startswith("chr") else "chr" + f[ic]
            vals = []
            for alt in f[ia].split(","):
                vals.append(annoter(d, c, int(f[ip]), f[ir], alt))
            out.write("\t".join(f) + "\t" + " ;; ".join(vals) + "\n")


if __name__ == "__main__":
    main()
