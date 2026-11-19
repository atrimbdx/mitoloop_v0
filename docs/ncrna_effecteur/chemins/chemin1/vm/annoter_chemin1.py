#!/usr/bin/env python3
"""Chemin 1, annotation sur la VM hors ligne : un variant d'ARN non codant nucléaire agit-il sur des gènes
mitochondriaux, via ses cibles (miARN, snoRNA, lncRNA, partenaires Pol III) ou via l'épissage (snRNA) ?

Entrée : VCF (.vcf / .vcf.gz) ou TSV (colonnes chr, pos, ref, alt), hg38.
Sortie : le TSV d'entrée (colonnes chr, pos, ref, alt + celles du TSV) avec une colonne NC_chemin1 =
         score|niveau|justification ; niveau dans fort / modere / faible / hypothese / aucun.

Aucun accès réseau. Les tables viennent du paquet préparé sur Curta (curta/prep_chemin1.py).
Les seuils sont fixés a priori dans SEUILS (justification dans README.md) ; ils ne doivent pas être
ajustés sur des variants de patients ni sur des sorties AlphaGenome / SpliceAI / SpliceAI2 / AlphaFold 3.
"""
import argparse
import csv
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from commun import (adn, appliquer_variant, arn, hypergeom_sf, lire_variants, meilleur, ouvrir, revcomp,
                    sortie)

SEUILS = {
    "p_fort": 1e-3,          # enrichissement mito des gènes touchés : niveau fort (miARN) / modere (snRNA)
    "k_fort": 5,             # nombre minimal de gènes mito touchés pour ce niveau
    "p_modere": 0.05,
    "k_modere": 3,
    "ddg_precurseur": 2.0,   # |ΔΔG| (kcal/mol) du repliement du précurseur de miARN
    "dist_clivage": 2,       # nt d'un site de clivage Drosha / Dicer (extrémités des bras matures)
    "ddg_appariement": 1.5,  # |ΔΔG| (kcal/mol) de l'appariement snRNA / site d'épissage
    "r_coexpression": 0.5,   # |r| de co-expression GTEx lncRNA / module OXPHOS
}
REGIONS_SNRNA = {"U1": ("ATACTTACCTG", 0, 0), "U6": ("ACAGAGA", 2, 2), "U2": ("GTAGTA", 3, 3)}


def lire_tsv(path):
    with ouvrir(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fmt_p(p):
    return "%.2g" % p


def score_p(p):
    return round(min(50.0, -math.log10(max(p, 1e-300))), 2)


class Paquet:
    def __init__(self, d, mito_path=None, snodb=None, boucles=None, coexpr=None):
        self.d = d
        mp = mito_path or os.path.join(d, "genes_mito.tsv")
        with ouvrir(mp) as fh:
            first = fh.readline().rstrip("\n").split("\t")
            rest = [l.rstrip("\n").split("\t")[0] for l in fh]
        self.mito = set(([] if first[0] == "gene" else [first[0]]) + [g for g in rest if g])
        self.elements = []
        for nom, cl in (("mirna.tsv.gz", "miARN_precurseur"), ("elements_ncrna.tsv.gz", None)):
            p = os.path.join(d, nom)
            if os.path.exists(p):
                for r in lire_tsv(p):
                    r["_type"] = cl or r["classe"]
                    if r["_type"] == "miARN":  # déjà couvert, avec bras matures, par mirna.tsv.gz
                        continue
                    self.elements.append(r)
        self.index = {}
        for e in self.elements:
            self.index.setdefault(e["chrom"], []).append(e)
        self._utr = self._sites = None
        self._u12 = self._part = None
        self.snodb = lire_tsv(snodb) if snodb else None
        self.boucles = lire_tsv(boucles) if boucles else None
        self.coexpr = lire_tsv(coexpr) if coexpr else None
        self._cache = {}

    def chevauche(self, chrom, a, b):
        return [e for e in self.index.get(chrom, []) if int(e["debut0"]) < b and a < int(e["fin"])]

    @property
    def utr(self):
        if self._utr is None:
            self._utr = {r["gene"]: r["sequence"] for r in lire_tsv(os.path.join(self.d, "utr3.tsv.gz"))}
        return self._utr

    @property
    def sites(self):
        if self._sites is None:
            self._sites = lire_tsv(os.path.join(self.d, "sites_epissage.tsv.gz"))
        return self._sites

    @property
    def u12(self):
        if self._u12 is None:
            self._u12 = []
            for nom in ("gene_mito_introns_U12.tsv", "gene_mito_introns_U12_intronIC.tsv"):
                p = os.path.join(self.d, nom)
                if os.path.exists(p):
                    self._u12 += lire_tsv(p)
            pg = os.path.join(self.d, "resume_U12_genome_entier.tsv")
            self.u12_resume = lire_tsv(pg)[0] if os.path.exists(pg) else None
        return self._u12

    @property
    def partenaires(self):
        if self._part is None:
            p = os.path.join(self.d, "partenaires_pol3_arnt_arnr.tsv")
            self._part = lire_tsv(p) if os.path.exists(p) else []
        return self._part


# ------------------------------------------------------------- outils variant
def transcrit_ref_alt(e, v):
    """Séquences REF / ALT de l'élément en orientation du transcrit, offset 0-based du variant
    dans le transcrit, concordance de l'allèle REF."""
    g0, g1 = int(e["debut0"]), int(e["fin"])
    seq_t = e.get("sequence_transcrit_flancs") or e.get("sequence_brin_plus", ".")
    if seq_t in (".", ""):
        return None
    plus = seq_t if ("sequence_brin_plus" in e or e["brin"] == "+") else revcomp(seq_t)
    alt_plus, ok = appliquer_variant(plus, g0, v["pos"], v["ref"], v["alt"])
    if alt_plus is None:
        return None
    if e["brin"] == "+":
        return plus, alt_plus, v["pos"] - 1 - g0, ok
    return revcomp(plus), revcomp(alt_plus), (g1 - 1) - (v["pos"] - 1) - (len(v["ref"]) - 1), ok


def mfe(seq):
    import RNA
    return RNA.fold(arn(seq))[1]


def duplex(a, b):
    import RNA
    return RNA.duplexfold(arn(a), arn(b)).energy


# ------------------------------------------------------------------- miARN
def cibles(paquet, mature):
    """Gènes portant au moins un site 8mer, 7mer-m8 ou 7mer-A1 (Bartel 2009 ; Agarwal 2015) dans l'UTR 3'."""
    m = adn(mature)
    if len(m) < 8:
        return set()
    k7m8, k7a1 = revcomp(m[1:8]), revcomp(m[1:7]) + "A"
    cle = (k7m8, k7a1)
    if cle not in paquet._cache:
        paquet._cache[cle] = {g for g, s in paquet.utr.items() if k7m8 in s or k7a1 in s}
    return paquet._cache[cle]


def enrichissement(paquet, genes):
    univers = set(paquet.utr)
    N, K = len(univers), len(paquet.mito & univers)
    gm = genes & paquet.mito
    return len(gm), len(genes), hypergeom_sf(len(gm), N, K, len(genes)), sorted(gm), K / N


def ajuste_longueur(paquet, genes):
    """Information seulement (pas de niveau) : nombre attendu de gènes mito dans `genes` en stratifiant
    l'univers par déciles de longueur d'UTR 3' (les UTR des gènes mito sont plus courts), et z-score."""
    if not hasattr(paquet, "_deciles"):
        lg = sorted((len(s), g) for g, s in paquet.utr.items())
        n = len(lg)
        paquet._deciles = {g: min(9, i * 10 // n) for i, (_, g) in enumerate(lg)}
        paquet._strates = [[0, 0] for _ in range(10)]
        for g, d in paquet._deciles.items():
            paquet._strates[d][0] += 1
            paquet._strates[d][1] += g in paquet.mito
    nd = [0] * 10
    for g in genes:
        if g in paquet._deciles:
            nd[paquet._deciles[g]] += 1
    E = V = 0.0
    for d in range(10):
        Nd, Kd = paquet._strates[d]
        if Nd > 1 and nd[d]:
            q = Kd / Nd
            E += nd[d] * q
            V += nd[d] * q * (1 - q) * (Nd - nd[d]) / (Nd - 1)
    k = len(genes & paquet.mito)
    return "obs/att ajusté longueur UTR %d/%.1f (z=%+.1f)" % (k, E, (k - E) / math.sqrt(V) if V > 0 else 0.0)


def lister(g, n=12):
    g = sorted(g)
    return ",".join(g[:n]) + (",…(+%d)" % (len(g) - n) if len(g) > n else "") if g else "aucun"


def bras(e):
    out = []
    if e["matures_nom:offset0:seq"] != ".":
        for x in e["matures_nom:offset0:seq"].split(";"):
            n, o, s = x.split(":")
            out.append((int(o), n, s))
    return sorted(out)


def annoter_mirna(paquet, e, v):
    t = transcrit_ref_alt(e, v)
    if t is None:
        return None
    ref_t, alt_t, off, ok = t
    nom = e["gene"]
    if not ok:
        return (0, "aucun", "%s : allèle REF discordant avec le génome de référence, non évalué" % nom)
    arms = bras(e)
    if not arms:
        return (0, "hypothese", "%s : précurseur sans bras mature miRBase localisé" % nom)
    L = max(len(v["ref"]), 1)
    res = []
    for o, n, m in arms:
        if off + L > o and off < o + len(m):
            pos_m = off - o + 1
            if 2 <= pos_m <= 8 or (pos_m < 2 and off + L - 1 - o + 1 >= 2):
                m_alt = alt_t[o:o + len(m)]
                ref_c, alt_c = cibles(paquet, m), cibles(paquet, m_alt)
                perdus, gagnes = ref_c - alt_c, alt_c - ref_c
                kp, np_, pp, gp, fond = enrichissement(paquet, perdus)
                kg, ng, pg, gg, _ = enrichissement(paquet, gagnes)
                kr, nr, pr, _, _ = enrichissement(paquet, ref_c)
                p, k = (pp, kp) if pp <= pg else (pg, kg)
                if p < SEUILS["p_fort"] and k >= SEUILS["k_fort"]:
                    niv = "fort"
                elif p < SEUILS["p_modere"] and k >= SEUILS["k_modere"]:
                    niv = "modere"
                elif kp + kg >= 1:
                    niv = "hypothese"  # cibles mito changées mais sans enrichissement
                else:
                    niv = "aucun"
                res.append((score_p(p), niv,
                            "%s %s graine pos %d (%s>%s) : cibles REF %d dont %d mito (P=%s) ; perdues %d dont "
                            "%d mito (P=%s ; %s : %s) ; gagnées %d dont %d mito (P=%s ; %s : %s) ; fond mito "
                            "%.1f %% ; REF : %s"
                            % (nom, n, pos_m, m[1:8], m_alt[1:8], nr, kr, fmt_p(pr), np_, kp, fmt_p(pp),
                               ajuste_longueur(paquet, perdus), lister(gp), ng, kg, fmt_p(pg),
                               ajuste_longueur(paquet, gagnes), lister(gg), 100 * fond,
                               ajuste_longueur(paquet, ref_c))))
    if res:
        return sorted(res, key=lambda r: -r[0])[0]
    # hors graine : repliement du précurseur et sites de clivage
    ddg = mfe(alt_t) - mfe(ref_t)
    coupes = []
    for o, n, m in arms:
        coupes += [o, o + len(m)]
    dist = min(min(abs(off - c), abs(off + L - c)) for c in coupes)
    dans_mature = any(off + L > o and off < o + len(m) for o, n, m in arms)
    if abs(ddg) >= SEUILS["ddg_precurseur"] or dist <= SEUILS["dist_clivage"]:
        tout = set()
        for o, n, m in arms:
            tout |= cibles(paquet, m)
        k, n_, p, gm, fond = enrichissement(paquet, tout)
        niv = "faible" if (p < SEUILS["p_modere"] and k >= SEUILS["k_modere"]) else "hypothese"
        return (score_p(p), niv,
                "%s : %s, ΔΔG précurseur %+.1f kcal/mol, %d nt d'un site Drosha/Dicer ; si perte de maturation : "
                "%d cibles dont %d mito (P=%s ; fond %.1f %% ; %s)"
                % (nom, "mature hors graine" if dans_mature else "précurseur hors bras matures", ddg, dist, n_, k,
                   fmt_p(p), 100 * fond, lister(gm)))
    return (0, "aucun", "%s : %s, ΔΔG %+.1f kcal/mol, %d nt des sites de clivage (sous les seuils)"
            % (nom, "mature hors graine" if dans_mature else "précurseur", ddg, dist))


# ---------------------------------------------------------------- épissage
def annoter_snrna_majeur(paquet, e, v, sc):
    t = transcrit_ref_alt(e, v)
    if t is None:
        return None
    ref_t, alt_t, off, ok = t
    nom = "%s (%s)" % (e["gene"], sc)
    if not ok:
        return (0, "aucun", "%s : allèle REF discordant, non évalué" % nom)
    if sc not in REGIONS_SNRNA:
        return (0, "hypothese", "%s : snRNA du spliceosome majeur hors région d'appariement aux sites "
                "(effet global possible, non spécifique des gènes mito)" % nom)
    motif, avant, apres = REGIONS_SNRNA[sc]
    i = ref_t.find(motif)
    if i < 0:
        return (0, "hypothese", "%s : motif %s absent du paralogue, appariement non évalué" % (nom, motif))
    a, b = max(0, i - avant), i + len(motif) + apres
    L = max(len(v["ref"]), 1)
    if not (off + L > a and off < b):
        return (0, "hypothese", "%s : variant hors région d'appariement (%d-%d), effet global possible"
                % (nom, a + 1, b))
    reg_ref, reg_alt = ref_t[a:b], alt_t[a:b + len(v["alt"]) - len(v["ref"])]
    col, cible_desc = {"U1": ("ss5_m3_p6", "site 5' (-3..+6)"), "U6": ("intron_p1_p12", "intron +1..+12"),
                       "U2": ("intron_fin50", "50 nt en 3' de l'intron (point de branchement)")}[sc]
    cache = {}
    touches, total = set(), set()
    for s in paquet.sites:
        q = s[col] if sc != "U2" else s[col][:42]
        total.add(s["gene"])
        if q not in cache:
            cache[q] = duplex(reg_alt, q) - duplex(reg_ref, q)
        if abs(cache[q]) >= SEUILS["ddg_appariement"]:
            touches.add(s["gene"])
    N, K = len(total), len(paquet.mito & total)
    gm = touches & paquet.mito
    p = hypergeom_sf(len(gm), N, K, len(touches))
    if p < SEUILS["p_fort"] and len(gm) >= SEUILS["k_fort"]:
        niv = "modere"
    elif gm:
        niv = "hypothese"  # gènes mito touchés sans enrichissement : effet global attendu
    else:
        niv = "aucun"
    return (score_p(p), niv,
            "%s : région d'appariement %d-%d (%s>%s) contre %s : %d gènes avec |ΔΔG|>=%.1f, dont %d mito "
            "(P=%s ; fond %.1f %% ; %s) ; effet attendu global sur l'épissage"
            % (nom, a + 1, b, arn(reg_ref), arn(reg_alt), cible_desc, len(touches), SEUILS["ddg_appariement"],
               len(gm), fmt_p(p), 100 * K / max(N, 1), lister(gm)))


def annoter_snrna_mineur(paquet, e, v, sc):
    u12 = paquet.u12
    gm = sorted(r["gene_hgnc"] for r in u12 if r.get("dans_ensemble_mito", "oui") != "non")
    gm_set = set(gm) & paquet.mito if paquet.mito else set(gm)
    txt = "%s (%s, spliceosome mineur) : %d gènes mito à intron U12 (%s)" % (e["gene"], sc, len(gm_set),
                                                                          lister(gm_set, 20))
    niv = "hypothese" if gm_set else "aucun"
    if getattr(paquet, "u12_resume", None):
        r = paquet.u12_resume
        p = float(r["P_hypergeom_enrichissement"])
        txt += " ; génome entier : %s/%s gènes U12 sont mito (fond %.1f %%, P enrichissement=%s)" % (
            r["genes_U12_mito"], r["genes_U12"], 100 * float(r["part_mito_fond"]), fmt_p(p))
        if p < SEUILS["p_modere"] and gm_set:
            niv = "faible"
    else:
        txt += " ; part génome entier non disponible (sortie intronIC absente du paquet)"
    return (len(gm_set), niv, txt + " ; effet attendu global, non spécifique des gènes mito")


# ---------------------------------------------------- snoRNA, lncRNA, Pol III
def annoter_snorna(paquet, e, v):
    if paquet.snodb is None:
        return (0, "aucun", "%s (snoRNA) : cibles snoDB non chargées (--snodb), non évalué" % e["gene"])
    cib = {r["cible_gene"] for r in paquet.snodb if r["snorna"] == e["gene"]}
    gm = cib & paquet.mito
    if gm:
        return (len(gm), "faible", "%s (snoRNA) : cibles connues (snoDB) dans des gènes mito : %s ; position "
                "du variant dans le guide non évaluée (v0)" % (e["gene"], lister(gm)))
    return (0, "aucun", "%s (snoRNA) : %d cibles snoDB, aucune mito" % (e["gene"], len(cib)))


def annoter_lncrna(paquet, e, v):
    res = []
    anti = {g for g in e.get("antisens_chevauchant_codant", ".").split(",") if g and g != "."}
    gm = anti & paquet.mito
    if gm:
        res.append((len(gm), "faible", "%s (lncRNA) : antisens chevauchant le(s) gène(s) mito %s (cis)"
                    % (e["gene"], lister(gm))))
    if paquet.boucles is not None:
        lb = {r["gene_mito"] for r in paquet.boucles if r["chrom"] == v["chrom"]
              and int(r["debut0"]) < v["pos"] <= int(r["fin"])} & paquet.mito
        if lb:
            res.append((len(lb), "hypothese", "%s (lncRNA) : variant dans une ancre de boucle mitoloop reliée à %s"
                        % (e["gene"], lister(lb))))
    if paquet.coexpr is not None:
        co = [r for r in paquet.coexpr if r["lncrna"] == e["gene"] and abs(float(r["r"])) >= SEUILS["r_coexpression"]]
        if co:
            res.append((round(max(abs(float(r["r"])) for r in co), 2), "hypothese",
                        "%s (lncRNA) : co-expression GTEx |r|>=%.1f avec %s" % (
                            e["gene"], SEUILS["r_coexpression"],
                            ", ".join("%s/%s r=%s" % (r["module"], r["tissu"], r["r"]) for r in co[:5]))))
    if not res:
        return (0, "aucun", "%s (lncRNA) : aucun lien cis/trans avec un gène mito dans les tables chargées"
                % e["gene"])
    return sorted(res, key=lambda r: ["fort", "modere", "faible", "hypothese", "aucun"].index(r[1]))[0]


def annoter_pol3(paquet, e, v):
    g = e["gene"]
    cle = "ARNt" if e["classe"] == "ARNt" else next(
        (k for k in ("RMRP", "RPPH1", "TERC", "RN7SK", "RN7SL", "RNA5S") if g.startswith(k)), None)
    if cle is None:
        return (0, "aucun", "%s (%s) : pas de partenaire protéique mito répertorié" % (g, e["sous_classe"]))
    parts = [r for r in paquet.partenaires if r["arn"] == cle and r["partenaire"] in paquet.mito
             and not r["agit_sur_arn_nucleaire"].startswith("non")]
    if not parts:
        return (0, "aucun", "%s : partenaires protéiques connus (%s) sans gène de l'ensemble mito" % (g, cle))
    return (len(parts), "hypothese", "%s : partenaires protéiques de l'ensemble mito : %s"
            % (g, "; ".join("%s (%s)" % (r["partenaire"], r["lien"]) for r in parts)))


def annoter(paquet, v):
    a = v["pos"] - 1
    b = a + max(len(v["ref"]), 1)
    res = []
    for e in paquet.chevauche(v["chrom"], a, b):
        t, sc = e["_type"], e.get("sous_classe", "")
        if t == "miARN_precurseur":
            r = annoter_mirna(paquet, e, v)
        elif t == "snRNA" and sc in ("U11", "U12", "U4atac", "U6atac"):
            r = annoter_snrna_mineur(paquet, e, v, sc)
        elif t == "snRNA":
            r = annoter_snrna_majeur(paquet, e, v, sc)
        elif t == "snoRNA":
            r = annoter_snorna(paquet, e, v)
        elif sc.startswith("PolIII_") or t in ("ARNt", "ARNr", "ARN_voute", "ARN_Y", "ARN_misc"):
            r = annoter_pol3(paquet, e, v)
        elif t == "lncRNA":
            r = annoter_lncrna(paquet, e, v)
        else:
            r = None
        if r:
            res.append(r)
    return meilleur(res)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", required=True, help="VCF ou TSV chr/pos/ref/alt (hg38)")
    ap.add_argument("--paquet", required=True, help="dossier des tables préparées sur Curta")
    ap.add_argument("--mito", help="liste de gènes mito (remplace genes_mito.tsv du paquet)")
    ap.add_argument("--snodb", help="TSV snorna/cible_gene (snoDB, préparé sur Curta)")
    ap.add_argument("--boucles", help="TSV gene_mito/chrom/debut0/fin (ancres de boucles mitoloop)")
    ap.add_argument("--coexpression", help="TSV lncrna/module/tissu/r (GTEx)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    paquet = Paquet(a.paquet, a.mito, a.snodb, a.boucles, a.coexpression)
    with open(a.out, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        entete = False
        for v in lire_variants(a.variants):
            if not entete:
                w.writerow(["chr", "pos", "ref", "alt"] + v["header"] + ["NC_chemin1"])
                entete = True
            w.writerow([v["chrom"], v["pos"], v["ref"], v["alt"]] + v["ligne"] + [annoter(paquet, v)])


if __name__ == "__main__":
    main()
