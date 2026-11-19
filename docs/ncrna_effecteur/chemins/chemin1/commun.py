"""Fonctions communes du chemin 1 (préparation Curta et annotation VM).

Bibliothèque standard uniquement, sauf ViennaRNA (module RNA) importé à la demande.
Coordonnées : 0-based demi-ouvertes en interne (comme refFlat/BED), 1-based pour les variants (VCF).
"""
import csv
import gzip
import math
import os
import re
import sys

csv.field_size_limit(sys.maxsize)

NIVEAUX = ["fort", "modere", "faible", "hypothese", "aucun"]  # du plus fort au plus faible
_COMP = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def ouvrir(path, mode="rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def revcomp(s):
    return s.translate(_COMP)[::-1]


def adn(s):
    return s.upper().replace("U", "T")


def arn(s):
    return s.upper().replace("T", "U")


def norm_chrom(c):
    c = str(c)
    return c if c.startswith("chr") else "chr" + c


# --------------------------------------------------------------------- génome
class Fasta:
    """Accès aléatoire à un FASTA indexé (.fai), sans dépendance."""

    def __init__(self, path):
        self.f = open(path, "rb")
        self.idx = {}
        with open(path + ".fai") as fh:
            for line in fh:
                name, length, offset, lb, lw = line.split("\t")[:5]
                self.idx[name] = (int(length), int(offset), int(lb), int(lw))

    def seq(self, chrom, start, end):
        """Séquence [start, end) en 0-based, majuscules, brin +."""
        length, offset, lb, lw = self.idx[chrom]
        start, end = max(0, start), min(end, length)
        if end <= start:
            return ""
        a = offset + (start // lb) * lw + start % lb
        b = offset + (end // lb) * lw + end % lb
        self.f.seek(a)
        return self.f.read(b - a).replace(b"\n", b"").replace(b"\r", b"").decode().upper()


# ---------------------------------------------------------------- annotation
def lire_refflat(path):
    """refFlat UCSC : geneName name chrom strand txStart txEnd cdsStart cdsEnd exonCount exonStarts exonEnds."""
    vus = set()
    with ouvrir(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 11 or p[1] in vus:
                continue
            vus.add(p[1])
            starts = [int(x) for x in p[9].strip(",").split(",") if x]
            ends = [int(x) for x in p[10].strip(",").split(",") if x]
            yield {"gene": p[0], "tx": p[1], "chrom": p[2], "strand": p[3],
                   "start": int(p[4]), "end": int(p[5]), "cds_start": int(p[6]), "cds_end": int(p[7]),
                   "exons": list(zip(starts, ends))}


def utr3_blocs(tx):
    """Blocs génomiques de l'UTR 3' d'un transcrit codant (ordre génomique)."""
    if tx["cds_start"] >= tx["cds_end"]:
        return []
    out = []
    for s, e in tx["exons"]:
        if tx["strand"] == "+":
            s2 = max(s, tx["cds_end"])
            if s2 < e:
                out.append((s2, e))
        else:
            e2 = min(e, tx["cds_start"])
            if s < e2:
                out.append((s, e2))
    return out


def introns(tx):
    ex = sorted(tx["exons"])
    return [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]


# ---------------------------------------------------------------------- HGNC
class Hgnc:
    """Mise à jour des symboles : approuvé > ancien symbole (unique) > alias (unique)."""

    def __init__(self, path):
        self.approuve, self.prev, self.alias, self.info = set(), {}, {}, {}
        prev_multi, alias_multi = {}, {}
        with ouvrir(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r.get("status", "Approved") != "Approved":
                    continue
                s = r["symbol"]
                self.approuve.add(s)
                self.info[s] = r
                for col, d in (("prev_symbol", prev_multi), ("alias_symbol", alias_multi)):
                    for a in (r.get(col) or "").strip('"').split("|"):
                        if a:
                            d.setdefault(a, set()).add(s)
        self.prev = {k: next(iter(v)) for k, v in prev_multi.items() if len(v) == 1}
        self.alias = {k: next(iter(v)) for k, v in alias_multi.items() if len(v) == 1}
        self.ambigus = {k for k, v in prev_multi.items() if len(v) > 1}

    def maj(self, sym):
        """Renvoie (symbole_à_jour, statut)."""
        if sym in self.approuve:
            return sym, "approuve"
        if sym in self.prev:
            return self.prev[sym], "ancien_symbole"
        if sym in self.alias:
            return self.alias[sym], "alias"
        if sym in self.ambigus:
            return sym, "ambigu"
        return sym, "introuvable"

    def locus_type(self, sym):
        r = self.info.get(sym)
        return r["locus_type"] if r else ""


def lire_liste_genes(path, hgnc=None, colonne=None):
    """Liste de gènes : un symbole par ligne, ou TSV avec une colonne Symbol / symbol / gene (MitoCarta)."""
    genes = set()
    with ouvrir(path) as fh:
        first = fh.readline()
        cols = first.rstrip("\n").split("\t")
        col = None
        for c in ([colonne] if colonne else []) + ["Symbol", "symbol", "gene", "Gene Symbol", "gene_symbol"]:
            if c in cols:
                col = cols.index(c)
                break
        lignes = fh if col is not None else [first] + list(fh)
        for line in lignes:
            p = line.rstrip("\n").split("\t")
            s = (p[col] if col is not None and col < len(p) else p[0]).strip()
            if s and not s.startswith("#"):
                genes.add(hgnc.maj(s)[0] if hgnc else s)
    return genes


# --------------------------------------------------------------- statistique
def _lchoose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def hypergeom_sf(k, N, K, n):
    """P(X >= k), X ~ Hypergéométrique(N population, K succès, n tirages)."""
    if k <= 0:
        return 1.0
    hi = min(K, n)
    if k > hi:
        return 0.0
    den = _lchoose(N, n)
    return min(1.0, sum(math.exp(_lchoose(K, i) + _lchoose(N - K, n - i) - den) for i in range(k, hi + 1)))


# ------------------------------------------------------------------ variants
def lire_variants(path):
    """VCF (.vcf, .vcf.gz) ou TSV avec colonnes chr/chrom/CHROM, pos/POS, ref/REF, alt/ALT.
    Les ALT multiples sont éclatés. Renvoie des dict avec la ligne d'origine."""
    with ouvrir(path) as fh:
        header = None
        for line in fh:
            if line.startswith("##"):
                continue
            p = line.rstrip("\n").split("\t")
            if header is None:
                header = [h.lstrip("#") for h in p]
                low = [h.lower() for h in header]
                try:
                    ic = next(low.index(x) for x in ("chrom", "chr", "chromosome") if x in low)
                    ip, ir, ia = low.index("pos"), low.index("ref"), low.index("alt")
                except (StopIteration, ValueError):
                    sys.exit("Colonnes chr/pos/ref/alt introuvables dans %s" % path)
                continue
            for alt in p[ia].split(","):
                yield {"chrom": norm_chrom(p[ic]), "pos": int(p[ip]), "ref": p[ir].upper(),
                       "alt": alt.upper(), "ligne": p, "header": header}


def appliquer_variant(seq_genome, g_start, pos, ref, alt):
    """Applique REF>ALT (pos 1-based) à une séquence génomique brin + commençant en g_start (0-based).
    Renvoie (seq_alt, ref_concorde)."""
    i = pos - 1 - g_start
    if i < 0 or i + len(ref) > len(seq_genome):
        return None, False
    ok = seq_genome[i:i + len(ref)] == ref
    return seq_genome[:i] + alt + seq_genome[i + len(ref):], ok


def sortie(score, niveau, justification):
    assert niveau in NIVEAUX
    j = justification.replace("|", "/").replace("\t", " ").replace("\n", " ")
    return "%s|%s|%s" % (score, niveau, j)


def meilleur(resultats):
    """Combine plusieurs (score, niveau, justification) : garde le meilleur niveau, concatène les justifications."""
    if not resultats:
        return sortie(0, "aucun", "variant hors des éléments du chemin 1")
    resultats = sorted(resultats, key=lambda r: (NIVEAUX.index(r[1]), -float(r[0])))
    s, n, _ = resultats[0]
    return sortie(s, n, " ; ".join(r[2] for r in resultats))
