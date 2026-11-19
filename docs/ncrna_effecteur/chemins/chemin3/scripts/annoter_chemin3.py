#!/usr/bin/env python3
"""Chemin 3, annotation sur la VM (hors ligne, bibliothèque standard seulement).

Entrée : variants hg38 (VCF, ou TSV avec colonnes chr pos ref alt) ;
table sORF produite par preparer_sorf.py (*.table_sorf.tsv.gz).
Sortie : TSV des variants + colonne NC_chemin3 = score|niveau|justification.

Le variant est interprété comme codant dans le sORF : perte d'initiation, stop
gagné, décalage du cadre, stop perdu, indel en phase, faux-sens, site d'épissage
d'un sORF épissé. score = poids_conséquence x poids_niveau_mito (règles a priori,
sorf_lib.py). Aucune donnée ne sort de la machine.
"""
import argparse
import bisect
import sys
from collections import defaultdict

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import sorf_lib as L  # noqa: E402


class Index:
    def __init__(self, orfs):
        self.par_chr = defaultdict(list)
        for o in orfs:
            self.par_chr[o["chrom"]].append((o["exons"][0][0], o["exons"][-1][1], o))
        for c in self.par_chr:
            self.par_chr[c].sort(key=lambda x: x[0])
        self.debuts = {c: [x[0] for x in v] for c, v in self.par_chr.items()}
        self.maxlen = {c: max(x[1] - x[0] for x in v) for c, v in self.par_chr.items()}

    def chevauche(self, c, a, b):
        v = self.par_chr.get(c, [])
        if not v:
            return []
        i = bisect.bisect_left(self.debuts[c], a - self.maxlen[c] - 3)
        out = []
        while i < len(v) and v[i][0] < b + 3:
            if v[i][1] + 3 > a:
                out.append(v[i][2])
            i += 1
        return out


def pos_tx(o, g):
    """Index (0-based, sens du transcrit) du nt génomique g, ou None hors exon."""
    pos = 0
    ex = o["exons"] if o["strand"] == "+" else list(reversed(o["exons"]))
    for a, b in ex:
        if a <= g < b:
            return pos + (g - a if o["strand"] == "+" else b - 1 - g)
        pos += b - a
    return None


def site_epissage(o, a, b):
    """Le variant [a, b) touche-t-il les 2 nt introniques d'une jonction interne ?"""
    ex = o["exons"]
    for k in range(len(ex) - 1):
        i0, i1 = ex[k][1], ex[k + 1][0]  # intron [i0, i1)
        for x, y in ((i0, i0 + 2), (i1 - 2, i1)):
            if a < y and b > x:
                return True
    return False


def consequence(o, chrom, pos1, ref, alt):
    """Retourne (conséquence, détail) ou (None, raison)."""
    a = pos1 - 1
    b = a + len(ref)
    cds = o.get("cds_nt", "")
    if not cds:
        return None, "séquence du sORF absente de la table (préparer avec --genome)"
    dans = [pos_tx(o, g) for g in range(a, b)]
    if all(d is None for d in dans):
        if len(o["exons"]) > 1 and site_epissage(o, a, b):
            return "splice_site_orf", "site canonique d'une jonction du sORF épissé"
        return None, "hors exons du sORF"
    if any(d is None for d in dans):
        lo, hi = o["exons"][0][0], o["exons"][-1][1]
        if any(d is None and lo <= g < hi for d, g in zip(dans, range(a, b))):
            return "splice_site_orf", "variant à cheval sur une jonction exon-intron du sORF"
        touches = [d for d in dans if d is not None]
        if 0 in touches:
            return "start_lost", "variant chevauchant le début du sORF"
        return "stop_lost", "variant chevauchant la fin du sORF (codon stop)"
    i, j = min(dans), max(dans) + 1
    r, al = ref.upper(), alt.upper()
    if o["strand"] == "-":
        r, al = L.revcomp(r), L.revcomp(al)
    if cds[i:j].upper() != r:
        return None, "allèle de référence discordant avec la table (%s attendu)" % cds[i:j]
    nouv = cds[:i] + al + cds[j:]
    n_aa = len(cds) // 3 - (1 if L.traduire(cds[-3:]) == "*" else 0)
    codon = i // 3 + 1
    if i < 3:
        if nouv[:3].upper() != cds[:3].upper():
            return "start_lost", "codon d'initiation %s>%s" % (cds[:3], nouv[:3])
    diff = len(al) - len(r)
    if diff % 3 != 0:
        frac = 1 - (codon - 1) / max(n_aa, 1)
        return "frameshift", "décalage au codon %d/%d (%.0f %% du peptide touché)" % (codon, n_aa, 100 * frac)
    p_ref = L.traduire(cds)
    p_alt = L.traduire(nouv)
    stop_ref = p_ref.find("*")
    stop_alt = p_alt.find("*")
    if stop_ref >= 0 and i <= 3 * stop_ref + 2 and j > 3 * stop_ref:
        if stop_alt == -1 or stop_alt > stop_ref + diff // 3:
            return "stop_lost", "codon stop %s perdu" % cds[3 * stop_ref:3 * stop_ref + 3]
        return "stop_retained", "stop conservé"
    if stop_alt != -1 and (stop_ref == -1 or stop_alt < stop_ref + diff // 3):
        frac = 1 - stop_alt / max(n_aa, 1)
        return "stop_gained", "stop au codon %d/%d (%.0f %% du peptide retiré)" % (stop_alt + 1, n_aa, 100 * frac)
    if diff != 0:
        return "inframe_indel", "indel en phase de %+d aa au codon %d" % (diff // 3, codon)
    if p_ref == p_alt:
        return "synonymous", "codon %d synonyme" % codon
    k = next(x for x in range(len(p_ref)) if p_ref[x] != p_alt[x])
    return "missense", "%s%d%s" % (p_ref[k], k + 1, p_alt[k])


def poids_csq(csq, detail):
    w = L.POIDS_CSQ.get(csq, 0.0)
    if csq in ("stop_gained", "frameshift"):
        try:
            pct = float(detail.split("(")[1].split(" %")[0]) / 100
        except (IndexError, ValueError):
            pct = 1.0
        if pct < L.FRACTION_DISTALE:
            w = L.POIDS_TRONCATURE_DISTALE
    return w


def annoter(idx, chrom, pos, ref, alt):
    c = L.norm_chrom(chrom)
    if c in ("M", "MT"):
        return "0|aucun|chrM hors périmètre (humanine, MOTS-c exclus)"
    a = int(pos) - 1
    meilleurs = []
    for o in idx.chevauche(c, a, a + len(ref)):
        csq, det = consequence(o, c, int(pos), ref, alt)
        if csq is None:
            if "hors exons" not in det:
                meilleurs.append((0.0, "aucun", "%s : %s" % (o["orf_id"], det)))
            continue
        w = poids_csq(csq, det)
        wm = L.POIDS_MITO.get(o["niveau_mito"], 0.0)
        s = round(w * wm, 3)
        niv = L.niveau_depuis_score(s) if s > 0 else "aucun"
        j = "%s dans %s (hôte %s, %s aa) : %s ; mito %s : %s" % (
            csq, o["orf_id"], o["gene_hote"] or "?", o["longueur_aa"], det,
            o["niveau_mito"], o["justification_mito"])
        meilleurs.append((s, niv, j))
    if not meilleurs:
        return "0|aucun|aucun sORF traduit chevauché"
    meilleurs.sort(key=lambda x: -x[0])
    s, niv, j = meilleurs[0]
    if len(meilleurs) > 1:
        j += " ; %d autre(s) sORF chevauché(s)" % (len(meilleurs) - 1)
    return "%s|%s|%s" % (s, niv, j.replace("|", "/").replace("\t", " "))


def lire_variants(chemin):
    with L.ouvrir(chemin) as f:
        tete = None
        for l in f:
            if l.startswith("##"):
                continue
            c = l.rstrip("\n").split("\t")
            if tete is None:
                tete = [x.lstrip("#") for x in c]
                low = [x.lower() for x in tete]
                if "chrom" in low:
                    ic, ip, ir, ia = low.index("chrom"), low.index("pos"), low.index("ref"), low.index("alt")
                else:
                    ic, ip, ir, ia = (low.index(x) for x in ("chr", "pos", "ref", "alt"))
                yield "tete", tete, (ic, ip, ir, ia)
                continue
            yield "ligne", c, None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variants", required=True)
    ap.add_argument("--table", required=True, help="*.table_sorf.tsv.gz produit sur Curta")
    ap.add_argument("--sortie", required=True)
    a = ap.parse_args()
    idx = Index([o for o in L.lire_table(a.table) if o["exons"]])
    with open(a.sortie, "w") as out:
        cols = None
        for t, c, k in lire_variants(a.variants):
            if t == "tete":
                cols = k
                out.write("\t".join(c + ["NC_chemin3"]) + "\n")
                continue
            ic, ip, ir, ia = cols
            res = [annoter(idx, c[ic], c[ip], c[ir], x) for x in c[ia].split(",")]
            out.write("\t".join(c + [";;".join(res)]) + "\n")


if __name__ == "__main__":
    main()
