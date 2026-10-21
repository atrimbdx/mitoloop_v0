#!/usr/bin/env python3
"""
ajout_mitoloop.py — VM. Ajoute aux tables du pipeline reannot (*.final.tsv, ou
*.spliceai2.tsv) les arguments MitoLoop sur l'effet possible d'un variant non
codant sur l'expression d'un gène mitochondrial. Le fichier d'origine n'est
pas modifié : le résultat est écrit à côté, avec le suffixe .mitoloop.tsv.

Usage : python3 ajout_mitoloop.py run1.final.tsv [run2.final.tsv ...]
        (un VCF est aussi accepté : run.vcf.gz -> run.mitoloop.tsv)
        (REANNOT_BASE, MITOLOOP_DIR, FASTA, SAMTOOLS comme dans les autres
        scripts ; --data, --fasta, --samtools pour les donner autrement)

Seuls les variants situés dans les régions régulatrices des gènes
mitochondriaux du paquet (promoteur, UTR, éléments ENCODE-rE2G, ancres de
boucles) sont annotés ; les autres reçoivent le statut hors_regions_mito.

Colonnes ajoutées :
  MITO_gene          gène(s) mitochondrial(aux) dont une région contient le
                     variant, avec la classe de région (GENE:promoter|...)
  MITO_bidirectionnel  partenaire du promoteur bidirectionnel du gène (autre
                     gène, TSS sur le brin opposé à moins de 1 kb : distance,
                     type) : un variant du promoteur peut toucher les deux
  MITO_L1_classe     transmitted (l'intervalle à 95 % de L1 est au-dessus du
                     seuil haut : une baisse de l'ARNm se retrouve sur la
                     protéine), uncertain, not_transmitted, not_estimable
                     (le gène varie trop peu entre tissus pour être mesuré)
  MITO_L1, MITO_L1_IC   L1 et son intervalle à 95 % (protéome GTEx, Jiang 2020)
  MITO_L1_note       sous-unité tenue par son complexe (sa protéine suit le
                     complexe plus que son ARNm), corrélation ARN-protéine
                     dans les fibroblastes de patients (Kopajtich)
  MITO_dist_TSS      distance au site d'initiation validé le plus proche
                     (EPDnew), orientée sur son brin (négatif = en amont)
  MITO_CAGE          le variant est sur un pic CAGE (FANTOM5)
  MITO_fixation      facteurs fixés au variant (ChIP-seq ReMap 2022) et
                     nombre de types cellulaires où le pic est retrouvé
  MITO_motifs        sites JASPAR touchés (SNV) : facteur, position dans le
                     motif, contenu en information de la position (bits),
                     variation du score du site (log2, alt - ref)
  MITO_CpG           CpG et méthylation : CpG détruit ou créé par le
                     variant ; o/e CpG et GC de la fenêtre de 500 pb centrée
                     sur le variant, référence > variant, et sa classe
                     (HCP, ICP, LCP de Weber et al. 2007) ; sites touchés de
                     facteurs dont la fixation est empêchée par la
                     méthylation de leur CpG (NRF1, ETS/GABPA, CREB, E2F...)
  MITO_chromatine    état de la chromatine au variant (Roadmap ChromHMM)
  MITO_phyloP        conservation de la base sur 241 mammifères (Zoonomia,
                     phyloP ; positif = conservée)
  MITO_gnomAD        contraintes du gène dans gnomAD v4.1 (transcrit MANE
                     Select) : LoF observés/attendus, o/e LoF, LOEUF, pLI,
                     o/e et z faux-sens
  MITO_GISMO         GISMO (perte du gène parmi 462 mammifères) et GISMO-mis
                     (contrainte faux-sens entre mammifères), avec leur
                     décile (1 = le plus contraint ; déciles 1 à 4 de GISMO
                     pour les gènes récessifs selon Liao et al. 2026)
  MITO_arguments     arguments réunis (liste) et MITO_n_arguments
  MITO_evenements    événements portés par l'allèle lui-même, tous outils
                     réunis, avec leur force (voir Priorité)
  MITO_priorite      1_fort, 2_modere, 3_contexte, 4_aucun (voir Priorité)
  MITO_statut        annote, hors_regions_mito, NA

Arguments (chacun issu d'une donnée publiée, sans pondération) :
  L1_transmis           le gène transmet une baisse d'ARNm à la protéine
  promoteur_minimal     de -50 à +10 d'un TSS EPDnew
  pic_CAGE              sur un pic CAGE FANTOM5
  fixation_<TF>         pic ReMap retrouvé dans >= --min-biosamples types
  site_<TF>             le SNV touche une position informative (>= --min-ic
                        bits) d'un site du même facteur, fixé au variant, et
                        baisse son score d'au moins --min-delta
  chromatine_active_<tissu>  état actif (TssA, TssAFlnk, Enh, EnhG, TxFlnk)
  conservation_mammiferes    phyloP >= --min-phylop (2,27 : seuil de
                        conservation significative, FDR 5 %, de Zoonomia)
  CpG_site_<TF>         le variant détruit un CpG d'un site du facteur fixé
                        au variant (ReMap), facteur dont la fixation est
                        empêchée par la méthylation : perte du site et
                        possible méthylation (Domcke et al. 2015, NRF1)
  bascule_ilot_CpG      la fenêtre de 500 pb perd sa classe (HCP -> ICP ou
                        LCP, ICP -> LCP) : promoteur exposé à la méthylation
                        de novo (Weber et al. 2007 ; Lienert et al. 2011)
Les contraintes gnomAD et GISMO décrivent le gène, pas le variant : elles sont
données en colonnes, sans entrer dans les arguments.

Priorité (MITO_evenements, MITO_priorite). Les annotations du gène et de la
position sont les mêmes pour un variant pathogène et un variant bénin au même
endroit ; seules celles de l'allèle les séparent. Les événements de l'allèle,
de tous les outils de la table, sont donc réunis et classés :
  fort     épissage : SpliceAI ou SpliceAI2 >= --spliceai-fort (0,5 ; 0,2 pour
                       modéré : seuils publiés par les auteurs de SpliceAI,
                       fixés à la main ; selon les conditions d'utilisation,
                       ses sorties ne servent ni à entraîner ni à ajuster)
           uORF/dORF : un champ MORFEE renseigné (création ou perte d'un
                       uTIS/uSTOP, dTIS/dSTOP)
           AlphaGenome : |AG_*| >= --alphagenome-fort (seulement si donné ;
                       seuil fixé à la main : selon leurs conditions
                       d'utilisation, les sorties d'AlphaGenome ne servent à
                       entraîner ni à ajuster aucun paramètre)
           site_<TF> d'au moins --site-fort-delta (2, log2) sur une position
                       d'au moins --site-fort-ic (1,5 bits)
           un événement modéré dans le promoteur minimal (-50 à +10 d'un TSS)
  modéré   SpliceAI ou SpliceAI2 >= --spliceai-modere (0,2) ; site_<TF> ;
           CpG_site_<TF> ; bascule_ilot_CpG
  1_fort      au moins un événement fort, dans un gène dont la classe L1 n'est
              pas not_transmitted
  2_modere    un événement modéré, ou fort dans un gène not_transmitted
  3_contexte  pas d'événement de l'allèle, mais au moins deux arguments de
              position (promoteur minimal, CAGE, fixation, chromatine active,
              conservation)
  4_aucun     le reste
Les variants hors des régions MitoLoop reçoivent aussi les événements des
outils de la table (épissage, MORFEE, AlphaGenome) et leur priorité.

Python 3 de la bibliothèque standard et samtools (environnement vep).
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import re
import shutil
import subprocess
import sys
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path

ACTIVE = re.compile(r"(^|_)(TssA|TssAFlnk|Enh|EnhG|TxFlnk)$")
# facteurs dont la fixation est empêchée par la méthylation du CpG de leur site
METH_SENSIBLE = {"NRF1", "GABPA", "ELK1", "ELK4", "ELF1", "ETS1", "ETV6", "CREB1",
                 "ATF1", "E2F1", "E2F4", "E2F6", "ZBTB33", "NFYA"}
W_CPG = 250                       # demi-fenêtre CpG (fenêtre de 500 pb)
CLASS_ORDER = ["promoter", "utr5", "element", "anchor", "utr3"]
COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def opener(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def base_dir() -> Path:
    b = os.environ.get("REANNOT_BASE")
    if not b:
        u = os.environ.get("USER", "")
        b = f"/scratch/{u}/reannot" if os.path.isdir(f"/scratch/{u}") else \
            os.path.expanduser("~/reannot")
    return Path(b)


class Intervals:
    """intervalles par chromosome (chrom, start 0-based, end, étiquette)"""

    def __init__(self, path):
        d = defaultdict(list)
        with opener(path) as fh:
            for line in fh:
                if not line.strip() or line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                d[chrom_key(f[0])].append((int(f[1]), int(f[2]), f[3:]))
        self.d = {}
        for c, v in d.items():
            v.sort(key=lambda x: x[0])
            st = [x[0] for x in v]
            maxlen = max(x[1] - x[0] for x in v)
            self.d[c] = (st, v, maxlen)

    def query(self, chrom, s, e):
        x = self.d.get(chrom_key(chrom))
        if not x:
            return []
        st, v, maxlen = x
        lo = bisect_left(st, s - maxlen)
        hi = bisect_right(st, e - 1)
        return [r for r in v[lo:hi] if r[1] > s]


def chrom_key(c: str) -> str:
    c = c[3:] if c.startswith("chr") else c
    return "MT" if c in ("M", "MT") else c


def read_motifs(path):
    """matrices JASPAR (format jaspar) : nom -> colonnes de log2(p / 0.25)"""
    mats, name, rows = {}, None, {}
    if not Path(path).exists():
        return mats

    def close():
        if name and len(rows) == 4:
            n = len(rows["A"])
            cols = []
            for i in range(n):
                cnt = [rows[b][i] for b in "ACGT"]
                tot = sum(cnt) + 1.0
                p = [(x + 0.25) / tot for x in cnt]
                cols.append([math.log2(q / 0.25) for q in p] +
                            [2 + sum(q * math.log2(q) for q in p)])
            mats[name] = cols

    with opener(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                close()
                f = line[1:].split()
                name, rows = (f[1] if len(f) > 1 else f[0]).upper(), {}
            elif line and line[0] in "ACGT":
                nums = re.findall(r"-?\d+(?:\.\d+)?", line[1:])
                rows[line[0]] = [float(x) for x in nums]
    close()
    return mats


def score(cols, seq):
    s = 0.0
    for c, b in zip(cols, seq):
        k = "ACGT".find(b)
        if k < 0:
            return None
        s += c[k]
    return s


def motif_hits(mats, tfs, seq, off, ref, alt, min_rel):
    """sites touchés par un SNV : seq autour du variant, off = index du
    variant dans seq ; (tf, position 1-based dans le motif, IC, delta)"""
    out = []
    alt_seq = seq[:off] + alt + seq[off + 1:]
    for tf in tfs:
        cols = mats.get(tf)
        if not cols:
            continue
        L = len(cols)
        smax = sum(max(c[:4]) for c in cols)
        smin = sum(min(c[:4]) for c in cols)
        best = None
        for strand in "+-":
            for st in range(off - L + 1, off + 1):
                if st < 0 or st + L > len(seq):
                    continue
                w, wa = seq[st:st + L], alt_seq[st:st + L]
                if strand == "-":
                    w = w.translate(COMP)[::-1]
                    wa = wa.translate(COMP)[::-1]
                    pos = st + L - 1 - off
                else:
                    pos = off - st
                s_ref, s_alt = score(cols, w), score(cols, wa)
                if s_ref is None or s_alt is None:
                    continue
                rel = (s_ref - smin) / (smax - smin) if smax > smin else 0
                if rel >= min_rel and (best is None or s_ref > best[0]):
                    best = (s_ref, pos, cols[pos][4], s_alt - s_ref)
        if best:
            out.append((tf, best[1] + 1, best[2], best[3]))
    return out


def fetch_seqs(fasta, samtools, regions):
    """séquences de la référence pour une liste de régions 'chrom:s-e'"""
    seqs = {}
    if not regions:
        return seqs
    fai = Path(str(fasta) + ".fai")
    names, fa = set(), None
    if not shutil.which(samtools):
        # sans samtools (Curta, environnement mitoloop) : pysam
        import pysam
        fa = pysam.FastaFile(str(fasta))
        names = set(fa.references)
    elif fai.exists():
        names = {ln.split("\t")[0] for ln in fai.read_text().splitlines()}
    req = []
    for r in regions:
        c, rest = r.rsplit(":", 1)
        k = chrom_key(c)
        name = k if k in names else ("chr" + k if "chr" + k in names else
                                     ("chrM" if k == "MT" and "chrM" in names
                                      else k))
        req.append((r, f"{name}:{rest}"))
    if fa is not None:
        for r, q in req:
            c, rest = q.rsplit(":", 1)
            s, e = (int(x) for x in rest.split("-"))
            try:
                seqs[r] = fa.fetch(c, s - 1, e).upper()
            except (KeyError, ValueError):
                seqs[r] = ""
        return seqs
    tmp = Path(f".mitoloop_regions_{os.getpid()}.txt")
    tmp.write_text("".join(q + "\n" for _r, q in req))
    try:
        res = subprocess.run([samtools, "faidx", str(fasta), "-r", str(tmp)],
                             capture_output=True, text=True, check=True)
    finally:
        tmp.unlink(missing_ok=True)
    cur, buf, out = None, [], {}
    for line in res.stdout.splitlines():
        if line.startswith(">"):
            if cur is not None:
                out[cur] = "".join(buf).upper()
            cur, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if cur is not None:
        out[cur] = "".join(buf).upper()
    for r, q in req:
        seqs[r] = out.get(q, "")
    return seqs


def load_genes(path):
    g = {}
    with opener(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            r = dict(zip(head, f))
            g[r["gene"]] = r
    return g


def main():
    p = argparse.ArgumentParser(description="Arguments MitoLoop pour les "
                                "variants non codants des gènes mitochondriaux")
    p.add_argument("tables", nargs="+")
    b = base_dir()
    p.add_argument("--data", type=Path, default=Path(
        os.environ.get("MITOLOOP_DIR", b / "data" / "mitoloop")))
    p.add_argument("--fasta", default=os.environ.get(
        "FASTA", str(b / "data/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")))
    p.add_argument("--samtools", default=os.environ.get("SAMTOOLS", "samtools"))
    p.add_argument("--min-biosamples", type=int, default=5)
    p.add_argument("--min-rel", type=float, default=0.8,
                   help="score relatif minimal d'un site JASPAR sur la référence")
    p.add_argument("--min-ic", type=float, default=1.0)
    p.add_argument("--min-delta", type=float, default=1.0)
    p.add_argument("--min-phylop", type=float, default=2.27)
    p.add_argument("--spliceai-fort", type=float, default=0.5)
    p.add_argument("--spliceai-modere", type=float, default=0.2)
    p.add_argument("--site-fort-delta", type=float, default=2.0)
    p.add_argument("--site-fort-ic", type=float, default=1.5)
    p.add_argument("--alphagenome-fort", type=float,
                   help="|score AlphaGenome| d'un événement fort (sans cette option, "
                        "AlphaGenome n'entre pas dans la priorité)")
    p.add_argument("--tabix", default=os.environ.get("TABIX", "tabix"))
    a = p.parse_args()
    D = a.data
    genes = load_genes(D / "genes.tsv")
    regions = Intervals(D / "regions.bed.gz")
    ann = {}
    for f in sorted((D / "annotations").glob("*.bed.gz")):
        ann[f.name[:-7]] = Intervals(f)
    mats = read_motifs(D / "motifs.jaspar")
    tfs = sorted(mats)
    chrom_ann = [k for k in ann if k.startswith("chromhmm_")]
    print(f"MitoLoop : {len(genes)} gènes, {len(ann)} annotations, "
          f"{len(mats)} motifs ({', '.join(tfs)})", file=sys.stderr)

    for path in a.tables:
        out = re.sub(r"\.tsv$", "", path) + ".mitoloop.tsv"
        if re.search(r"\.vcf(\.gz)?$", path):
            # VCF : une ligne par allèle alternatif
            head, rows = ["chrom", "pos", "id", "ref", "alt"], []
            with opener(path) as fh:
                for ln in fh:
                    if ln.startswith("#"):
                        continue
                    f = ln.rstrip("\n").split("\t")
                    # les allèles symboliques (<DEL>, SV) ne sont pas annotés
                    rows += [[f[0], f[1], f[2], f[3], x] for x in f[4].split(",")
                             if re.fullmatch(r"[ACGTNacgtn]+", x)]
            out = re.sub(r"\.vcf(\.gz)?$", "", path) + ".mitoloop.tsv"
        else:
            with open(path) as fh:
                head = fh.readline().rstrip("\n").split("\t")
                rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
        col = {h.lower().lstrip("#"): i for i, h in enumerate(head)}
        for k in ("chrom", "pos", "ref", "alt"):
            if k not in col:
                sys.exit(f"{path} : colonne {k} introuvable dans l'en-tête")
        C, P, R, A = col["chrom"], col["pos"], col["ref"], col["alt"]
        G = col.get("gene")

        # 1. régions et séquences à lire (SNV dans les régions seulement)
        hits_of, need = {}, set()
        for i, f in enumerate(rows):
            if not f[P].isdigit():
                continue
            pos, ref = int(f[P]), f[R]
            h = regions.query(f[C], pos - 1, pos - 1 + max(len(ref), 1))
            if not h:
                continue
            hits_of[i] = h
            if re.fullmatch(r"[ACGTacgt]+", ref) and re.fullmatch(r"[ACGTacgt]+", f[A]):
                need.add(f"{f[C]}:{max(1, pos - W_CPG)}-{pos + len(ref) + W_CPG}")
        seqs = fetch_seqs(a.fasta, a.samtools, sorted(need)) if need else {}
        phy = fetch_phylop(D / "annotations" / "phylop.bedgraph.gz", a.tabix,
                           {(rows[i][C], int(rows[i][P])) for i in hits_of})

        # 2. annotation
        extra = ["MITO_gene", "MITO_bidirectionnel", "MITO_L1_classe", "MITO_L1",
                 "MITO_L1_IC", "MITO_L1_note", "MITO_gnomAD", "MITO_GISMO",
                 "MITO_dist_TSS", "MITO_CAGE", "MITO_fixation", "MITO_motifs",
                 "MITO_CpG", "MITO_chromatine", "MITO_phyloP", "MITO_arguments",
                 "MITO_n_arguments", "MITO_evenements", "MITO_priorite",
                 "MITO_statut"]
        lowhead = [h.lower() for h in head]
        n_stat = defaultdict(int)
        with open(out, "w") as fo:
            fo.write("\t".join(head + extra) + "\n")
            for i, f in enumerate(rows):
                if not f[P].isdigit():
                    fo.write("\t".join(f + ["."] * (len(extra) - 1) + ["NA"]) + "\n")
                    n_stat["NA"] += 1
                    continue
                if i not in hits_of:
                    ev, pr = priority(dict(zip(lowhead, f)), {}, [], a)
                    fo.write("\t".join(f + ["."] * (len(extra) - 3) +
                                       [ev, pr, "hors_regions_mito"]) + "\n")
                    n_stat["hors_regions_mito"] += 1
                    continue
                v = annotate(f, C, P, R, A, G, hits_of[i], genes, ann,
                             chrom_ann, mats, tfs, seqs, a, phy)
                v["MITO_evenements"], v["MITO_priorite"] = priority(
                    dict(zip(lowhead, f)), v, v.pop("_events"), a)
                fo.write("\t".join(f + [v[k] for k in extra[:-1]] + ["annote"]) + "\n")
                n_stat["annote"] += 1
        for s_, n in sorted(n_stat.items()):
            print(f"           {s_:18s} {n}", file=sys.stderr)
        print(f"           Résultat : {out}", file=sys.stderr)


def fetch_phylop(path, tabix, positions):
    """(chrom, pos) -> phyloP, lu dans le bedGraph indexé (tabix, ou pysam)"""
    out = {}
    if not positions or not Path(path).exists():
        return out
    by_key = {}
    for c, p_ in positions:
        by_key[(chrom_key(c), p_)] = (c, p_)
    lines = []
    if shutil.which(tabix):
        tmp = Path(f".mitoloop_pos_{os.getpid()}.txt")
        names = None
        try:
            res = subprocess.run([tabix, "-l", str(path)], capture_output=True,
                                 text=True, check=True)
            names = {chrom_key(x): x for x in res.stdout.split()}
        except subprocess.CalledProcessError:
            names = {}
        tmp.write_text("".join(f"{names.get(k, k)}\t{p_}\n"
                               for k, p_ in sorted(by_key)))
        try:
            res = subprocess.run([tabix, "-R", str(tmp), str(path)],
                                 capture_output=True, text=True, check=True)
            lines = res.stdout.splitlines()
        finally:
            tmp.unlink(missing_ok=True)
    else:
        import pysam
        t = pysam.TabixFile(str(path))
        names = {chrom_key(x): x for x in t.contigs}
        for k, p_ in by_key:
            if k in names:
                lines += list(t.fetch(names[k], p_ - 1, p_))
    for ln in lines:
        f = ln.split("\t")
        k, s, e = chrom_key(f[0]), int(f[1]), int(f[2])
        for p_ in range(s + 1, e + 1):
            if (k, p_) in by_key:
                out[by_key[(k, p_)]] = float(f[3])
    return out


def cpg_stats(w):
    """o/e CpG (Gardiner-Garden et Frommer), GC et classe de Weber et al.
    2007 d'une fenêtre"""
    n = sum(w.count(b) for b in "ACGT")
    if n < 100:
        return None
    c, g, cg = w.count("C"), w.count("G"), w.count("CG")
    oe = cg * n / (c * g) if c and g else 0.0
    gc = (c + g) / n
    cls = "HCP" if oe > 0.75 and gc > 0.55 else "LCP" if oe <= 0.48 else "ICP"
    return oe, gc, cls


def cpg_context(seq, off, ref, alt):
    """CpG détruit ou créé, et fenêtre de 500 pb avant et après le variant ;
    seq : référence autour du variant, off : indice du variant"""
    ref, alt = ref.upper(), alt.upper()
    if seq[off:off + len(ref)] != ref:
        return None
    alt_seq = seq[:off] + alt + seq[off + len(ref):]
    lo, hi = max(0, off - 1), off + len(ref) + 1
    n_ref = seq[lo:hi].count("CG")
    n_alt = alt_seq[lo:off + len(alt) + 1].count("CG")
    a = max(0, off - W_CPG)
    r = cpg_stats(seq[a:off + len(ref) + W_CPG])
    v = cpg_stats(alt_seq[a:off + len(alt) + W_CPG])
    return {"change": "detruit" if n_alt < n_ref else "cree" if n_alt > n_ref else ".",
            "ref": r, "alt": v}


def gene_constraints(gi):
    """colonnes MITO_gnomAD et MITO_GISMO"""
    gn = []
    if gi.get("gnomad_lof_obs", "") != "" or gi.get("gnomad_loeuf", ""):
        gn.append(f"LoF={gi.get('gnomad_lof_obs', '')}/{gi.get('gnomad_lof_exp', '')}")
        for k, lab in (("gnomad_lof_oe", "oe_LoF"), ("gnomad_loeuf", "LOEUF"),
                       ("gnomad_pli", "pLI"), ("gnomad_mis_oe", "oe_mis"),
                       ("gnomad_mis_z", "z_mis")):
            if gi.get(k, ""):
                gn.append(f"{lab}={gi[k]}")
    gm = []
    for k, lab in (("gismo", "GISMO"), ("gismo_mis", "GISMO-mis")):
        if gi.get(k, ""):
            d = gi.get(k + "_decile", "")
            gm.append(f"{lab}={gi[k]}" + (f"(d{d})" if d else ""))
    return ";".join(gn) or ".", ";".join(gm) or "."


def annotate(f, C, P, R, A, G, hits, genes, ann, chrom_ann, mats, tfs, seqs, a,
             phy=None):
    chrom, pos, ref, alt = f[C], int(f[P]), f[R], f[A]
    by_gene = defaultdict(list)
    for _s, _e, lab in hits:
        by_gene[lab[0]].append(lab[1])
    gl = sorted(by_gene, key=lambda g: (G is None or f[G] != g,
                                        min(CLASS_ORDER.index(c) if c in CLASS_ORDER
                                            else 9 for c in by_gene[g]), g))
    main_g = gl[0]
    gi = genes.get(main_g, {})
    args = []
    v = {"MITO_gene": "|".join(f"{g}:{','.join(sorted(set(by_gene[g])))}"
                               for g in gl)}
    cls = gi.get("l1_class", "not_estimable") or "not_estimable"
    v["MITO_L1_classe"] = cls
    v["MITO_L1"] = gi.get("l1", ".") or "."
    lo, hi = gi.get("ci_low", ""), gi.get("ci_high", "")
    v["MITO_L1_IC"] = f"[{lo};{hi}]" if lo and hi else "."
    notes = []
    try:
        if float(gi.get("protein_coupling", "nan")) - \
                float(gi.get("rna_coupling", "nan")) > 0.2:
            notes.append(f"sous-unité tenue par {gi.get('complex', 'son complexe')}")
    except ValueError:
        pass
    if gi.get("fibroblast_r", ""):
        notes.append(f"r_fibroblastes={gi['fibroblast_r']}")
    v["MITO_L1_note"] = ";".join(notes) or "."
    v["MITO_bidirectionnel"] = "|".join(f"{g}<>{genes[g]['bidirectional']}" for g in gl
                                        if genes.get(g, {}).get("bidirectional")) or "."
    events = []                       # (force, texte) des événements de l'allèle
    if cls == "transmitted":
        args.append("L1_transmis")
    v["MITO_gnomAD"], v["MITO_GISMO"] = gene_constraints(gi)

    s0, e0 = pos - 1, pos - 1 + max(len(ref), 1)
    # site d'initiation validé (EPDnew) le plus proche, orienté
    best = None
    if "epd" in ann:
        # annotations/epd.bed.gz : la base du TSS (0-based) et son brin
        for tss, _e, lab in ann["epd"].query(chrom, s0 - 1000, e0 + 1000):
            strand = lab[0] if lab else "+"
            d = (pos - 1 - tss) if strand == "+" else (tss - (pos - 1))
            if best is None or abs(d) < abs(best):
                best = d
    v["MITO_dist_TSS"] = "." if best is None else str(best)
    if best is not None and -50 <= best <= 10:
        args.append("promoteur_minimal")
    cage = "cage" in ann and bool(ann["cage"].query(chrom, s0, e0))
    v["MITO_CAGE"] = "oui" if cage else "non" if "cage" in ann else "."
    if cage:
        args.append("pic_CAGE")
    bound = {}
    for k in ann:
        if not k.startswith("remap_"):
            continue
        tf = k[6:].upper()
        for _s, _e, lab in ann[k].query(chrom, s0, e0):
            try:
                n = int(lab[0])
            except (ValueError, IndexError):
                n = 1
            bound[tf] = max(bound.get(tf, 0), n)
    v["MITO_fixation"] = ",".join(f"{t}({n})" for t, n in sorted(bound.items())) or "."
    for t, n in sorted(bound.items()):
        if n >= a.min_biosamples:
            args.append(f"fixation_{t}")
    mot = []
    key = f"{chrom}:{max(1, pos - W_CPG)}-{pos + len(ref) + W_CPG}"
    seq = seqs.get(key, "")
    off = pos - max(1, pos - W_CPG)
    ref_ok = bool(seq) and seq[off:off + len(ref)] == ref.upper()
    sites = []
    if len(ref) == 1 and len(alt) == 1 and mats:
        if ref_ok:
            for tf, mp, ic, dlt in motif_hits(mats, tfs, seq, off, ref.upper(),
                                              alt.upper(), a.min_rel):
                mot.append(f"{tf}:p{mp}/IC{ic:.1f}/d{dlt:+.1f}")
                sites.append(tf)
                if ic >= a.min_ic and dlt <= -a.min_delta and tf in bound:
                    args.append(f"site_{tf}")
                    strong = ic >= a.site_fort_ic and dlt <= -a.site_fort_delta
                    events.append(("fort" if strong else "modere",
                                   f"site_{tf}(IC{ic:.1f},d{dlt:+.1f})"))
        elif seq:
            mot.append("REF_discordante")
    v["MITO_motifs"] = ",".join(mot) or "."
    # CpG et méthylation
    cpg = cpg_context(seq, off, ref, alt) if ref_ok else None
    if cpg:
        parts = [f"CpG_{cpg['change']}"] if cpg["change"] != "." else []
        r_, a_ = cpg["ref"], cpg["alt"]
        if r_ and a_:
            parts.append(f"oe={r_[0]:.2f}>{a_[0]:.2f};GC={r_[1]:.2f}>{a_[1]:.2f};{r_[2]}>{a_[2]}")
            order = {"HCP": 2, "ICP": 1, "LCP": 0}
            if order[a_[2]] < order[r_[2]]:
                args.append("bascule_ilot_CpG")
                events.append(("modere", f"bascule_ilot_CpG({r_[2]}>{a_[2]})"))
        if cpg["change"] == "detruit":
            sens = [t for t in dict.fromkeys(sites) if t in METH_SENSIBLE]
            if sens:
                parts.append("sites=" + ",".join(sens))
            for t in sens:
                if t in bound:
                    args.append(f"CpG_site_{t}")
                    events.append(("modere", f"CpG_site_{t}"))
        v["MITO_CpG"] = ";".join(parts) or "."
    else:
        v["MITO_CpG"] = "REF_discordante" if seq and not ref_ok else "."
    states = []
    for k in chrom_ann:
        tissue = k[len("chromhmm_"):]
        for _s, _e, lab in ann[k].query(chrom, s0, s0 + 1):
            st = lab[0]
            states.append(f"{tissue}:{st}")
            if ACTIVE.search(st):
                args.append(f"chromatine_active_{tissue}")
    v["MITO_chromatine"] = ",".join(states) or "."
    x = (phy or {}).get((chrom, pos))
    v["MITO_phyloP"] = "." if x is None else f"{x:.2f}"
    if x is not None and x >= a.min_phylop:
        args.append("conservation_mammiferes")
    v["MITO_arguments"] = ",".join(dict.fromkeys(args)) or "."
    v["MITO_n_arguments"] = str(len(dict.fromkeys(args)))
    v["_events"] = events
    return v


POSITION_ARGS = ("promoteur_minimal", "pic_CAGE", "fixation_", "chromatine_active_",
                 "conservation_mammiferes")


def fnum(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if v == v else None


def priority(row, v, events, a):
    """événements de l'allèle (MitoLoop et outils de la table) et priorité ;
    row : colonnes de la table en minuscules"""
    ev = list(events)
    # épissage : SpliceAI (VEP) et SpliceAI2
    for tool, keys in (("SpliceAI", ("spliceai_ds_max", "spliceai_ds_ag", "spliceai_ds_al",
                                     "spliceai_ds_dg", "spliceai_ds_dl")),
                       ("SpliceAI2", ("spliceai2_max", "spliceai2_ds_ag", "spliceai2_ds_al",
                                      "spliceai2_ds_dg", "spliceai2_ds_dl"))):
        x = [fnum(row.get(k)) for k in keys]
        x = [y for y in x if y is not None]
        if x and max(x) >= a.spliceai_modere:
            ev.append(("fort" if max(x) >= a.spliceai_fort else "modere",
                       f"{tool}={max(x):.2f}"))
    # MORFEE : création ou perte d'un uTIS / uSTOP / dTIS / dSTOP
    mf = [k[len("morfee_"):] for k, x in row.items() if k.startswith("morfee_")
          and x not in (None, "", ".", "NA", "na", "-")]
    if mf:
        ev.append(("fort", "MORFEE:" + "+".join(mf)))
    # AlphaGenome : seulement avec un seuil
    if a.alphagenome_fort is not None:
        x = [(k, fnum(y)) for k, y in row.items() if k.startswith("ag_")]
        x = [(k, y) for k, y in x if y is not None]
        if x:
            k, y = max(x, key=lambda t: abs(t[1]))
            if abs(y) >= a.alphagenome_fort:
                ev.append(("fort", f"AlphaGenome:{k[3:]}={y:.3g}"))
    args = (v.get("MITO_arguments") or ".").split(",")
    if "promoteur_minimal" in args:      # un événement du promoteur minimal compte comme fort
        ev = [("fort" if s == "modere" and not t.startswith(("SpliceAI", "MORFEE"))
               else s, t) for s, t in ev]
    forts = [t for s, t in ev if s == "fort"]
    moderes = [t for s, t in ev if s == "modere"]
    n_pos = sum(1 for x in args if x.startswith(POSITION_ARGS))
    if forts and v.get("MITO_L1_classe") != "not_transmitted":
        pr = "1_fort"
    elif forts or moderes:
        pr = "2_modere"
    elif n_pos >= 2:
        pr = "3_contexte"
    else:
        pr = "4_aucun"
    txt = ",".join([f"{t}[fort]" for t in forts] + [f"{t}[modere]" for t in moderes])
    return txt or ".", pr


if __name__ == "__main__":
    main()
