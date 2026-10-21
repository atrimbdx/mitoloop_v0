#!/usr/bin/env python3
"""
mpra_variants.py — MCIA. Variants non codants dont l'effet sur l'activité d'un
promoteur, d'un enhancer ou d'un UTR a été mesuré expérimentalement (MPRA,
mutagenèse à saturation), pour tester les annotations de MitoLoop qui portent
sur l'allèle (classe C : variation du score d'un site JASPAR, site_<TF>,
CpG, MORFEE) contre celles qui portent sur la position
(classe B : phyloP, CAGE, fixation ReMap, chromatine). SpliceAI, SpliceAI2 et
AlphaGenome n'y sont pas évalués : selon leurs conditions d'utilisation, leurs
sorties ne servent ni à entraîner ni à ajuster des seuils.

Sources (téléchargées dans --datadir si absentes ; une source illisible est
sautée, avec un message)
  kircher    Kircher et al. 2019, Nat Commun 10:3583 : mutagenèse à saturation
             de 20 promoteurs et enhancers (TERT, LDLR, F9, HBB, PKLR...) ;
             toutes les substitutions et délétions d'une base ; Coefficient
             (log2) et pValue ; GRCh38 donné. Significatif : p < 1e-5 (seuil
             du plugin VEP des auteurs) ; sans effet : p > 0,05.
  sure       van Arensbergen et al. 2019, Nat Genet 51:1160 : SuRE, 5,9
             millions de SNP en K562 et HepG2 (hg19, converti par la chaîne
             UCSC hg19ToHg38). Significatifs : les raQTL publiés (FDR 5 %) ;
             sans effet : SNP de la table complète couverts par 10 à 999
             fragments pour chaque allèle, p de Wilcoxon > 0,5, dans un
             élément actif (signal au moins égal au 5e centile de celui des
             raQTL de la lignée), tirés au hasard (--ratio-negatifs par raQTL).
             Effet : log2((alt + 1) / (ref + 1)) des signaux moyens.
  mpravardb  MPRAVarDB (Bioinformatics 2024) : 18 études MPRA de variants
             (eQTL, GWAS...), en hg38 par le bigBed de l'UCSC. Significatif :
             FDR < 0,05 ; sans effet : FDR >= 0,5 (ou p > 0,5 sans FDR).
  sample     Sample et al. 2019, Nat Biotechnol 37:803 : charge en ribosomes
             (MRL) de 5'UTR humains et de leurs SNV ; effet = MRL alt - MRL
             du 5'UTR de référence ; pas de significativité.
  mqtl       mQTL du muscle squelettique de GTEx (Oliva et al. 2023, Nat Genet
             55:112 ; puces EPIC, hg38) : variant de tête de chaque CpG,
             effet = pente sur la méthylation (normalisée). Significatif :
             q < 0,05 ; sans effet : q > 0,5. Seuls les variants à plus de
             50 pb (ils peuvent fausser l'hybridation de la sonde) et à
             moins de --mqtl-distance pb (1000) du CpG mesuré sont gardés,
             coordonnées des CpG par le manifeste EPIC hg38 de Zhou et al.
             Teste les annotations de méthylation (MITO_CpG) ; « element »
             = le CpG mesuré.
  --table    toute autre table (Griesemer 2021 3'UTR, Tewhey 2016, Abell 2022,
             tableaux d'un article...) : TSV avec en-tête chrom, pos, ref, alt
             et, si connus, effet, p, signif (1/0), element, cellule ;
             --table-hg19 pour une table en hg19 (convertie). Nom de source :
             le nom du fichier.

Sous-commandes
  build     sources -> --out (TSV) : chrom pos ref alt source element cellule
            effet p signif ; base de référence contrôlée sur --fasta (GRCh38),
            les discordances sont écartées et comptées ; --max variants par
            source au plus (positions tirées au hasard, tous les allèles
            d'une position gardés ensemble).
  regions   BED d'évaluation : une région « element » autour de chaque
            variant, pour faire avec mitoloop_paquet.py un paquet qui annote
            tous les variants (et pas seulement ceux des gènes
            mitochondriaux) ; écrit aussi la liste de gènes et un l1.tsv vide
            demandés par mitoloop_paquet.py.
  entree    variants situés dans les régions MitoLoop (regions.bed de
            vcf_annotate.py) -> entrée de reannot_complet.sh sur la VM
            (VEP, SpliceAI, MORFEE, AlphaGenome, MitoLoop).
  evaluer   tables annotées (ajout_mitoloop.py, ou sortie de la VM) contre
            les effets mesurés : par source et lignée, AUC significatifs
            contre sans effet de chaque annotation, AUC appariée par
            position (paires d'allèles d'une même base : les annotations de
            position y valent 0,5 par construction, seules celles de l'allèle
            peuvent s'en écarter), corrélation de Spearman avec |effet| et,
            pour la variation du score des sites, avec l'effet signé.

Usage
  python3 scripts/mpra_variants.py build --datadir data/mpra --fasta data/GRCh38.fa \\
      --out results/mpra/mpra_variants.tsv
  python3 scripts/mpra_variants.py regions --variants results/mpra/mpra_variants.tsv \\
      --out results/mpra/eval
  python3 scripts/mpra_variants.py entree --variants results/mpra/mpra_variants.tsv \\
      --regions results/reannot/regions.bed --out results/mpra/mpra_mito.tsv
  python3 scripts/mpra_variants.py evaluer --variants results/mpra/mpra_variants.tsv \\
      --annote results/mpra/mpra_variants.mitoloop.tsv --out results/mpra/evaluation.tsv
Bibliothèque standard ; pyBigWig pour le bigBed de MPRAVarDB.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import random
import re
import sys
import urllib.request
from bisect import bisect_left, bisect_right
from collections import defaultdict
from pathlib import Path

URLS = {
    "kircher": "https://raw.githubusercontent.com/kircherlab/MPRA_SaturationMutagenesis/"
               "master/data/elements.tsv.gz",
    "sure_k562": "https://raw.githubusercontent.com/Mele-Lab/2023_BloodTraitsMPRA_CellGenomics/"
                 "main/Hypertension_paper/data/Nature_Gen_SNPs/k562.sign.id.LP190708.txt.gz",
    "sure_hepg2": "https://raw.githubusercontent.com/Mele-Lab/2023_BloodTraitsMPRA_CellGenomics/"
                  "main/Hypertension_paper/data/Nature_Gen_SNPs/hepg2.sign.id.LP190708.txt.gz",
    "sure_table": "https://osf.io/vxfk3/download",
    "chain": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    "mpravardb": "https://hgdownload.soe.ucsc.edu/gbdb/hg38/mpra/mpravardb/mpravardb.bb",
    "sample": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3130nnn/GSM3130443/suppl/"
              "GSM3130443_designed_library.csv.gz",
    "mqtl": "https://storage.googleapis.com/egtex/methylation/epic-arrays/mQTLs/"
            "MuscleSkeletal.regular.perm.fdr.txt",
    "epic": "https://raw.githubusercontent.com/zhou-lab/InfiniumAnnotationV1/main/Anno/EPIC/"
            "EPIC.hg38.manifest.tsv.gz",
}
FILES = {"kircher": "kircher2019_elements.tsv.gz", "sure_k562": "sure_k562.sign.txt.gz",
         "sure_hepg2": "sure_hepg2.sign.txt.gz", "sure_table": "SuRE_SNP_table_LP190708.txt.gz",
         "chain": "hg19ToHg38.over.chain.gz", "mpravardb": "mpravardb.bb",
         "sample": "GSM3130443_designed_library.csv.gz",
         "mqtl": "MuscleSkeletal.regular.perm.fdr.txt", "epic": "EPIC.hg38.manifest.tsv.gz"}
SOURCES = ("kircher", "sure", "mpravardb", "sample", "mqtl")
COLS = ["chrom", "pos", "ref", "alt", "source", "element", "cellule", "effet", "p", "signif"]
COMP = str.maketrans("ACGTN", "TGCAN")
ACTIVE = re.compile(r"(^|_)(TssA|TssAFlnk|Enh|EnhG|TxFlnk)$")
csv.field_size_limit(1 << 30)


def log(msg):
    print(msg, file=sys.stderr)


def opener(p):
    p = str(p)
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def download(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    log(f"téléchargement {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (mitoloop)"})
    with urllib.request.urlopen(req, timeout=300) as r, open(tmp, "wb") as fh:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            fh.write(b)
    tmp.rename(dest)
    return dest


def chrom_ok(c):
    c = str(c)
    c = c[3:] if c.lower().startswith("chr") else c
    c = "M" if c in ("MT", "M") else c
    return f"chr{c}" if re.fullmatch(r"\d{1,2}|X|Y|M", c) else None


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if v == v else None


class Fasta:
    """bases d'un FASTA indexé (.fai à côté), fichier non compressé"""

    def __init__(self, path):
        self.fh = open(path, "rb")
        self.idx = {}
        for line in open(str(path) + ".fai"):
            f = line.split("\t")
            self.idx[f[0]] = tuple(int(x) for x in f[1:5])

    def name(self, c):
        k = c[3:] if c.startswith("chr") else c
        for n in (c, k, "chr" + k, "MT" if k == "M" else None, "chrM" if k == "MT" else None):
            if n and n in self.idx:
                return n
        return None

    def seq(self, c, start, end):
        """bases start..end (1-based, incluses) ; '' hors du génome"""
        n = self.name(c)
        if n is None:
            return ""
        length, off, lb, lw = self.idx[n]
        if start < 1 or end > length or end < start:
            return ""
        b0 = off + (start - 1) // lb * lw + (start - 1) % lb
        b1 = off + (end - 1) // lb * lw + (end - 1) % lb
        self.fh.seek(b0)
        return self.fh.read(b1 - b0 + 1).decode().replace("\n", "").replace("\r", "").upper()


class Chain:
    """chaîne UCSC (hg19ToHg38.over.chain.gz) : conversion d'une base ; une
    base hors de tout bloc, ou couverte par les blocs de plusieurs chaînes,
    n'est pas convertie"""

    def __init__(self, path):
        blocks = defaultdict(list)       # chrom t -> (début t, fin t, chrom q, début q, brin q, taille q)
        tc = qc = qstrand = None
        tpos = qpos = qsize = 0
        with opener(path) as fh:
            for line in fh:
                f = line.split()
                if not f:
                    continue
                if f[0] == "chain":
                    tc, tpos = f[2], int(f[5])
                    qc, qsize, qstrand, qpos = f[7], int(f[8]), f[9], int(f[10])
                    continue
                size = int(f[0])
                blocks[tc].append((tpos, tpos + size, qc, qpos, qstrand, qsize))
                if len(f) == 3:
                    tpos += size + int(f[1])
                    qpos += size + int(f[2])
        self.blocks, self.starts, self.maxend = {}, {}, {}
        for c, v in blocks.items():
            v.sort()
            m, me = 0, []
            for b in v:
                m = max(m, b[1])
                me.append(m)
            self.blocks[chrom_ok(c) or c] = v
            self.starts[chrom_ok(c) or c] = [b[0] for b in v]
            self.maxend[chrom_ok(c) or c] = me

    def point(self, chrom, pos):
        """base 1-based -> (chrom, base 1-based, brin) ou None"""
        c = chrom_ok(chrom) or chrom
        v = self.blocks.get(c)
        if not v:
            return None
        p0 = pos - 1
        j = bisect_right(self.starts[c], p0) - 1
        hits = []
        while j >= 0 and self.maxend[c][j] > p0:
            ts, te, qc, qs, st, qsize = v[j]
            if ts <= p0 < te:
                q = qs + p0 - ts
                hits.append((chrom_ok(qc) or qc, (q if st == "+" else qsize - 1 - q) + 1, st))
            j -= 1
        return hits[0] if len(hits) == 1 else None


def lift(chain, chrom, pos, ref, alt):
    """variant hg19 -> hg38 (allèles complémentés sur une chaîne -)"""
    h = chain.point(chrom, pos)
    if h is None:
        return None
    if len(ref) > 1:
        h2 = chain.point(chrom, pos + len(ref) - 1)
        if h2 is None or h2[0] != h[0] or abs(h2[1] - h[1]) != len(ref) - 1:
            return None
    c, p, st = h
    if st == "-":
        ref, alt = ref.translate(COMP)[::-1], alt.translate(COMP)[::-1]
        p = p - len(ref) + 1
        if len(ref) != len(alt):          # indel sur le brin - : ancrage perdu
            return None
    return c, p, ref, alt


def row(chrom, pos, ref, alt, source, element=".", cellule=".", effet=None, p=None, signif="."):
    return {"chrom": chrom, "pos": int(pos), "ref": ref.upper(), "alt": alt.upper(),
            "source": source, "element": element or ".", "cellule": cellule or ".",
            "effet": "." if effet is None else f"{effet:.4g}",
            "p": "." if p is None else f"{p:.3g}", "signif": signif}


# ------------------------------------------------------------- sources

KIRCHER_CELLS = {"TERT-HEK": "HEK293T", "TERT-GBM": "SF7996", "TERT-GAa": "glioblastome",
                 "TERT-GSc": "glioblastome"}


def src_kircher(d, fasta):
    out, n_del = [], 0
    with opener(d / FILES["kircher"]) as fh:
        for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            if r.get("Release") != "GRCh38":
                continue
            p, eff = num(r.get("pValue")), num(r.get("Coefficient"))
            if p is None or eff is None:
                continue
            c, pos, ref, alt = chrom_ok(r["Chrom"]), int(r["Pos"]), r["Ref"], r["Alt"]
            if c is None:
                continue
            if alt == "-":                    # délétion d'une base : ancrée sur la base précédente
                b = fasta.seq(c, pos - 1, pos - 1)
                if not b:
                    continue
                pos, ref, alt = pos - 1, b + ref, b
                n_del += 1
            sig = "1" if p < 1e-5 else "0" if p > 0.05 else "."
            out.append(row(c, pos, ref, alt, "kircher2019", r["Element"],
                           KIRCHER_CELLS.get(r["Element"], "."), eff, p, sig))
    log(f"kircher : {len(out)} variants ({n_del} délétions)")
    return out


def sure_rows(path):
    with opener(path) as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def sure_effect(r, cell):
    a, b = num(r.get(f"{cell}.ref.mean")), num(r.get(f"{cell}.alt.mean"))
    return None if a is None or b is None else math.log2((b + 1) / (a + 1))


def src_sure(d, fasta, chain, ratio, seed):
    pos_rows, ids, act = [], defaultdict(set), {}
    for cell in ("k562", "hepg2"):
        f = d / FILES[f"sure_{cell}"]
        if not f.exists():
            log(f"sure : liste des raQTL {cell} absente")
            continue
        sig = list(sure_rows(f))
        ids[cell] = {r["SNP_ID"] for r in sig}
        mx = sorted(max(num(r[f"{cell}.ref.mean"]) or 0, num(r[f"{cell}.alt.mean"]) or 0) for r in sig)
        act[cell] = mx[len(mx) // 20] if mx else 0
        pos_rows += [(r, cell, "1") for r in sig]
    if not pos_rows:
        return []
    # raQTL : base 0 ou 1 de SNPabspos, choisie sur la base de référence hg38
    off = 0
    test = pos_rows[:3000]
    score = {}
    for o in (0, 1):
        ok = 0
        for r, _c, _s in test:
            x = lift(chain, r["chr"], int(r["SNPabspos"]) + o, r["ref"], r["alt"])
            ok += bool(x) and fasta.seq(x[0], x[1], x[1] + len(x[2]) - 1) == x[2]
        score[o] = ok
    off = max(score, key=score.get)
    log(f"sure : SNPabspos + {off} (base de référence retrouvée : {score[off]}/{len(test)})")
    neg = []
    tab = d / FILES["sure_table"]
    if tab.exists():
        rng = random.Random(seed)
        k = {c: ratio * len(v) for c, v in ids.items()}
        seen = defaultdict(int)
        res = {c: [] for c in ids}
        for r in sure_rows(tab):
            rc, ac = num(r.get("ref.element.count")), num(r.get("alt.element.count"))
            if rc is None or ac is None or not (10 <= rc <= 999 and 10 <= ac <= 999):
                continue
            for cell in ids:
                p = num(r.get(f"{cell}.wilcox.p.value"))
                if p is None or p <= 0.5 or r["SNP_ID"] in ids[cell]:
                    continue
                a, b = num(r.get(f"{cell}.ref.mean")) or 0, num(r.get(f"{cell}.alt.mean")) or 0
                if max(a, b) < act[cell]:
                    continue
                seen[cell] += 1                         # tirage par réservoir
                if len(res[cell]) < k[cell]:
                    res[cell].append(r)
                else:
                    j = rng.randrange(seen[cell])
                    if j < k[cell]:
                        res[cell][j] = r
        for cell, v in res.items():
            neg += [(r, cell, "0") for r in v]
            log(f"sure : {cell} {len(ids[cell])} raQTL, {len(v)} SNP sans effet tirés sur {seen[cell]}")
    else:
        log("sure : table complète absente : pas de SNP sans effet")
    out, lost = [], 0
    for r, cell, s in pos_rows + neg:
        x = lift(chain, r["chr"], int(r["SNPabspos"]) + off, r["ref"], r["alt"])
        if x is None:
            lost += 1
            continue
        out.append(row(*x, "sure2019", r["SNP_ID"], cell.upper().replace("HEPG2", "HepG2"),
                       sure_effect(r, cell), num(r.get(f"{cell}.wilcox.p.value")), s))
    log(f"sure : {len(out)} variants ({lost} non convertis en GRCh38)")
    return out


def bigbed_rows(path):
    """lignes du bigBed (ou d'un BED texte avec en-tête) : dictionnaires
    nommés par l'autoSql"""
    with open(path, "rb") as fh:
        magic = fh.read(4)
    if magic in (b"\xeb\xf2\x89\x87", b"\x87\x89\xf2\xeb"):
        import pyBigWig
        bb = pyBigWig.open(str(path))
        names = re.findall(r"^\s*\S+(?:\[\S*\])?\s+(\w+)\s*;", bb.SQL().decode()
                           if isinstance(bb.SQL(), bytes) else bb.SQL(), re.M)
        for c, size in bb.chroms().items():
            for s, e, rest in bb.entries(c, 0, size) or []:
                yield dict(zip(names, [c, str(s), str(e)] + rest.split("\t")))
    else:
        with opener(path) as fh:
            yield from csv.DictReader((ln.lstrip("#") for ln in fh), delimiter="\t")


def pick(r, *keys):
    low = {k.lower(): v for k, v in r.items()}
    for k in keys:
        for name, v in low.items():
            if name == k or (k.endswith("*") and name.startswith(k[:-1])):
                return v
    return None


def src_mpravardb(d):
    out = []
    for r in bigbed_rows(d / FILES["mpravardb"]):
        c = chrom_ok(pick(r, "chrom", "chr"))
        ref, alt = (pick(r, "ref") or "").upper(), (pick(r, "alt") or "").upper()
        if c is None or not re.fullmatch(r"[ACGT]+", ref) or not re.fullmatch(r"[ACGT]+", alt):
            continue
        pos = int(pick(r, "chromstart", "start")) + 1
        fdr, p = num(pick(r, "fdr")), num(pick(r, "pvalue", "pval*"))
        sig = "1" if fdr is not None and fdr < 0.05 else \
            "0" if (fdr is not None and fdr >= 0.5) or (fdr is None and p is not None and p > 0.5) \
            else "."
        out.append(row(c, pos, ref, alt, "mpravardb:" + (pick(r, "study", "mpra_study") or "?"),
                       pick(r, "rsid", "name") or ".", pick(r, "cellline", "cell*") or ".",
                       num(pick(r, "log2fc")), p if fdr is None else fdr, sig))
    log(f"mpravardb : {len(out)} variants")
    return out


def src_sample(d):
    ref_rl, snv = {}, []
    with opener(d / FILES["sample"]) as fh:
        for r in csv.DictReader(fh):
            lib, rl = r.get("library", ""), num(r.get("rl"))
            if rl is None:
                continue
            if lib == "snv":
                snv.append(r)
            elif lib.startswith("human"):
                ref_rl[r.get("utr", "")] = rl
    out, miss = [], 0
    for r in snv:
        m = re.match(r"NC_0+(\d+)\.\d+:g\.(\d+)([ACGT])>([ACGT])", r.get("id", ""))
        base = ref_rl.get(r.get("mother", ""))
        if not m or base is None:
            miss += 1
            continue
        c = {"23": "X", "24": "Y"}.get(m.group(1), m.group(1))
        out.append(row(f"chr{c}", m.group(2), m.group(3), m.group(4), "sample2019",
                       r.get("info1") or ".", "HEK293T", num(r["rl"]) - base))
    log(f"sample : {len(out)} SNV de 5'UTR ({miss} sans 5'UTR de référence ou sans position)")
    return out


def src_mqtl(d, max_dist, min_dist=50):
    pos = {}
    with opener(d / FILES["epic"]) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            c = chrom_ok(r.get("CpG_chrm", ""))
            if c and r.get("CpG_beg", "NA") != "NA":
                pos[r["Probe_ID"]] = (c, int(r["CpG_beg"]) + 1)     # C du CpG, 1-based
    n_far = 0
    rows = []
    with opener(d / FILES["mqtl"]) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            q, cg = num(r.get("qval")), pos.get(r["cpg_id"].split(":")[0])
            m = re.fullmatch(r"(chr\w+)_(\d+)_([ACGT]+)_([ACGT]+)_b38", r.get("variant_id", ""))
            if q is None or cg is None or not m or m.group(1) != cg[0]:
                continue
            dist = int(m.group(2)) - cg[1]
            if not (min_dist < abs(dist) <= max_dist):
                n_far += 1
                continue
            sig = "1" if q < 0.05 else "0" if q > 0.5 else "."
            rows.append(row(m.group(1), m.group(2), m.group(3), m.group(4), "mqtl_muscle_gtex",
                            f"{r['cpg_id']}({dist:+d})", "muscle", num(r.get("slope")), q, sig))
    log(f"mqtl : {len(rows)} paires variant-CpG à {min_dist}-{max_dist} pb ({n_far} hors distance ; "
        f"{len(pos)} sondes EPIC)")
    return rows


def src_table(path, build, chain):
    out, lost = [], 0
    with opener(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r = {k.lower(): v for k, v in r.items()}
            c = chrom_ok(r.get("chrom", ""))
            if c is None or not r.get("pos", "").isdigit():
                continue
            v = (c, int(r["pos"]), r["ref"].upper(), r["alt"].upper())
            if build == "hg19":
                v = lift(chain, *v)
                if v is None:
                    lost += 1
                    continue
            s = r.get("signif", ".")
            out.append(row(*v, Path(path).name.split(".")[0], r.get("element", "."),
                           r.get("cellule", "."), num(r.get("effet")), num(r.get("p")),
                           s if s in ("0", "1") else "."))
    log(f"{Path(path).name} : {len(out)} variants ({lost} non convertis)")
    return out


def cmd_build(a):
    d = a.datadir
    d.mkdir(parents=True, exist_ok=True)
    fasta = Fasta(a.fasta)
    want = set(a.sources)
    need = {"kircher": ["kircher"], "sure": ["sure_k562", "sure_hepg2", "sure_table", "chain"],
            "mpravardb": ["mpravardb"], "sample": ["sample"], "mqtl": ["mqtl", "epic"]}
    if a.table_hg19:
        need.setdefault("table", ["chain"])
        want.add("table")
    for s in sorted(want):
        for k in need.get(s, []):
            try:
                download(URLS[k], d / FILES[k])
            except Exception as e:  # source indisponible : sautée
                log(f"{k} : téléchargement impossible ({type(e).__name__}: {e})")
    chain = Chain(d / FILES["chain"]) if (d / FILES["chain"]).exists() else None
    rows = []
    for s, f in (("kircher", lambda: src_kircher(d, fasta)),
                 ("sure", lambda: src_sure(d, fasta, chain, a.ratio_negatifs, a.seed)),
                 ("mpravardb", lambda: src_mpravardb(d)), ("sample", lambda: src_sample(d)),
                 ("mqtl", lambda: src_mqtl(d, a.mqtl_distance))):
        if s not in want:
            continue
        if s == "sure" and chain is None:
            log("sure : chaîne hg19ToHg38 absente : source sautée")
            continue
        try:
            rows += f()
        except Exception as e:  # format inattendu : source sautée
            log(f"{s} : illisible ({type(e).__name__}: {e})")
    for t in a.table:
        rows += src_table(t, "hg38", chain)
    for t in a.table_hg19:
        rows += src_table(t, "hg19", chain)
    # base de référence contrôlée sur GRCh38 ; doublons écartés
    ok, bad, seen = [], defaultdict(int), set()
    for r in rows:
        k = (r["chrom"], r["pos"], r["ref"], r["alt"], r["source"], r["element"], r["cellule"])
        if k in seen:
            continue
        seen.add(k)
        if fasta.seq(r["chrom"], r["pos"], r["pos"] + len(r["ref"]) - 1) != r["ref"]:
            bad[r["source"]] += 1
            continue
        ok.append(r)
    for s, n in sorted(bad.items()):
        log(f"{s} : {n} variants écartés (base de référence discordante sur GRCh38)")
    # --max : positions tirées au hasard dans chaque source
    rng = random.Random(a.seed)
    by_src = defaultdict(lambda: defaultdict(list))
    for r in ok:
        by_src[r["source"]][(r["chrom"], r["pos"])].append(r)
    final = []
    for s, groups in sorted(by_src.items()):
        keys = sorted(groups)
        rng.shuffle(keys)
        n = 0
        for k in keys:
            if a.max and n >= a.max:
                break
            final += groups[k]
            n += len(groups[k])
        sig = defaultdict(int)
        for k in keys:
            for r in groups[k]:
                sig[r["signif"]] += 1
        log(f"{s} : {sum(sig.values())} variants (significatifs {sig['1']}, sans effet {sig['0']}, "
            f"intermédiaires ou sans test {sig['.']})" + (f" ; {n} gardés" if a.max and n < sum(sig.values()) else ""))
    final.sort(key=lambda r: (r["chrom"], r["pos"], r["alt"], r["source"]))
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in final:
            fh.write("\t".join(str(r[c]) for c in COLS) + "\n")
    log(f"{a.out} : {len(final)} variants")


# ------------------------------------------------------------- regions, entree

def read_variants(path):
    with opener(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def cmd_regions(a):
    a.out.mkdir(parents=True, exist_ok=True)
    win = defaultdict(list)
    for r in read_variants(a.variants):
        p = int(r["pos"])
        win[r["chrom"]].append((p - 1 - a.flank, p - 1 + len(r["ref"]) + a.flank))
    n = 0
    with open(a.out / "regions_eval.bed", "w") as fh:
        fh.write("#chrom\tstart\tend\tgene\tclass\tsource\ttss\tstrand\tdetail\n")
        for c in sorted(win):
            cur = None
            for s, e in sorted(win[c]):
                s = max(0, s)
                if cur and s <= cur[1]:
                    cur[1] = max(cur[1], e)
                    continue
                if cur:
                    fh.write(f"{c}\t{cur[0]}\t{cur[1]}\tMPRA\telement\tmpra\t.\t.\t\n")
                    n += 1
                cur = [s, e]
            if cur:
                fh.write(f"{c}\t{cur[0]}\t{cur[1]}\tMPRA\telement\tmpra\t.\t.\t\n")
                n += 1
    (a.out / "genes_eval.txt").write_text("MPRA\n")
    (a.out / "l1_vide.tsv").write_text("gene\tgene_id\tset\tstatus\tscore\tci_low\tci_high\t"
                                       "complex\tprotein_coupling\trna_coupling\tl1_class\n")
    log(f"{a.out}/regions_eval.bed : {n} régions")


def cmd_entree(a):
    regs = defaultdict(list)
    with opener(a.regions) as fh:
        for r in csv.DictReader((ln.lstrip("#") for ln in fh), delimiter="\t"):
            if r.get("class") in ("cds", "exon_nc", "gene_body"):
                continue
            c = chrom_ok(r["chrom"]) or r["chrom"]
            regs[c].append((int(r["start"]), int(r["end"]), r["gene"], r["class"]))
    idx = {}
    for c, v in regs.items():
        v.sort()
        idx[c] = ([x[0] for x in v], v, max(x[1] - x[0] for x in v))
    n = 0
    with open(a.out, "w") as fh:
        fh.write("\t".join(["transcrit", "variant_entree", "gene", "classe"] + COLS) + "\n")
        for r in read_variants(a.variants):
            x = idx.get(r["chrom"])
            if not x:
                continue
            st, v, ml = x
            p0 = int(r["pos"]) - 1
            hit = [h for h in v[bisect_left(st, p0 - ml):bisect_right(st, p0)] if p0 < h[1]]
            if not hit:
                continue
            genes = ",".join(sorted({h[2] for h in hit}))
            classes = ",".join(sorted({h[3] for h in hit}))
            fh.write("\t".join([".", ".", genes, classes] + [r[c] for c in COLS]) + "\n")
            n += 1
    log(f"{a.out} : {n} variants dans les régions MitoLoop")


# ------------------------------------------------------------- evaluer

def ranks(x):
    order = sorted(range(len(x)), key=lambda i: x[i])
    r = [0.0] * len(x)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def auc(pos, neg):
    if not pos or not neg:
        return None
    r = ranks(pos + neg)
    return (sum(r[:len(pos)]) - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def spearman(x, y):
    if len(x) < 10:
        return None
    rx, ry = ranks(x), ranks(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    sx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    sy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return sxy / (sx * sy) if sx and sy else None


def features(r):
    """annotations d'un variant annoté : {nom: (classe, valeur)} ; classe B =
    position, C = allèle"""
    f = {}
    x = num(r.get("MITO_phyloP"))
    if x is not None:
        f["phyloP"] = ("B", x)
    if r.get("MITO_CAGE") in ("oui", "non"):
        f["CAGE"] = ("B", float(r["MITO_CAGE"] == "oui"))
    fix = r.get("MITO_fixation", ".")
    if fix is not None:
        f["n_facteurs_fixes"] = ("B", float(0 if fix in (".", "") else len(fix.split(","))))
    ch = r.get("MITO_chromatine", ".")
    if ch is not None:
        f["chromatine_active"] = ("B", float(any(ACTIVE.search(s.split(":")[-1])
                                                  for s in ch.split(",") if s != ".")))
    mot = r.get("MITO_motifs")
    if mot is not None and mot != "REF_discordante":
        ds, ics = [], []
        for m in re.finditer(r"/IC([-\d.]+)/d([-+\d.]+)", mot):
            ics.append(float(m.group(1)))
            ds.append(float(m.group(2)))
        f["perte_site"] = ("C", -min(ds + [0.0]))
        f["delta_site_abs"] = ("C", max([abs(d) for d in ds] + [0.0]))
        f["delta_site_signe"] = ("C", max(ds, key=abs) if ds else 0.0)   # le site le plus changé
        f["site_dans_motif"] = ("B", float(bool(ds)))
        f["IC_position"] = ("B", max(ics + [0.0]))
    args = r.get("MITO_arguments")
    if args is not None:
        f["site_TF"] = ("C", float("site_" in args))
        f["n_arguments_position"] = ("B", float(sum(1 for x in args.split(",")
                                                    if x != "." and not x.startswith("site_")
                                                    and x != "L1_transmis")))
    cpg = r.get("MITO_CpG")
    if cpg is not None and cpg != "REF_discordante":
        f["CpG_detruit_ou_cree"] = ("C", float(cpg.startswith(("CpG_detruit", "CpG_cree"))))
        m = re.search(r"oe=([\d.]+)>([\d.]+)", cpg)
        if m:
            f["CpG_oe_fenetre"] = ("B", float(m.group(1)))
            f["CpG_oe_perte"] = ("C", float(m.group(1)) - float(m.group(2)))
    if args is not None:
        f["CpG_site_TF"] = ("C", float("CpG_site_" in args))
    ev = r.get("MITO_evenements")
    if ev is not None:                 # événements propres à MitoLoop : fort 2, modéré 1
        own = [e for e in ev.split(",") if e.startswith(("site_", "CpG_site_", "bascule_"))]
        f["evenements_MitoLoop"] = ("C", float(sum(2 if e.endswith("[fort]") else 1 for e in own)))
    # SpliceAI, SpliceAI2 (colonnes SpliceAI*) et AlphaGenome (AG_*) ne sont pas
    # évalués ici : selon leurs conditions d'utilisation, leurs sorties ne
    # servent ni à entraîner ni à ajuster des seuils
    for k, v in r.items():
        if k and k.startswith("MORFEE_") and k.endswith("creation"):
            f[k] = ("C", float(v not in (None, "", ".", "NA")))
    return f


def cmd_evaluer(a):
    meas = defaultdict(list)
    for r in read_variants(a.variants):
        meas[(r["chrom"], r["pos"], r["ref"], r["alt"])].append(r)
    feats = {}
    for path in a.annote:
        with opener(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
                r = {k: v for k, v in r.items() if k}
                low = {k.lower().lstrip("#"): v for k, v in r.items()}
                if low.get("mito_statut", "annote") != "annote" or not low.get("pos", "").isdigit():
                    continue
                k = (chrom_ok(low["chrom"]) or low["chrom"], low["pos"], low["ref"].upper(), low["alt"].upper())
                if k in meas:
                    feats.setdefault(k, {}).update(features(r))
    log(f"{len(feats)} variants annotés retrouvés sur {len(meas)}")
    # groupes : source x lignée, et la source entière
    groups = defaultdict(list)
    for k, ms in meas.items():
        if k not in feats:
            continue
        for m in ms:
            src = m["source"]
            groups[(src, m["cellule"])].append((k, m))
            groups[(src, "toutes")].append((k, m))
    names = sorted({n for f in feats.values() for n in f})
    out = []
    for (src, cell), items in sorted(groups.items()):
        if cell != "toutes" and len({c for s, c in groups if s == src}) <= 2:
            continue                       # une seule lignée : la ligne « toutes » suffit
        for n in names:
            cls = next(f[n][0] for f in feats.values() if n in f)
            pos_v, neg_v, eff, val = [], [], [], []
            by_pos = defaultdict(lambda: ([], []))
            for k, m in items:
                if n not in feats[k]:
                    continue
                v = feats[k][n][1]
                if m["signif"] == "1":
                    pos_v.append(v)
                    by_pos[(k[0], k[1], m["element"], m["cellule"])][0].append(v)
                elif m["signif"] == "0":
                    neg_v.append(v)
                    by_pos[(k[0], k[1], m["element"], m["cellule"])][1].append(v)
                e = num(m["effet"])
                if e is not None:
                    eff.append(e)
                    val.append(v)
            # paires d'allèles d'une même base : significatif contre sans effet
            pairs = wins = 0.0
            for pv, nv in by_pos.values():
                for x in pv:
                    for y in nv:
                        pairs += 1
                        wins += 1 if x > y else 0.5 if x == y else 0
            A = auc(pos_v, neg_v)
            sa = spearman(val, [abs(e) for e in eff])
            ss = spearman(val, eff) if n == "delta_site_signe" else None
            out.append([src, cell, n, cls, len(pos_v), len(neg_v),
                        "." if A is None else f"{A:.3f}", int(pairs),
                        f"{wins / pairs:.3f}" if pairs else ".", len(eff),
                        "." if sa is None else f"{sa:.3f}", "." if ss is None else f"{ss:.3f}"])
    head = ["source", "cellule", "annotation", "classe", "n_significatifs", "n_sans_effet",
            "auc", "n_paires_meme_position", "auc_meme_position", "n_effet",
            "spearman_effet_abs", "spearman_effet_signe"]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("\t".join(head) + "\n")
        for r in out:
            fh.write("\t".join(str(x) for x in r) + "\n")
    log(f"{a.out} : {len(out)} lignes")
    for r in out:
        if r[1] == "toutes" and r[2] in ("perte_site", "site_TF", "phyloP", "evenements_MitoLoop"):
            log(f"  {r[0]:28s} {r[2]:18s} AUC {r[6]:>6s}  même position {r[8]:>6s} "
                f"({r[7]} paires)  Spearman |effet| {r[10]}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("--datadir", type=Path, default=Path("data/mpra"))
    b.add_argument("--fasta", required=True, help="GRCh38 indexé (.fai), non compressé")
    b.add_argument("--out", type=Path, required=True)
    b.add_argument("--sources", nargs="*", default=list(SOURCES), choices=SOURCES)
    b.add_argument("--table", type=Path, action="append", default=[], help="table en GRCh38")
    b.add_argument("--table-hg19", type=Path, action="append", default=[], help="table en hg19")
    b.add_argument("--ratio-negatifs", type=int, default=3)
    b.add_argument("--mqtl-distance", type=int, default=1000,
                   help="distance maximale variant-CpG des mQTL gardés (pb)")
    b.add_argument("--max", type=int, default=0, help="variants par source au plus (0 : tous)")
    b.add_argument("--seed", type=int, default=1)
    r = sub.add_parser("regions")
    r.add_argument("--variants", type=Path, required=True)
    r.add_argument("--out", type=Path, required=True)
    r.add_argument("--flank", type=int, default=50)
    e = sub.add_parser("entree")
    e.add_argument("--variants", type=Path, required=True)
    e.add_argument("--regions", type=Path, required=True)
    e.add_argument("--out", type=Path, required=True)
    v = sub.add_parser("evaluer")
    v.add_argument("--variants", type=Path, required=True)
    v.add_argument("--annote", nargs="+", required=True)
    v.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    {"build": cmd_build, "regions": cmd_regions, "entree": cmd_entree,
     "evaluer": cmd_evaluer}[a.cmd](a)


if __name__ == "__main__":
    main()
