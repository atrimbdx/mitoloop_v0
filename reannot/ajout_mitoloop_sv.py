#!/usr/bin/env python3
"""
ajout_mitoloop_sv.py — VM. Arguments MitoLoop pour les variants structuraux
(délétions, duplications, inversions, insertions, points de cassure) : ce
que le SV retire, copie ou sépare d'un gène mitochondrial et de ses régions
régulatrices. Le fichier d'origine n'est pas modifié : le résultat est écrit
à côté, avec le suffixe .mitoloop.tsv.

Usage : python3 ajout_mitoloop_sv.py PREFIXE.sv.vcf [--snv run.final.mitoloop.tsv]
        python3 ajout_mitoloop_sv.py PREFIXE.annotsv.tsv     (sortie d'AnnotSV)
        (REANNOT_BASE, MITOLOOP_DIR comme ajout_mitoloop.py ; --data, --gtf
        pour les donner autrement)

Entrées acceptées
  VCF de SV (SVTYPE et END dans INFO, ou allèle symbolique <DEL>, <DUP>...),
  celui de preparer_entree.py --sv ; sortie TSV d'AnnotSV (SV_chrom,
  SV_start, SV_end, SV_type ; seules les lignes Annotation_mode = full) ; ou
  toute table avec en-tête chrom, start, end, type. Positions en GRCh38,
  début et fin en base 1 incluses.

Données : le paquet MitoLoop (genes.tsv, regions.bed.gz, annotations/
tss_types.bed.gz, annotations/gnomad_sv.bed.gz s'il est présent) et le GTF
Ensembl du pipeline (structure des gènes : transcrit MANE Select, sinon
Ensembl canonique, sinon le plus long).

Colonnes ajoutées
  MITO_SV_genes        gènes mitochondriaux touchés et ce qui l'est : exons
                       (k/n du transcrit), TSS, gène entier, régions
                       (promoter, utr5, utr3, element, anchor), points de
                       cassure dans le gène
  MITO_SV_effet        effet de chaque gène (voir Événements)
  MITO_L1_classe, MITO_type_promoteur, MITO_gnomAD, MITO_GISMO   par gène
                       (GENE:valeur, séparés par |)
  MITO_SV_gnomAD       SV de gnomAD-SV v4.1 le plus semblable : même type,
                       recouvrement réciproque >= 50 % ; fréquence, type,
                       recouvrement
  MITO_second_allele   avec --snv : variants de la table SNV (même patient)
                       dans le même gène, et leur priorité MitoLoop ; en
                       récessif, un SV est souvent le second allèle
  MITO_evenements, MITO_priorite, MITO_statut (annote, hors_genes_mito, NA)

Événements
  fort     DEL du gène entier, d'exons ou du TSS ; DEL du promoteur ou du
           5'UTR ; DUP interne au gène (deux points de cassure dans le gène,
           au moins un exon copié : cadre de lecture rompu) ; INV ou point de
           cassure (BND) qui coupe le gène ou son promoteur ; INS dans un
           exon ou le promoteur
  modéré   DEL d'un élément distal (rE2G), d'une ancre de boucle ou du 3'UTR ;
           DUP partielle (un seul point de cassure dans le gène) ; INV qui
           sépare le gène de ses éléments ou ancres (gène d'un côté du point
           de cassure, éléments de l'autre) ; BND près du gène (régions)
           Comme pour les SNV, un événement modéré d'un élément distal compte
           comme fort si le promoteur du gène est de type tissu ou
           développement.
  noté     DUP du gène entier (dosage, peu pertinent en récessif), DEL ou DUP
           seulement intronique ou d'un élément d'un gène de ménage : pas
           d'événement, mais 3_contexte
Priorité : 1_fort (événement fort, classe L1 autre que not_transmitted),
2_modere, 3_contexte (le SV touche le gène ou ses régions sans événement),
4_aucun ; un SV fréquent dans gnomAD-SV (AF >= --af-max, 1 %) descend d'un
niveau, avec la fréquence dans MITO_evenements.

Python 3 de la bibliothèque standard.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ajout_mitoloop import (Intervals, base_dir, chrom_key, gene_constraints,  # noqa: E402
                            load_genes, opener)

TYPES = {"DEL": "DEL", "DELETION": "DEL", "LOSS": "DEL", "DUP": "DUP", "DUPLICATION": "DUP",
         "GAIN": "DUP", "TANDEMDUP": "DUP", "INV": "INV", "INVERSION": "INV", "INS": "INS",
         "INSERTION": "INS", "BND": "BND", "TRA": "BND", "CTX": "BND", "CNV": "CNV", "CPX": "CPX"}
PROX = {"promoter", "utr5"}
DISTAL = {"element", "anchor"}
ORDER = ["1_fort", "2_modere", "3_contexte", "4_aucun"]


def sv_type(x):
    x = (x or "").strip().upper().strip("<>")
    x = re.split(r"[:_]", x)[0] if x not in TYPES else x
    return TYPES.get(x, x or "?")


def read_svs(path):
    """lignes (en-tête, champs, chrom, début 1-based, fin, type) d'un VCF, d'une
    sortie AnnotSV ou d'une table chrom/start/end/type"""
    if re.search(r"\.vcf(\.gz)?$", str(path)):
        head, rows = ["chrom", "pos", "id", "ref", "alt", "end", "svtype"], []
        with opener(path) as fh:
            for ln in fh:
                if ln.startswith("#"):
                    continue
                f = ln.rstrip("\n").split("\t")
                info = dict(x.split("=", 1) if "=" in x else (x, "") for x in f[7].split(";"))
                t = sv_type(info.get("SVTYPE") or f[4])
                end = int(info["END"]) if info.get("END", "").lstrip("-").isdigit() else \
                    int(f[1]) + max(len(f[3]) - 1, 0)
                rows.append(([f[0], f[1], f[2], f[3], f[4], str(end), t], f[0], int(f[1]), end, t))
        return head, rows
    with opener(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    sep = "\t" if lines[0].count("\t") >= lines[0].count(",") else ","
    head = next(csv.reader([lines[0]], delimiter=sep))
    low = {h.lower().lstrip("#"): i for i, h in enumerate(head)}

    def col(*names):
        return next((low[n] for n in names if n in low), None)
    C, S, E, T = (col("sv_chrom", "chrom", "chr", "chromosome"),
                  col("sv_start", "start", "debut", "pos"), col("sv_end", "end", "fin"),
                  col("sv_type", "type", "svtype", "type_sv"))
    M = col("annotation_mode")
    if None in (C, S, E, T):
        sys.exit(f"{path} : colonnes chrom, start, end, type introuvables ({head[:10]})")
    rows = []
    for f in csv.reader(lines[1:], delimiter=sep):
        if M is not None and f[M] != "full":
            continue
        try:
            rows.append((f, f[C], int(float(f[S])), int(float(f[E])), sv_type(f[T])))
        except (ValueError, IndexError):
            rows.append((f, None, None, None, None))
    return head, rows


def gtf_models(gtf, wanted):
    """gène -> (chrom, brin, début tx 0-based, fin, TSS 0-based, exons [(s, e)])
    du transcrit MANE Select, sinon Ensembl canonique, sinon le plus long"""
    tx = {}
    with opener(gtf) as fh:
        for ln in fh:
            if ln[0] == "#":
                continue
            f = ln.split("\t", 9)
            if len(f) < 9 or f[2] not in ("transcript", "exon"):
                continue
            g = re.search(r'gene_name "([^"]+)"', f[8])
            if not g or g.group(1) not in wanted:
                continue
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            s, e = int(f[3]) - 1, int(f[4])
            if f[2] == "transcript":
                rank = 2 if "MANE_Select" in f[8] else 1 if "Ensembl_canonical" in f[8] else 0
                tx[t] = [g.group(1), f[0], f[6], s, e, rank, []]
            elif t in tx:
                tx[t][6].append((s, e))
    best = {}
    for g, c, st, s, e, rank, ex in tx.values():
        k = (rank, e - s)
        if g not in best or k > best[g][0]:
            best[g] = (k, (c, st, s, e, s if st == "+" else e - 1, sorted(ex)))
    return {g: v[1] for g, v in best.items()}


def recip(a, b, c, d):
    o = min(b, d) - max(a, c)
    return o / max(b - a, d - c, 1) if o > 0 else 0.0


def gnomad_match(gn, chrom, s0, e0, t):
    if gn is None or t not in ("DEL", "DUP", "INV", "CNV", "INS"):
        return None
    fam = {"DEL": {"DEL", "CNV"}, "DUP": {"DUP", "CNV"}, "CNV": {"DEL", "DUP", "CNV"},
           "INV": {"INV", "CPX"}, "INS": {"INS"}}[t]
    best = None
    for s, e, lab in gn.query(chrom, s0 - 1, max(e0, s0 + 1) + 1):
        if sv_type(lab[0]) not in fam:
            continue
        ro = 1.0 if t == "INS" and abs(s - s0) <= 100 else recip(s0, e0, s, e)
        if ro >= 0.5 and (best is None or ro > best[0]):
            best = (ro, lab)
    if not best:
        return None
    try:
        af = float(best[1][1])
    except (ValueError, IndexError):
        af = None
    return af, best[1][0], best[0]


def gene_effects(m, regs, t, s0, e0):
    """événements (force, texte) et parties touchées d'un gène ; s0, e0 :
    intervalle 0-based demi-ouvert du SV"""
    c, strand, ts, te, tss, ex = m if m else (None, "+", None, None, None, [])
    parts, ev = [], []
    bps = [s0] if t in ("INS", "BND") else [s0, e0 - 1]
    hit_regs = {k for rs, re_, k in regs if rs < e0 and s0 < re_} if t not in ("INS", "BND") \
        else {k for rs, re_, k in regs if rs <= s0 < re_}
    bp_regs = {k for rs, re_, k in regs for b in bps if rs <= b < re_}
    in_tx = [b for b in bps if ts is not None and ts < b < te - 1]
    ex_hit = [x for x in ex if x[0] < e0 and s0 < x[1]] if t not in ("INS", "BND") else \
        [x for x in ex if x[0] <= s0 < x[1]]
    whole = ts is not None and s0 <= ts and te <= e0 and t not in ("INS", "BND")
    tss_in = tss is not None and s0 <= tss < e0 and t not in ("INS", "BND")
    if whole:
        parts.append("gene_entier")
    elif ex_hit:
        parts.append(f"exons {len(ex_hit)}/{len(ex)}")
    if tss_in and not whole:
        parts.append("TSS")
    if in_tx:
        parts.append("cassure_dans_le_gene")
    parts += sorted(hit_regs)
    if t in ("DEL", "CNV"):
        if whole or ex_hit or tss_in:
            ev.append(("fort", "DEL:" + ("gene_entier" if whole else
                                         ",".join(p for p in parts if p.startswith(("exons", "TSS"))))))
        elif hit_regs & PROX:
            ev.append(("fort", "DEL:promoteur"))
        elif hit_regs & (DISTAL | {"utr3"}):
            ev.append(("modere", "DEL:" + ",".join(sorted(hit_regs & (DISTAL | {"utr3"})))))
    elif t == "DUP":
        if len(in_tx) == 2 and ex_hit:
            ev.append(("fort", f"DUP:interne exons {len(ex_hit)}/{len(ex)}"))
        elif len(in_tx) == 1:
            ev.append(("modere", "DUP:partielle"))
    elif t in ("INV", "CPX"):
        if in_tx or bp_regs & PROX:
            ev.append(("fort", "INV:coupe_" + ("le_gene" if in_tx else "le_promoteur")))
        elif tss is not None:
            inside = lambda p: s0 <= p < e0                     # noqa: E731
            sep = sorted({k for rs, re_, k in regs if k in DISTAL
                          and inside((rs + re_) // 2) != inside(tss)})
            if sep:
                ev.append(("modere", "INV:separe_" + ",".join(sep)))
    elif t == "BND":
        if in_tx or bp_regs & PROX:
            ev.append(("fort", "BND:coupe_" + ("le_gene" if in_tx else "le_promoteur")))
        elif bp_regs:
            ev.append(("modere", "BND:" + ",".join(sorted(bp_regs))))
    elif t == "INS":
        if ex_hit or bp_regs & PROX:
            ev.append(("fort", "INS:" + ("exon" if ex_hit else "promoteur")))
        elif bp_regs:
            ev.append(("modere", "INS:" + ",".join(sorted(bp_regs))))
    return parts, ev


def snv_index(paths):
    """gène -> variants de la table SNV et leur priorité MitoLoop"""
    out = defaultdict(list)
    for p in paths or []:
        with opener(p) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                low = {k.lower(): v for k, v in r.items() if k}
                genes = {low.get("gene", "")} | {x.split(":")[0] for x in
                                                 (low.get("mito_gene") or "").split("|")}
                tag = f"{low.get('chrom', '?')}:{low.get('pos', '?')}{low.get('ref', '')}>{low.get('alt', '')}"
                pr = low.get("mito_priorite", "")
                for g in genes - {"", ".", "NA"}:
                    out[g].append(tag + (f"[{pr}]" if pr and pr not in (".", "NA") else ""))
    return out


def main():
    p = argparse.ArgumentParser(description="Arguments MitoLoop pour les variants structuraux")
    p.add_argument("tables", nargs="+")
    b = base_dir()
    p.add_argument("--data", type=Path, default=Path(
        os.environ.get("MITOLOOP_DIR", b / "data" / "mitoloop")))
    rel = (b / "VEP_RELEASE").read_text().strip() if (b / "VEP_RELEASE").exists() else "116"
    p.add_argument("--gtf", type=Path, default=Path(
        os.environ.get("GTF", b / f"data/ref/Homo_sapiens.GRCh38.{rel}.gtf.gz")))
    p.add_argument("--snv", action="append", default=[],
                   help="table SNV annotée du même patient (second allèle), répétable")
    p.add_argument("--af-max", type=float, default=0.01)
    p.add_argument("--flank", type=int, default=1_000_000,
                   help="distance au gène au-delà de laquelle un SV n'est pas considéré")
    a = p.parse_args()
    D = a.data
    genes = load_genes(D / "genes.tsv")
    regs = defaultdict(list)                         # gène -> (début, fin, classe)
    with opener(D / "regions.bed.gz") as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) >= 5 and not ln.startswith("#"):
                regs[f[3]].append((chrom_key(f[0]), int(f[1]), int(f[2]), f[4]))
    models = gtf_models(a.gtf, set(genes)) if a.gtf.exists() else {}
    if not models:
        print(f"MitoLoop SV : GTF {a.gtf} absent ou sans gène du paquet : régions seules",
              file=sys.stderr)
    types = {}
    if (D / "annotations" / "tss_types.bed.gz").exists():
        for s, _e, lab in [(int(f[1]), 0, f[3:]) for f in (ln.rstrip("\n").split("\t") for ln in
                           opener(D / "annotations" / "tss_types.bed.gz"))]:
            types.setdefault(lab[0], lab[2])
    gn = Intervals(D / "annotations" / "gnomad_sv.bed.gz") \
        if (D / "annotations" / "gnomad_sv.bed.gz").exists() else None
    snv = snv_index(a.snv)
    # domaine de chaque gène : transcrit et régions
    dom = {}
    for g in genes:
        m = models.get(g)
        xs = [(c, s, e) for c, s, e, _k in regs.get(g, [])]
        if m:
            xs.append((chrom_key(m[0]), m[2], m[3]))
        if xs:
            dom[g] = (xs[0][0], min(x[1] for x in xs), max(x[2] for x in xs))
    print(f"MitoLoop SV : {len(genes)} gènes ({len(models)} avec leur structure), "
          f"gnomAD-SV {'présent' if gn else 'absent'}", file=sys.stderr)

    extra = ["MITO_SV_genes", "MITO_SV_effet", "MITO_L1_classe", "MITO_type_promoteur",
             "MITO_gnomAD", "MITO_GISMO", "MITO_SV_gnomAD", "MITO_second_allele",
             "MITO_evenements", "MITO_priorite", "MITO_statut"]
    for path in a.tables:
        head, rows = read_svs(path)
        out = re.sub(r"\.(tsv|csv|txt|vcf(\.gz)?)$", "", path) + ".mitoloop.tsv"
        n = defaultdict(int)
        with open(out, "w") as fo:
            fo.write("\t".join(head + extra) + "\n")
            for f, chrom, start, end, t in rows:
                if chrom is None:
                    fo.write("\t".join(f + ["."] * (len(extra) - 1) + ["NA"]) + "\n")
                    n["NA"] += 1
                    continue
                c = chrom_key(chrom)
                s0, e0 = start - 1, max(end, start)
                touched = []
                for g, (dc, ds, de) in dom.items():
                    if dc != c:
                        continue
                    near = ds - a.flank < e0 and s0 < de + a.flank
                    if not near or not (s0 < de and ds < e0 or t in ("INV", "CPX", "BND")
                                        and (ds <= s0 < de or ds <= e0 - 1 < de)):
                        continue
                    m = models.get(g)
                    rg = [(s, e, k) for rc, s, e, k in regs.get(g, []) if rc == c]
                    parts, ev = gene_effects(m, rg, t, s0, e0)
                    if parts or ev:
                        touched.append((g, parts, ev))
                if not touched:
                    fo.write("\t".join(f + ["."] * (len(extra) - 1) + ["hors_genes_mito"]) + "\n")
                    n["hors_genes_mito"] += 1
                    continue
                forts, moderes, best = [], [], "3_contexte"
                for g, parts, ev in touched:
                    gi = genes.get(g, {})
                    ty = types.get(g) or gi.get("type_promoteur", "")
                    up = ty.split("+")[0] in ("tissu", "developpement") and "menage" not in ty
                    for s, txt in ev:
                        if s == "modere" and up and any(k in txt for k in ("element", "anchor")):
                            s = "fort"
                        tier = "1_fort" if s == "fort" and gi.get("l1_class") != "not_transmitted" \
                            else "2_modere"
                        (forts if s == "fort" else moderes).append(f"{g}:{txt}")
                        best = min(best, tier, key=ORDER.index)
                gm = gnomad_match(gn, chrom, s0, e0, t)
                gtxt = "." if gm is None else \
                    f"AF={'NA' if gm[0] is None else f'{gm[0]:.3g}'}({gm[1]},RO={gm[2]:.2f})"
                evs = [f"{x}[fort]" for x in forts] + [f"{x}[modere]" for x in moderes]
                if gm and gm[0] is not None and gm[0] >= a.af_max:
                    best = ORDER[min(ORDER.index(best) + 1, 3)]
                    evs.append(f"gnomAD-SV_frequent(AF={gm[0]:.3g})")
                gs = [g for g, _p, _e in touched]
                second = [f"{g}:" + ",".join(snv[g][:5]) for g in gs if snv.get(g)]
                cons = [gene_constraints(genes.get(g, {})) for g in gs]
                v = {"MITO_SV_genes": "|".join(f"{g}:{','.join(p) or 'regions'}"
                                               for g, p, _e in touched),
                     "MITO_SV_effet": "|".join(f"{g}:{','.join(t_ for _s, t_ in e)}"
                                               for g, _p, e in touched if e) or ".",
                     "MITO_L1_classe": "|".join(f"{g}:{genes.get(g, {}).get('l1_class') or 'not_estimable'}"
                                                for g in gs),
                     "MITO_type_promoteur": "|".join(f"{g}:{types.get(g) or genes.get(g, {}).get('type_promoteur') or '.'}"
                                                     for g in gs),
                     "MITO_gnomAD": "|".join(f"{g}:{x[0]}" for g, x in zip(gs, cons)),
                     "MITO_GISMO": "|".join(f"{g}:{x[1]}" for g, x in zip(gs, cons)),
                     "MITO_SV_gnomAD": gtxt,
                     "MITO_second_allele": "|".join(second) or ".",
                     "MITO_evenements": ",".join(evs) or ".",
                     "MITO_priorite": best}
                fo.write("\t".join(f + [v[k] for k in extra[:-1]] + ["annote"]) + "\n")
                n["annote"] += 1
        for k, x in sorted(n.items()):
            print(f"           {k:18s} {x}", file=sys.stderr)
        print(f"           Résultat : {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
