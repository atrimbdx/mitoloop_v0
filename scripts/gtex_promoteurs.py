#!/usr/bin/env python3
"""
gtex_promoteurs.py - for every gene of a list, the promoter each tissue uses,
from the GTEx v8 transcript-level medians (RSEM), compared with the promoter
of the MANE transcript that the AlphaGenome campaign uses.

The transcripts of each gene are read once from the GTF and grouped by TSS
(TSSs closer than --groupe bp form one promoter). For each gene, the GTEx
portal API gives the median TPM of every transcript in each tissue
(medianTranscriptExpression, gtex_v8); answers are cached in --cache, so an
interrupted run resumes without asking again. For each gene and tissue the
script reports the TPM summed per promoter, the share of the MANE promoter,
the dominant promoter and its share.

The share that a promoter deletion of the campaign actually removes is
part_couverte: the share of the gene's TPM carried by transcripts whose TSS
lies in --fenetre (AMONT:AVAL bp around the MANE TSS, gene strand; default
2500:100, the calibration deletion). Grouping by 1 kb can merge a TSS 1 kb
downstream with the MANE one; part_couverte does not.

Flag (column alerte): "promoteur_dominant_autre" when the dominant promoter of
the tissue is not the MANE one; "hors_deletion" when part_couverte is under
--couverte-min (default 0.5) in an expressed tissue; "peu_exprime" under --tpm-min total TPM;
"transcrits_absents" when more than --absents of the GTEx TPM comes from
transcripts absent from the GTF (GENCODE v26 against the GTF's version).
RSEM splits the reads of isoforms that share most of their exons: shares are
estimates.

Usage
  python3 scripts/gtex_promoteurs.py --genes genes.txt --gtf data/gencode.gtf.gz \
      --out results/gtex_promoteurs.tsv
Output columns: gene, type, tissu, tpm_total, part_mane, tss_mane,
promoteur_dominant (TSS of the most expressed transcript of the dominant
promoter, relative to the MANE TSS, gene strand), part_dominant,
transcrits_dominant, part_absents, alerte, part_couverte.
"""

import argparse
import gzip
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from pathlib import Path

API = "https://gtexportal.org/api/v2"
# types cellulaires de la campagne AlphaGenome -> tissu GTEx correspondant
TISSUS = ("muscle=Muscle_Skeletal,cardiomyocyte=Heart_Left_Ventricle,"
          "neurone=Brain_Frontal_Cortex_BA9,rein=Kidney_Cortex")


def get(path, essais=4, **params):
    url = f"{API}{path}?{urllib.parse.urlencode(params, doseq=True)}"
    for k in range(essais):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception as exc:                # réseau ou limite du portail : on attend
            if k == essais - 1:
                raise
            print(f"  {path} : {type(exc).__name__}, nouvel essai", file=sys.stderr)
            time.sleep(5 * (k + 1))


def lire_gtf(gtf, genes):
    """gène -> {ENST sans version: (nom, tss, brin, mane)}, et gène -> gene_id sans version"""
    tx, gid = defaultdict(dict), {}
    op = gzip.open if str(gtf).endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            m = re.search(r'gene_name "([^"]+)"', f[8])
            if not m or m.group(1) not in genes:
                continue
            g = m.group(1)
            t = re.search(r'transcript_id "([^".]+)', f[8]).group(1)
            nom = re.search(r'transcript_name "([^"]+)"', f[8])
            tss = int(f[3]) if f[6] == "+" else int(f[4])
            tx[g][t] = (nom.group(1) if nom else t, tss, f[6], "MANE_Select" in f[8])
            gid.setdefault(g, re.search(r'gene_id "([^".]+)', f[8]).group(1))
    return tx, gid


def grouper(tss_list, groupe):
    """TSS -> identifiant de promoteur (premier TSS du groupe)"""
    out, cur = {}, None
    for t in sorted(set(tss_list)):
        if cur is None or t - cur[-1] > groupe:
            cur = [t]
        else:
            cur.append(t)
        out[t] = cur
    return {t: min(g) for t, g in out.items()}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", required=True, help="gènes, séparés par des blancs")
    p.add_argument("--gtf", required=True)
    p.add_argument("--tissus", default=TISSUS, help="type=tissueSiteDetailId GTEx, séparés par des virgules")
    p.add_argument("--groupe", type=int, default=1000)
    p.add_argument("--tpm-min", type=float, default=1.0)
    p.add_argument("--absents", type=float, default=0.3, help="part des TPM de transcrits absents du GTF")
    p.add_argument("--fenetre", default="2500:100", help="AMONT:AVAL autour du TSS MANE (délétion)")
    p.add_argument("--couverte-min", type=float, default=0.5, dest="couverte_min")
    p.add_argument("--cache", default="results/gtex_cache")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    genes = set(Path(a.genes).read_text().split())
    amont, aval = (int(x) for x in a.fenetre.split(":"))
    tissus = dict(x.split("=", 1) for x in a.tissus.split(",") if "=" in x)
    tx, gid = lire_gtf(a.gtf, genes)
    cache = Path(a.cache)
    cache.mkdir(parents=True, exist_ok=True)
    manque = sorted(genes - set(tx))
    if manque:
        print(f"{len(manque)} gènes absents du GTF : {', '.join(manque[:8])}", file=sys.stderr)
    rows, n_alerte = [], defaultdict(int)
    for k, g in enumerate(sorted(tx), 1):
        fc = cache / f"{g}.json"
        if fc.exists():
            data = json.loads(fc.read_text())
        else:
            ref = get("/reference/gene", geneId=gid[g], gencodeVersion="v26", genomeBuild="GRCh38/hg38")
            if not ref.get("data"):
                fc.write_text(json.dumps({"data": [], "absent": True}))
                continue
            data = get("/expression/medianTranscriptExpression", gencodeId=ref["data"][0]["gencodeId"],
                       tissueSiteDetailId=list(tissus.values()), datasetId="gtex_v8")
            fc.write_text(json.dumps(data))
        if k % 100 == 0:
            print(f"  {k}/{len(tx)} gènes", file=sys.stderr)
        t_g = tx[g]
        mane = next((v for v in t_g.values() if v[3]), None)
        if mane is None:
            continue
        brin, tss_m = mane[2], mane[1]
        prom = grouper([v[1] for v in t_g.values()], a.groupe)
        p_mane = prom[tss_m]
        sens = 1 if brin == "+" else -1
        par = defaultdict(lambda: defaultdict(float))
        noms = defaultdict(lambda: defaultdict(list))
        top = defaultdict(dict)                       # tissu -> promoteur -> (tpm, tss) du transcrit le plus exprimé
        couv = defaultdict(float)
        abs_ = defaultdict(float)
        for r in data.get("data", []):
            t = str(r.get("transcriptId", "")).split(".")[0]
            ti = r.get("tissueSiteDetailId", "")
            v = float(r.get("median", 0) or 0)
            if t in t_g:
                pr = prom[t_g[t][1]]
                par[ti][pr] += v
                noms[ti][pr].append(t_g[t][0])
                if v > top[ti].get(pr, (-1, 0))[0]:
                    top[ti][pr] = (v, t_g[t][1])
                if -amont <= (t_g[t][1] - tss_m) * sens <= aval:
                    couv[ti] += v
            else:
                abs_[ti] += v
        for typ, ti in tissus.items():
            tot = sum(par[ti].values())
            tot_all = tot + abs_[ti]
            if tot:
                dom = max(par[ti], key=par[ti].get)
                part_m, part_d = par[ti].get(p_mane, 0.0) / tot, par[ti][dom] / tot
            else:
                dom, part_m, part_d = p_mane, float("nan"), float("nan")
            rel = (top[ti].get(dom, (0, dom))[1] - tss_m) * sens
            part_c = couv[ti] / tot if tot else float("nan")
            alerte = []
            if tot_all < a.tpm_min:
                alerte.append("peu_exprime")
            elif dom != p_mane:
                alerte.append("promoteur_dominant_autre")
            if tot_all >= a.tpm_min and tot and part_c < a.couverte_min:
                alerte.append("hors_deletion")
            if tot_all and abs_[ti] / tot_all > a.absents:
                alerte.append("transcrits_absents")
            for x in alerte:
                n_alerte[(typ, x)] += 1
            rows.append([g, typ, ti, f"{tot_all:.2f}", f"{part_m:.3f}", str(tss_m), f"{rel:+d}",
                         f"{part_d:.3f}", ",".join(sorted(set(noms[ti][dom]))) or ".",
                         f"{abs_[ti] / tot_all:.3f}" if tot_all else "nan", ",".join(alerte) or ".",
                         f"{part_c:.3f}"])
    with open(a.out, "w") as fo:
        fo.write("gene\ttype\ttissu\ttpm_total\tpart_mane\ttss_mane\tpromoteur_dominant\t"
                 "part_dominant\ttranscrits_dominant\tpart_absents\talerte\tpart_couverte\n")
        for r in rows:
            fo.write("\t".join(r) + "\n")
    print(f"{len({r[0] for r in rows})} gènes -> {a.out}", file=sys.stderr)
    for (typ, x), n in sorted(n_alerte.items()):
        print(f"  {typ:14s} {x:26s} {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
