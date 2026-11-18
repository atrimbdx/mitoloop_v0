#!/usr/bin/env python3
"""
gtex_isoformes.py - which promoter of a gene each tissue uses, from the GTEx
v8 transcript-level expression (RSEM, median TPM per tissue).

The transcripts of the gene are read from the GTF and grouped by TSS (TSSs
closer than --groupe bp form one promoter). The median TPM of each
transcript in each tissue comes from the GTEx portal API
(medianTranscriptExpression, dataset gtex_v8); transcripts are matched on
the Ensembl id without version. For each tissue the script reports the TPM
summed per promoter and its share of the gene.

RSEM splits the reads of isoforms that share most of their exons: the
share of each promoter is an estimate, to be read with that in mind.

Usage
  python3 scripts/gtex_isoformes.py --gene ABAT --gtf data/gencode.gtf.gz
  (--tissus "Brain_Frontal_Cortex_BA9,Kidney_Cortex,..." ; --json reponse.json
   reads a saved API answer instead of the network)
Output: printed table, and with --out a TSV (tissu, promoteur, tss,
transcrits, tpm, part).
"""

import argparse
import gzip
import json
import re
import sys
import urllib.parse
import urllib.request
from collections import defaultdict

API = "https://gtexportal.org/api/v2"
TISSUS = ("Brain_Frontal_Cortex_BA9,Brain_Cortex,Brain_Caudate_basal_ganglia,Brain_Cerebellum,"
          "Muscle_Skeletal,Heart_Left_Ventricle,Kidney_Cortex,Liver")


def get(path, **params):
    url = f"{API}{path}?{urllib.parse.urlencode(params, doseq=True)}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def transcrits(gtf, gene):
    """ENST sans version -> (nom, TSS, brin, gene_id)"""
    out = {}
    op = gzip.open if str(gtf).endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if len(f) < 9 or f[2] != "transcript" or f'gene_name "{gene}"' not in f[8]:
                continue
            tid = re.search(r'transcript_id "([^".]+)', f[8]).group(1)
            nom = re.search(r'transcript_name "([^"]+)"', f[8])
            gid = re.search(r'gene_id "([^".]+)', f[8]).group(1)
            tss = int(f[3]) if f[6] == "+" else int(f[4])
            out[tid] = (nom.group(1) if nom else tid, tss, f[6], gid)
    if not out:
        sys.exit(f"{gene} absent de {gtf}")
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gene", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--tissus", default=TISSUS, help="tissueSiteDetailId GTEx, séparés par des virgules")
    p.add_argument("--groupe", type=int, default=1000, help="TSS plus proches que N pb : un promoteur")
    p.add_argument("--json", default="", help="réponse de l'API enregistrée (au lieu du réseau)")
    p.add_argument("--out", default="")
    a = p.parse_args()
    tx = transcrits(a.gtf, a.gene)
    # promoteurs : TSS triés, regroupés
    tss = sorted({v[1] for v in tx.values()})
    prom, cur = {}, None
    for t in tss:
        if cur is None or t - cur[-1] > a.groupe:
            cur = [t]
        else:
            cur.append(t)
        prom[t] = cur
    nom_prom = {t: f"TSS_{min(g):,}" for t, g in prom.items()}
    tissus = [x for x in a.tissus.split(",") if x]
    if a.json:
        data = json.load(open(a.json))
    else:
        gid = next(iter(tx.values()))[3]
        ref = get("/reference/gene", geneId=gid, gencodeVersion="v26", genomeBuild="GRCh38/hg38")
        genes = ref.get("data", [])
        if not genes:
            sys.exit(f"{gid} introuvable dans GTEx (gencode v26)")
        gencode = genes[0]["gencodeId"]
        print(f"{a.gene} : {gencode} ; requête des médianes par transcrit...", file=sys.stderr)
        data = get("/expression/medianTranscriptExpression", gencodeId=gencode,
                   tissueSiteDetailId=tissus, datasetId="gtex_v8")
    rows = data.get("data", data if isinstance(data, list) else [])
    par = defaultdict(lambda: defaultdict(float))
    detail = defaultdict(lambda: defaultdict(list))
    absents = set()
    for r in rows:
        tid = str(r.get("transcriptId", "")).split(".")[0]
        t = r.get("tissueSiteDetailId", "")
        if tid not in tx:
            absents.add(tid)
            continue
        v = float(r.get("median", 0) or 0)
        pr = nom_prom[tx[tid][1]]
        par[t][pr] += v
        detail[t][pr].append((tx[tid][0], v))
    if absents:
        print(f"  {len(absents)} transcrit(s) GTEx absent(s) du GTF, ignoré(s) : "
              f"{', '.join(sorted(absents)[:6])}", file=sys.stderr)
    proms = sorted(set(nom_prom.values()), key=lambda s: int(s[4:].replace(",", "")))
    print(f"{'tissu':32s} " + " ".join(f"{p_:>16s}" for p_ in proms) + "   total TPM")
    lignes = []
    for t in tissus:
        tot = sum(par[t].values())
        cells = []
        for p_ in proms:
            v = par[t].get(p_, 0.0)
            cells.append(f"{v:7.1f} ({100 * v / tot:3.0f}%)" if tot else f"{'.':>16s}")
            lignes.append((t, p_, v, v / tot if tot else float("nan"),
                           ",".join(f"{n}:{x:.1f}" for n, x in detail[t][p_])))
        print(f"{t:32s} " + " ".join(f"{c:>16s}" for c in cells) + f"   {tot:7.1f}")
    print("\npromoteurs (TSS du GTF, regroupés à moins de "
          f"{a.groupe} pb) : " + " ; ".join(
              f"{p_} = {', '.join(sorted(v[0] for v in tx.values() if nom_prom[v[1]] == p_))}"
              for p_ in proms))
    if a.out:
        with open(a.out, "w") as fo:
            fo.write("tissu\tpromoteur\ttpm\tpart\ttranscrits\n")
            for t, p_, v, s, dtl in lignes:
                fo.write(f"{t}\t{p_}\t{v:.3f}\t{s:.3f}\t{dtl}\n")


if __name__ == "__main__":
    main()
