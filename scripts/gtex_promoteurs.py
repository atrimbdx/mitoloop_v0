#!/usr/bin/env python3
"""
gtex_promoteurs.py - for every gene of a list, the promoter each tissue uses,
from the GTEx v8 transcript-level medians (RSEM), compared with the promoter
of the MANE transcript that the AlphaGenome campaign uses.

The transcripts of each gene are read from the GENCODE v26 GTF of GTEx v8
(--gtf-gtex; the TPMs belong to this annotation) and grouped by TSS (a group
stops when a TSS is more than --groupe bp from the FIRST TSS of the group, so a
group spans at most 1 kb). The project GTF (--gtf) gives only the MANE TSS and
the length of the MANE exon 1 (coordinates are GRCh38 in both). Without
--gtf-gtex, --gtf is used for everything (former behaviour). For each gene, the
GTEx portal API gives the median TPM of every transcript in each tissue
(medianTranscriptExpression, gtex_v8); answers are cached in --cache with the
tissues asked, and a gene whose cache lacks a tissue is asked again (the
network is used, never AlphaGenome). Shares are ratios of sums of per-transcript
medians, not per-sample shares; RSEM splits the reads of isoforms that share
most of their exons (and does not model the fragment GC bias, which can move
reads to 5'-truncated isoforms): shares are estimates.

Denominator: every GTEx transcript of the gene (tpm_total), including those
absent from the GTF (part_absents_v26, a version check expected near 0 with
--gtf-gtex).

Shares removed by the campaign's promoter deletions, with the real bounds
(downstream side capped at MANE exon 1 - 10 bp, as aval_max of ag2_balayage.py):
part_couverte_calib, TSS in [-AMONT, +min(AVAL, exon1 - 10)] around the MANE TSS
(--fenetre, default 2500:100, the calibration deletion); part_couverte_balayage,
the same with --fenetre-balayage (250:100, the scanning target). part_couverte
repeats part_couverte_calib (former name).

Tags of the GTF (--gtf-gtex): part_mrna_start_nf (share of the TPM carried by
transcripts tagged mRNA_start_NF), tags_dominant (tags of the most expressed
transcript of the dominant group: mRNA_start_NF, basic, tsl<N>).

Independent TSS positions (--cage FANTOM5 peaks, --epd EPDnew; never used to
judge an AlphaGenome score): cage_mane, cage_dominant, epd_mane, epd_dominant
("oui" when a peak or promoter lies within +/- --dist-pic bp of the TSS, same
strand when the BED has one), meme_pic_cage (MANE and dominant TSS in the same
CAGE peak). statut_tss, rule fixed beforehand, on the most expressed transcript
of the dominant group: "mane" when its TSS lies in the calibration deletion
(part_couverte_calib bounds); otherwise "dominant_confirme" when a CAGE peak lies
within +/- --dist-pic of its TSS and it is not tagged mRNA_start_NF, and
"dominant_non_confirme" when not (a TSS 233 bp downstream, in the MANE group of
1 kb, is thus judged on its own). EPD is reported but is not an
independent confirmation (it is derived from CAGE). Only a dominant_confirme
TSS may change a campaign target (tissus-A3).

Alerts (column alerte): "promoteur_dominant_autre" (dominant promoter not the
MANE one), "hors_deletion" (part_couverte_calib < --couverte-min in an expressed
tissue), "peu_exprime" (under --tpm-min total TPM), "transcrits_absents" (more
than --absents of the TPM from transcripts absent from the GTF),
"dominant_tronque_5p" (most expressed transcript of the dominant group, outside
the deletion, tagged mRNA_start_NF).

--tss-dominants OUT writes, for the campaign types (--types-campagne), the
genes x types with tpm_total >= --tpm-min, part_couverte_calib < 0.5 and
statut_tss = dominant_confirme: gene, type, chrom, tss_dom (0-based), brin,
fin_exon1_dom (0-based: exclusive end of exon 1 of the dominant transcript on
+, its start on -), rel (TSS relative to the MANE TSS, gene strand). Read by
ag2_balayage.py cibles and calibrer --tss-dominants.

Usage
  python3 scripts/gtex_promoteurs.py --genes genes.txt --gtf data/gencode.gtf.gz --gtf-gtex data/gencode.v26.annotation.gtf.gz --cage hg38_fair+new_CAGE_peaks_phase1and2.bed.gz --epd Hs_EPDnew_006_hg38.bed --out results/gtex_promoteurs.tsv --tss-dominants results/ag2/tss_dominants.tsv
Output columns: gene, type, tissu, tpm_total, part_mane, tss_mane,
promoteur_dominant (TSS of the most expressed transcript of the dominant
promoter, relative to the MANE TSS, gene strand), part_dominant,
transcrits_dominant, part_absents, alerte, part_couverte, part_couverte_calib,
part_couverte_balayage, part_mrna_start_nf, tags_dominant, cage_mane,
cage_dominant, epd_mane, epd_dominant, meme_pic_cage, statut_tss, chrom, brin,
tss_dominant_0, fin_exon1_dom.
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
# types cellulaires de la campagne AlphaGenome -> tissu GTEx correspondant ; fibroblaste
# et foie : tissus de confirmation clinique (RNA-seq de fibroblastes) et d'atteinte
# hépatique, GTEx seulement (aucune condition AlphaGenome)
TISSUS = ("muscle=Muscle_Skeletal,cardiomyocyte=Heart_Left_Ventricle,"
          "neurone=Brain_Frontal_Cortex_BA9,rein=Kidney_Cortex,"
          "fibroblaste=Cells_Cultured_fibroblasts,foie=Liver")
# tissus des caches écrits avant que le cache ne garde la liste demandée
TISSUS_V1 = {"Muscle_Skeletal", "Heart_Left_Ventricle", "Brain_Frontal_Cortex_BA9", "Kidney_Cortex"}
TYPES_CAMPAGNE = "muscle,cardiomyocyte,neurone,rein"


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


def attr(f9, cle):
    m = re.search(rf'{cle} "?([^";]+)"?;', f9)
    return m.group(1) if m else ""


def lire_gtf(gtf, genes, par_id=None):
    """gène -> {ENST sans version: dict(nom, tss 1-based, brin, mane, chrom, tags,
    fin_exon1 0-based, exon1)}, et gène -> gene_id sans version. tags : mRNA_start_NF,
    basic, tsl<N> ; fin_exon1 : fin exclue de l'exon 1 sur le brin +, début 0-based sur
    le brin - ; exon1 : longueur de l'exon 1. par_id (ENSG sans version -> nom) : gènes
    appariés par gene_id plutôt que par nom (GENCODE v26 : gènes renommés depuis,
    ex. ATP5A1 -> ATP5F1A)"""
    tx, gid = defaultdict(dict), {}
    ex1 = {}
    op = gzip.open if str(gtf).endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.split("\t")
            if len(f) < 9 or f[2] not in ("transcript", "exon"):
                continue
            if par_id is not None:
                mi = re.search(r'gene_id "([^".]+)', f[8])
                if not mi or mi.group(1) not in par_id:
                    continue
                g = par_id[mi.group(1)]
            else:
                m = re.search(r'gene_name "([^"]+)"', f[8])
                if not m or m.group(1) not in genes:
                    continue
                g = m.group(1)
            t = re.search(r'transcript_id "([^".]+)', f[8]).group(1)
            if f[2] == "exon":
                if re.search(r'exon_number "?1"?;', f[8]):
                    ex1[t] = (int(f[4]) if f[6] == "+" else int(f[3]) - 1, int(f[4]) - int(f[3]) + 1)
                continue
            nom = re.search(r'transcript_name "([^"]+)"', f[8])
            tags = set()
            if 'tag "mRNA_start_NF"' in f[8]:
                tags.add("mRNA_start_NF")
            if 'tag "basic"' in f[8]:
                tags.add("basic")
            tsl = attr(f[8], "transcript_support_level")
            if tsl and tsl != "NA":
                tags.add(f"tsl{tsl}")
            tss = int(f[3]) if f[6] == "+" else int(f[4])
            tx[g][t] = {"nom": nom.group(1) if nom else t, "tss": tss, "brin": f[6],
                        "mane": "MANE_Select" in f[8], "chrom": f[0], "tags": tags}
            gid.setdefault(g, re.search(r'gene_id "([^".]+)', f[8]).group(1))
    for g in tx:
        for t, v in tx[g].items():
            v["fin_exon1"], v["exon1"] = ex1.get(t, (None, None))
    return tx, gid


def grouper(tss_list, groupe):
    """TSS -> identifiant de promoteur (premier TSS du groupe) ; un groupe s'arrête quand
    un TSS est à plus de `groupe` pb du PREMIER TSS du groupe (étendue <= groupe)"""
    out, cur = {}, None
    for t in sorted(set(tss_list)):
        if cur is None or t - cur[0] > groupe:
            cur = [t]
        else:
            cur.append(t)
        out[t] = cur
    return {t: min(g) for t, g in out.items()}


def lire_pics(path):
    """BED (pics CAGE, promoteurs EPD) -> chrom -> [(début, fin, brin, nom)]"""
    out = defaultdict(list)
    if not path:
        return out
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for ln in fh:
            if ln.startswith(("#", "track", "browser")):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            out[f[0]].append((int(f[1]), int(f[2]), f[5] if len(f) > 5 and f[5] in "+-" else ".",
                              f[3] if len(f) > 3 else f"{f[0]}:{f[1]}"))
    return out


def pic_proche(pics, chrom, tss1, brin, dist):
    """nom du premier pic (même brin si le BED en a un) à +/- dist pb du TSS (1-based), ou None"""
    t0 = tss1 - 1
    for s, e, b, nm in pics.get(chrom, ()):
        if b in (".", brin) and s - dist <= t0 < e + dist:
            return nm
    return None


def couverte(tx_g, tss_m, sens, amont, aval):
    """ENST dont le TSS est dans [-amont, +aval] autour du TSS MANE (brin du gène)"""
    return {t for t, v in tx_g.items() if -amont <= (v["tss"] - tss_m) * sens <= aval}


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", required=True, help="gènes, séparés par des blancs")
    p.add_argument("--gtf", required=True, help="GTF du projet : TSS et exon 1 du transcrit MANE")
    p.add_argument("--gtf-gtex", default="", help="GTF GENCODE v26 de GTEx v8 (transcrits quantifiés)")
    p.add_argument("--tissus", default=TISSUS, help="type=tissueSiteDetailId GTEx, séparés par des virgules")
    p.add_argument("--groupe", type=int, default=1000)
    p.add_argument("--tpm-min", type=float, default=1.0)
    p.add_argument("--absents", type=float, default=0.3, help="part des TPM de transcrits absents du GTF")
    p.add_argument("--fenetre", default="2500:100", help="AMONT:AVAL autour du TSS MANE (délétion de calibration)")
    p.add_argument("--fenetre-balayage", default="250:100", help="AMONT:AVAL de la cible de balayage")
    p.add_argument("--couverte-min", type=float, default=0.5, dest="couverte_min")
    p.add_argument("--cage", default="", help="pics CAGE FANTOM5 (BED), positions indépendantes des TSS")
    p.add_argument("--epd", default="", help="promoteurs EPDnew (BED), rapportés seulement")
    p.add_argument("--dist-pic", type=int, default=100)
    p.add_argument("--cache", default="results/gtex_cache")
    p.add_argument("--out", required=True)
    p.add_argument("--tss-dominants", default="", help="sortie : TSS dominants confirmés hors délétion")
    p.add_argument("--types-campagne", default=TYPES_CAMPAGNE)
    a = p.parse_args()
    genes = set(Path(a.genes).read_text().split())
    amont, aval = (int(x) for x in a.fenetre.split(":"))
    amont_b, aval_b = (int(x) for x in a.fenetre_balayage.split(":"))
    tissus = dict(x.split("=", 1) for x in a.tissus.split(",") if "=" in x)
    tx_p, gid = lire_gtf(a.gtf, genes)
    tx = lire_gtf(a.gtf_gtex, genes, {v: k for k, v in gid.items()})[0] if a.gtf_gtex else tx_p
    cage, epd = lire_pics(a.cage), lire_pics(a.epd)
    cache = Path(a.cache)
    cache.mkdir(parents=True, exist_ok=True)
    manque = sorted(genes - set(tx_p))
    if manque:
        print(f"{len(manque)} gènes absents du GTF : {', '.join(manque[:8])}", file=sys.stderr)
    rows, n_alerte, dominants = [], defaultdict(int), []
    for k, g in enumerate(sorted(tx_p), 1):
        fc = cache / f"{g}.json"
        data = json.loads(fc.read_text()) if fc.exists() else None
        if data is not None and not data.get("absent"):
            vus = set(data.get("tissus", TISSUS_V1))
            if not set(tissus.values()) <= vus:
                data = None                          # tissu absent du cache : redemandé
        if data is None:
            ref = get("/reference/gene", geneId=gid[g], gencodeVersion="v26", genomeBuild="GRCh38/hg38")
            if not ref.get("data"):
                fc.write_text(json.dumps({"data": [], "absent": True}))
                continue
            data = get("/expression/medianTranscriptExpression", gencodeId=ref["data"][0]["gencodeId"],
                       tissueSiteDetailId=list(tissus.values()), datasetId="gtex_v8")
            data["tissus"] = sorted(tissus.values())
            fc.write_text(json.dumps(data))
        if k % 100 == 0:
            print(f"  {k}/{len(tx_p)} gènes", file=sys.stderr)
        mane = next((v for v in tx_p[g].values() if v["mane"]), None)
        if mane is None:
            continue
        t_g = tx.get(g, {})
        brin, tss_m = mane["brin"], mane["tss"]
        sens = 1 if brin == "+" else -1
        ex1_m = mane["exon1"]
        cap = (lambda av: av if ex1_m is None else max(0, min(av, ex1_m - 10)))
        c_calib = couverte(t_g, tss_m, sens, amont, cap(aval))
        c_bal = couverte(t_g, tss_m, sens, amont_b, cap(aval_b))
        prom = grouper([v["tss"] for v in t_g.values()] + [tss_m], a.groupe)
        p_mane = prom[tss_m]
        par = defaultdict(lambda: defaultdict(float))
        noms = defaultdict(lambda: defaultdict(list))
        top = defaultdict(dict)                       # tissu -> promoteur -> (tpm, ENST) du plus exprimé
        couv, couv_b, abs_, nf = (defaultdict(float) for _ in range(4))
        tot_all = defaultdict(float)
        for r in data.get("data", []):
            t = str(r.get("transcriptId", "")).split(".")[0]
            ti = r.get("tissueSiteDetailId", "")
            v = float(r.get("median", 0) or 0)
            tot_all[ti] += v
            if t in t_g:
                pr = prom[t_g[t]["tss"]]
                par[ti][pr] += v
                noms[ti][pr].append(t_g[t]["nom"])
                if v > top[ti].get(pr, (-1, None))[0]:
                    top[ti][pr] = (v, t)
                couv[ti] += v * (t in c_calib)
                couv_b[ti] += v * (t in c_bal)
                nf[ti] += v * ("mRNA_start_NF" in t_g[t]["tags"])
            else:
                abs_[ti] += v
        pm_mane = pic_proche(cage, mane["chrom"], tss_m, brin, a.dist_pic)
        e_mane = pic_proche(epd, mane["chrom"], tss_m, brin, a.dist_pic)
        for typ, ti in tissus.items():
            tot = tot_all[ti]
            if par[ti] and tot > 0:            # tissu où tous les transcrits sont à 0 : non évaluable
                dom = max(par[ti], key=par[ti].get)
                part_m, part_d = par[ti].get(p_mane, 0.0) / tot, par[ti][dom] / tot
            else:
                dom, part_m, part_d = p_mane, float("nan"), float("nan")
            t_dom = top[ti].get(dom, (0, None))[1]
            vd = t_g.get(t_dom) if t_dom else None
            tss_d = vd["tss"] if vd else tss_m
            rel = (tss_d - tss_m) * sens
            part_c = couv[ti] / tot if tot else float("nan")
            part_cb = couv_b[ti] / tot if tot else float("nan")
            tags_d = ",".join(sorted(vd["tags"])) if vd and vd["tags"] else "."
            pc_d = pic_proche(cage, vd["chrom"], tss_d, brin, a.dist_pic) if vd else None
            ep_d = pic_proche(epd, vd["chrom"], tss_d, brin, a.dist_pic) if vd else None
            # TSS du transcrit le plus exprimé du groupe dominant couvert par la délétion de
            # calibration : « mane » (même s'il n'est pas le TSS MANE lui-même)
            hors_d = vd is not None and t_dom not in c_calib
            if not hors_d:
                statut = "mane"
            elif pc_d and vd and "mRNA_start_NF" not in vd["tags"]:
                statut = "dominant_confirme"
            else:
                statut = "dominant_non_confirme"
            alerte = []
            if tot < a.tpm_min:
                alerte.append("peu_exprime")
            elif dom != p_mane:
                alerte.append("promoteur_dominant_autre")
            if tot >= a.tpm_min and tot and part_c < a.couverte_min:
                alerte.append("hors_deletion")
            if tot and abs_[ti] / tot > a.absents:
                alerte.append("transcrits_absents")
            if hors_d and "mRNA_start_NF" in vd["tags"]:
                alerte.append("dominant_tronque_5p")
            for x in alerte:
                n_alerte[(typ, x)] += 1
            oui = (lambda x: "oui" if x else ("non" if a.cage else "."))
            fin1 = vd["fin_exon1"] if vd and vd["fin_exon1"] is not None else None
            rows.append([g, typ, ti, f"{tot:.2f}", f"{part_m:.3f}", str(tss_m), f"{rel:+d}",
                         f"{part_d:.3f}", ",".join(sorted(set(noms[ti][dom]))) or ".",
                         f"{abs_[ti] / tot:.3f}" if tot else "nan", ",".join(alerte) or ".",
                         f"{part_c:.3f}", f"{part_c:.3f}", f"{part_cb:.3f}",
                         f"{nf[ti] / tot:.3f}" if tot else "nan", tags_d,
                         oui(pm_mane), oui(pc_d), ("oui" if e_mane else "non") if a.epd else ".",
                         ("oui" if ep_d else "non") if a.epd else ".",
                         ("oui" if (pm_mane and pc_d == pm_mane) else "non") if a.cage else ".",
                         statut, mane["chrom"], brin, str(tss_d - 1), "." if fin1 is None else str(fin1)])
            if (typ in a.types_campagne.split(",") and tot >= a.tpm_min and part_c == part_c
                    and part_c < a.couverte_min and statut == "dominant_confirme" and fin1 is not None):
                dominants.append([g, typ, mane["chrom"], str(tss_d - 1), brin, str(fin1), str(rel)])
    cols = ["gene", "type", "tissu", "tpm_total", "part_mane", "tss_mane", "promoteur_dominant",
            "part_dominant", "transcrits_dominant", "part_absents", "alerte", "part_couverte",
            "part_couverte_calib", "part_couverte_balayage", "part_mrna_start_nf", "tags_dominant",
            "cage_mane", "cage_dominant", "epd_mane", "epd_dominant", "meme_pic_cage", "statut_tss",
            "chrom", "brin", "tss_dominant_0", "fin_exon1_dom"]
    with open(a.out, "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            fo.write("\t".join(r) + "\n")
    print(f"{len({r[0] for r in rows})} gènes -> {a.out}", file=sys.stderr)
    for (typ, x), n in sorted(n_alerte.items()):
        print(f"  {typ:14s} {x:26s} {n}", file=sys.stderr)
    if a.tss_dominants:
        with open(a.tss_dominants, "w") as fo:
            fo.write("gene\ttype\tchrom\ttss_dom\tbrin\tfin_exon1_dom\trel\n")
            for r in dominants:
                fo.write("\t".join(r) + "\n")
        par_t = defaultdict(int)
        for r in dominants:
            par_t[r[1]] += 1
        print(f"TSS dominants confirmés hors délétion : {len(dominants)} ("
              + ", ".join(f"{t} {n}" for t, n in sorted(par_t.items())) + f") -> {a.tss_dominants}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
