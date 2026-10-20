#!/usr/bin/env python3
"""
positifs_mitoloop.py — MCIA (accès internet). Prépare un TSV de variants non
codants pathogènes connus dans les gènes MitoLoop, à passer ensuite sur la VM
comme n'importe quelle entrée :
    bash reannot_complet.sh positifs_mitoloop.tsv positifs

Sources
  litterature  variants publiés avec effet démontré sur l'expression (liste
               ci-dessous, PMID à l'appui), plus ceux de --litterature-extra
  clinvar      ClinVar GRCh38 : Pathogenic / Likely_pathogenic (hors
               Conflicting), gène dans --genes, conséquences seulement non
               codantes (5'UTR, 3'UTR, amont/aval ; introniques avec
               --avec-introniques)
  tables       --tables : tableaux supplémentaires (xlsx, tsv, csv), par
               exemple ceux de Kremer 2017, Yépez 2022 ou Kopajtich (équipe
               Prokisch), à télécharger à la main. Chaque ligne citant une
               notation c. (avec NM_/ENST, ou un symbole de gène dans une autre
               cellule) est lue ; seules les notations non codantes sont gardées
               (c.-N, c.*N ; introniques avec --avec-introniques).
               Les données Zenodo de Prokisch (KNOWN_MUTATION) ne donnent que
               le gène, pas le variant : elles ne servent pas ici.

Les notations c. sont converties en position hg38 par l'API Ensembl
(variant_recoder, qui vérifie aussi l'allèle de référence), sur le transcrit
MANE Select du gène quand le transcrit n'est pas donné (colonne
transcrit_suppose = oui). Les réponses sont gardées dans --cache.

Sortie (TSV avec en-tête ; reconnu par preparer_entree.py par les colonnes
chrom/pos/ref/alt et transcrit ENST MANE) :
  chrom pos ref alt transcrit gene categorie source etoiles clinvar_id
  transcrit_suppose notation_publiee preuve
notation_publiee est écrite « NM_x (c.y) » pour ne pas être prise pour la
colonne variant. Les variants non convertis vont dans PREFIXE.non_convertis.tsv.

Usage
  python3 scripts/positifs_mitoloop.py --genes genes.txt [--out positifs_mitoloop.tsv]
      [--renamed data/genes_renamed.tsv] [--tables S1.xlsx S2.tsv ...]
      [--litterature-extra extra.tsv] [--avec-introniques] [--etoiles-min 0]
      [--clinvar clinvar.vcf.gz] [--mane MANE.summary.txt.gz] [--datadir data/positifs]
--litterature-extra : TSV gene, notation (NM_...:c... ou c...), source, preuve.
Bibliothèque standard uniquement.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "reannot"))

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
MANE_DIR = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/"
RECODER = "https://rest.ensembl.org/variant_recoder/homo_sapiens?fields=vcf_string"

# Variants publiés, effet sur l'expression démontré (PubMed)
LITTERATURE = [
    ("ETFA", "c.-85G>A", "PMID:36579410 Am J Med Genet A (doi:10.1002/ajmg.a.63104)",
     "MADD ; homozygote ; ETFA diminuée (lymphoblastes), promoteur diminué (luciférase)"),
    ("RARS2", "c.-2A>G", "PMID:25809939 J Hum Genet 2015 (doi:10.1038/jhg.2015.31)",
     "PCH6 ; homozygote ; ARNm diminué, activité promotrice diminuée"),
    ("BCS1L", "c.-122G>T", "PMID:40660675 Ann Clin Transl Neurol 2025 (doi:10.1002/acn3.70108)",
     "en trans de p.L280F ; expression résiduelle de l'allèle normal (fibroblastes)"),
]

# Conséquences ClinVar (MC) acceptées
MC_NONCODANT = {"5_prime_UTR_variant": "5UTR_promoteur", "3_prime_UTR_variant": "3UTR",
                "2KB_upstream_variant": "5UTR_promoteur", "upstream_gene_variant": "5UTR_promoteur",
                "500B_downstream_variant": "3UTR", "downstream_gene_variant": "3UTR",
                "non-coding_transcript_variant": "non_codant"}
MC_INTRON = {"intron_variant": "intronique"}
ETOILES = {"practice_guideline": 4, "reviewed_by_expert_panel": 3,
           "criteria_provided,_multiple_submitters,_no_conflicts": 2,
           "criteria_provided,_single_submitter": 1}

C_RE = re.compile(r"c\.(?P<pos>[-*]?\d+(?:[-+]\d+)?(?:_[-*]?\d+(?:[-+]\d+)?)?)"
                  r"(?P<chg>(?:[ACGT]+>[ACGT]+|del[ACGT]*ins[ACGT]+|del[ACGT]*|dup[ACGT]*|ins[ACGT]+))")
TX_RE = re.compile(r"\b((?:NM_|ENST)\d+(?:\.\d+)?)")


def log(msg):
    print(msg, file=sys.stderr)


def download(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    log(f"téléchargement {url}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as fh:
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            fh.write(b)
    tmp.rename(dest)
    return dest


def mane_file(datadir: Path) -> Path:
    with urllib.request.urlopen(MANE_DIR, timeout=60) as r:
        html = r.read().decode("utf-8", "replace")
    names = sorted(set(re.findall(r"MANE\.GRCh38\.v[\d.]+\.summary\.txt\.gz", html)))
    if not names:
        sys.exit(f"fichier MANE summary introuvable dans {MANE_DIR} (donner --mane)")
    return download(MANE_DIR + names[-1], datadir / names[-1])


def read_mane(path: Path):
    """symbole -> (NM versionné, ENST sans version) des MANE Select"""
    out = {}
    with gzip.open(path, "rt") as fh:
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ix["MANE_status"]] != "MANE Select":
                continue
            out[f[ix["symbol"]]] = (f[ix["RefSeq_nuc"]], f[ix["Ensembl_nuc"]].split(".")[0])
    return out


def read_genes(path: Path, renamed: Path | None):
    genes = {l.strip().split("\t")[0] for l in open(path) if l.strip() and not l.startswith("#")}
    if renamed and renamed.exists():
        for r in csv.DictReader(open(renamed), delimiter="\t"):
            if r["gtf_symbol"] in genes or r["symbol"] in genes:
                genes |= {r["symbol"], r["gtf_symbol"]}
    return genes


def categorie_c(c: str, introns: bool):
    """catégorie d'une notation c. non codante, None si codante ou rejetée"""
    m = C_RE.search(c)
    if not m:
        return None
    first = m.group("pos").split("_")[0]
    offset = re.search(r"\d([-+]\d+)$", first)
    if offset:
        return "intronique" if introns and abs(int(offset.group(1))) > 20 else None
    if first.startswith("-"):
        return "5UTR_promoteur"
    if first.startswith("*"):
        return "3UTR"
    return None


class Recoder:
    """notation HGVS -> [(chrom, pos, ref, alt)] par Ensembl, avec cache"""

    def __init__(self, cache: Path):
        self.cache = cache
        self.mem = json.loads(cache.read_text()) if cache.exists() else {}

    def __call__(self, hgvs_list):
        todo = [h for h in dict.fromkeys(hgvs_list) if h not in self.mem]
        for i in range(0, len(todo), 100):
            batch = todo[i:i + 100]
            try:
                res = self._post(batch)
            except Exception as e:                      # réseau, API
                log(f"  Ensembl variant_recoder indisponible ({e}) ; {len(batch)} notations non converties")
                continue
            for h in batch:
                self.mem[h] = res.get(h, [])
            self.cache.write_text(json.dumps(self.mem, indent=0))
            time.sleep(0.5)
        return {h: [tuple(x) for x in self.mem.get(h, [])] for h in hgvs_list}

    @staticmethod
    def _post(batch):
        req = urllib.request.Request(RECODER, data=json.dumps({"ids": batch}).encode(),
                                     headers={"Content-Type": "application/json",
                                              "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=120) as r:
            return parse_recoder(json.loads(r.read()))


def parse_recoder(data):
    """réponse de variant_recoder -> {entrée: [[chrom, pos, ref, alt], ...]}"""
    out = {}

    def walk(o):
        if isinstance(o, dict):
            if "input" in o and "vcf_string" in o:
                for s in o["vcf_string"]:
                    p = s.split("-")
                    if len(p) == 4 and p[1].isdigit():
                        out.setdefault(o["input"], [])
                        if p not in out[o["input"]]:
                            out[o["input"]].append(p)
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
    walk(data)
    return out


def chrom_ok(c):
    c = c[3:] if c.startswith("chr") else c
    return c if re.fullmatch(r"[0-9]{1,2}|X|Y|MT|M", c) else None


def clinvar_rows(path, genes, introns, min_stars):
    keep = dict(MC_NONCODANT, **(MC_INTRON if introns else {}))
    n = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##fileDate"):
                log(f"ClinVar {line.strip().split('=')[1]}")
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
            sig = info.get("CLNSIG", "")
            if not re.match(r"(Pathogenic|Likely_pathogenic)", sig) or "Conflicting" in sig:
                continue
            gi = [g.split(":")[0] for g in info.get("GENEINFO", "").split("|") if g]
            g = [x for x in gi if x in genes]
            if not g or not re.fullmatch(r"[ACGT]+", f[3]) or not re.fullmatch(r"[ACGT]+", f[4]):
                continue
            mc = [m.split("|")[1] for m in info.get("MC", "").split(",") if "|" in m]
            if not mc or any(m not in keep for m in mc):
                continue
            stars = ETOILES.get(info.get("CLNREVSTAT", ""), 0)
            if stars < min_stars:
                continue
            c = chrom_ok(f[0])
            if not c:
                continue
            n += 1
            cat = sorted({keep[m] for m in mc})[0]
            yield dict(chrom=c, pos=f[1], ref=f[3], alt=f[4], gene=g[0], categorie=cat,
                       source="clinvar", etoiles=str(stars), clinvar_id=f[2],
                       preuve=f"{sig} ; {info.get('CLNDN', '.').replace('_', ' ')}"[:300])
    log(f"ClinVar : {n} variants retenus")


def table_rows(path: Path, mane, introns):
    from preparer_entree import read_table
    sheets = [None]
    if path.suffix.lower() == ".xlsx":
        import zipfile
        import xml.etree.ElementTree as ET
        from preparer_entree import NS
        wb = ET.fromstring(zipfile.ZipFile(path).read("xl/workbook.xml"))
        sheets = [s.get("name") for s in wb.find("m:sheets", NS)]
    for sh in sheets:
        for i, row in enumerate(read_table(path, sh)):
            cells = [str(x).strip() for x in row]
            sym = next((x for x in cells if x in mane), None) or next(
                (g for x in cells for g in re.findall(r"\(([A-Za-z0-9-]+)\)", x) if g in mane), None)
            for cell in cells:
                for m in C_RE.finditer(cell):
                    cat = categorie_c(m.group(0), introns)
                    if not cat:
                        continue
                    tx = TX_RE.search(cell) or next(
                        (TX_RE.search(x) for x in cells if TX_RE.search(x)), None)
                    yield dict(gene=sym or ".", c=m.group(0), tx=tx.group(1) if tx else None,
                               categorie=cat, source=f"{path.name}" + (f"[{sh}]" if sh else ""),
                               preuve=f"ligne {i + 1}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--genes", type=Path, required=True)
    p.add_argument("--renamed", type=Path, default=Path("data/genes_renamed.tsv"))
    p.add_argument("--out", type=Path, default=Path("positifs_mitoloop.tsv"))
    p.add_argument("--datadir", type=Path, default=Path("data/positifs"))
    p.add_argument("--clinvar", type=Path)
    p.add_argument("--mane", type=Path)
    p.add_argument("--cache", type=Path)
    p.add_argument("--tables", type=Path, nargs="*", default=[])
    p.add_argument("--litterature-extra", type=Path)
    p.add_argument("--avec-introniques", action="store_true")
    p.add_argument("--etoiles-min", type=int, default=0)
    p.add_argument("--sans-clinvar", action="store_true")
    a = p.parse_args()
    a.datadir.mkdir(parents=True, exist_ok=True)

    genes = read_genes(a.genes, a.renamed)
    mane = read_mane(a.mane or mane_file(a.datadir))
    log(f"{len(genes)} gènes ; {sum(g in mane for g in genes)} avec un MANE Select")
    recode = Recoder(a.cache or a.datadir / "variant_recoder.json")

    # notations c. -> à convertir
    todo = [dict(gene=g, c=c, tx=None, categorie=categorie_c(c, True), source="litterature",
                 preuve=f"{ref} ; {note}") for g, c, ref, note in LITTERATURE]
    if a.litterature_extra:
        for r in csv.reader(open(a.litterature_extra), delimiter="\t"):
            if len(r) < 2 or r[0].startswith("#") or r[0].lower() == "gene":
                continue
            tx = TX_RE.search(r[1])
            todo.append(dict(gene=r[0], c=C_RE.search(r[1]).group(0), tx=tx.group(1) if tx else None,
                             categorie=categorie_c(r[1], True) or "autre",
                             source=r[2] if len(r) > 2 else "litterature",
                             preuve=r[3] if len(r) > 3 else "."))
    for t in a.tables:
        rows = [r for r in table_rows(t, mane, a.avec_introniques)
                if r["gene"] in genes or (r["gene"] == "." and r["tx"])]
        log(f"{t.name} : {len(rows)} notations non codantes")
        todo += rows

    rows, bad = [], []
    for r in todo:
        r["suppose"] = "non" if r["tx"] else "oui"
        if not r["tx"]:
            if r["gene"] not in mane:
                bad.append((r, "gène sans MANE Select"))
                continue
            r["tx"] = mane[r["gene"]][0]
        r["hgvs"] = f"{r['tx']}:{r['c']}"
    conv = recode([r["hgvs"] for r in todo if "hgvs" in r])
    # version du NM inconnue d'Ensembl : nouvel essai sans version
    retry = {r["hgvs"]: re.sub(r"^(NM_\d+)\.\d+", r"\1", r["hgvs"]) for r in todo
             if "hgvs" in r and not conv.get(r["hgvs"]) and re.match(r"NM_\d+\.", r["hgvs"])}
    if retry:
        conv2 = recode(list(retry.values()))
        for h, h2 in retry.items():
            conv[h] = conv2.get(h2) or []
    for r in todo:
        if "hgvs" not in r:
            continue
        hits = conv.get(r["hgvs"]) or []
        if not hits:
            bad.append((r, "non converti (Ensembl)"))
            continue
        for c, pos, ref, alt in hits:
            rows.append(dict(chrom=c, pos=pos, ref=ref, alt=alt, gene=r["gene"],
                             categorie=r["categorie"], source=r["source"], etoiles=".",
                             clinvar_id=".", transcrit_suppose=r["suppose"],
                             notation_publiee=f"{r['tx']} ({r['c']})", preuve=r["preuve"]))

    if not a.sans_clinvar:
        cv = a.clinvar or download(CLINVAR_URL, a.datadir / "clinvar.vcf.gz")
        rows += [dict(r, transcrit_suppose=".", notation_publiee=".")
                 for r in clinvar_rows(cv, genes, a.avec_introniques, a.etoiles_min)]

    # fusion par position (sources réunies)
    merged = {}
    for r in rows:
        r["chrom"] = chrom_ok(r["chrom"]) or r["chrom"]
        k = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
        if k in merged:
            m = merged[k]
            for f in ("source", "etoiles", "clinvar_id", "notation_publiee", "preuve"):
                vals = [v for v in (m[f].split(" | ") + [r[f]]) if v != "."]
                m[f] = " | ".join(dict.fromkeys(vals)) or "."
        else:
            merged[k] = dict(r)
    order = lambda k: (int(k[0]) if k[0].isdigit() else 100 + ord(k[0][0]), k[1])
    cols = ["chrom", "pos", "ref", "alt", "transcrit", "gene", "categorie", "source",
            "etoiles", "clinvar_id", "transcrit_suppose", "notation_publiee", "preuve"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for k in sorted(merged, key=order):
            r = merged[k]
            r["chrom"] = "chr" + r["chrom"]
            r["transcrit"] = mane.get(r["gene"], (None, "."))[1]
            fh.write("\t".join(str(r.get(c, ".")).replace("\t", " ") for c in cols) + "\n")
    nc = a.out.with_name(a.out.name.replace(".tsv", "") + ".non_convertis.tsv")
    with open(nc, "w") as fh:
        fh.write("gene\tnotation\tsource\traison\n")
        for r, why in bad:
            fh.write(f"{r['gene']}\t{r.get('hgvs', r['c'])}\t{r['source']}\t{why}\n")
    by = {}
    for r in merged.values():
        for s in r["source"].split(" | "):
            by[s.split("[")[0]] = by.get(s.split("[")[0], 0) + 1
    log(f"{len(merged)} variants -> {a.out} ({', '.join(f'{s} {n}' for s, n in sorted(by.items()))})")
    if bad:
        log(f"{len(bad)} notations non converties -> {nc}")


if __name__ == "__main__":
    main()
