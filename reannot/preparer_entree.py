#!/usr/bin/env python3
"""
preparer_entree.py — VM. Écrit l'entrée de reannot_alphagenome.sh (TSV sans
en-tête : transcrit Ensembl, variant en notation g. / c. / n.) à partir d'un
tableur Excel (.xlsx), d'un TSV ou d'un CSV avec en-tête.

Les colonnes sont reconnues par leur contenu, sur les 200 premières lignes :
  transcrit   valeurs ENST... (version facultative)
  variant     notations c. / n. / g. (chr2:g.10543284T>A, NM_...:c.541-6T>G
              est accepté : la partie avant ':' est retirée sauf pour g.)
ou, à défaut de colonne variant, les colonnes chrom, pos, ref, alt (noms
usuels : chrom/chr/chromosome, pos/position/start, ref, alt), converties en
notation g. (substitution, del, ins, delins). Un VCF (.vcf, .vcf.gz) est lu
de la même façon (une ligne par allèle alternatif).

Sans colonne transcrit, avec --gtf (le GTF Ensembl du pipeline), le transcrit
est le MANE Select qui chevauche le variant, sinon celui dont le TSS est le
plus proche (à moins de --max-dist pb) : pour un variant régulateur, le gène
voisin. Sinon il est laissé à "." (VEP annote le variant, mais la sélection
finale du transcrit ne se fait pas : statut transcrit_absent_du_cache).

Notation c. / n. sans transcrit, avec une colonne chromosome (1..22, X, Y, MT,
avec ou sans chr) et une colonne position (GRCh38), reconnues par leur contenu
ou leur nom : le transcrit est le MANE Select à cette position (--gtf). Comme
la notation c. a pu être écrite sur un autre transcrit, les positions saisies
sont gardées (SORTIE.positions.tsv) et comparées, après la conversion, à la
position obtenue : `--controle-positions SORTIE.positions.tsv FINAL.tsv` ajoute
à FINAL.tsv la colonne controle_position (concordante, discordante
(saisie chrN:pos), ou -). reannot_complet.sh le fait de lui-même.

Usage : python3 preparer_entree.py entree.xlsx|.tsv|.csv|.vcf.gz sortie.tsv
        [--gtf Homo_sapiens.GRCh38.NNN.gtf.gz] [--feuille NOM]
        python3 preparer_entree.py --controle-positions sortie.positions.tsv run.final.tsv
Bibliothèque standard uniquement (le .xlsx est lu comme une archive XML).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}
TX = re.compile(r"^ENST\d{11}(\.\d+)?$")
HGVS = re.compile(r"^(?:[^:\s]+:)?(?:(?:chr)?[0-9XYM]{1,2}T?:)?[cng]\.\S+$")
G_HGVS = re.compile(r"^(chr)?[0-9XYM]+T?:g\.")
CHROM = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$", re.I)
POS = re.compile(r"^\d{2,9}(\.0+)?$")


def chrom_name(c):
    c = c.strip()
    c = c[3:] if c.lower().startswith("chr") else c
    return "chrM" if c.upper() in ("M", "MT") else "chr" + c.upper()


def controle_positions(pos_file: Path, final: Path):
    """ajoute à FINAL.tsv la colonne controle_position : position obtenue par la
    conversion contre position saisie (substitution : la même base ; autre
    variant : à 2 pb près, l'ancrage VCF d'un indel)"""
    given = {}
    with open(pos_file) as fh:
        for line in fh:
            t, v, c, p = line.rstrip("\n").split("\t")
            given.setdefault((t, v), (chrom_name(c), int(p)))
    lines = final.read_text().splitlines()
    head = lines[0].split("\t")
    col = {h: i for i, h in enumerate(head)}
    n = {"concordante": 0, "discordante": 0, "-": 0}
    out = ["\t".join(head + ["controle_position"])]
    for line in lines[1:]:
        f = line.split("\t")
        g = given.get((f[col["transcrit"]], f[col["variant_entree"]]))
        res = "-"
        if g and f[col["pos"]].isdigit():
            d = abs(int(f[col["pos"]]) - g[1])
            ok = chrom_name(f[col["chrom"]]) == g[0] and (
                d == 0 or (">" not in f[col["variant_entree"]] and d <= 2))
            res = "concordante" if ok else f"discordante (saisie {g[0]}:{g[1]})"
        n[res.split(" ")[0]] += 1
        out.append("\t".join(f + [res]))
    final.write_text("\n".join(out) + "\n")
    print(f"contrôle des positions saisies : {n['concordante']} concordantes, "
          f"{n['discordante']} discordantes, {n['-']} sans position")
    if n["discordante"]:
        print("  discordantes : la notation c./n. ne porte sans doute pas sur le MANE Select "
              "(donner le transcrit dans une colonne), ou les positions ne sont pas en GRCh38")


def col_index(ref: str) -> int:
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0):
        n = n * 26 + ord(ch) - 64
    return n - 1


def read_xlsx(path: Path, sheet: str | None):
    z = zipfile.ZipFile(path)
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
            shared.append("".join(t.text or "" for t in si.iter(f"{{{NS['m']}}}t")))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels}
    sheets = [(s.get("name"), target[s.get(f"{{{NS['r']}}}id")])
              for s in wb.find("m:sheets", NS)]
    if sheet:
        hit = [t for n, t in sheets if n == sheet]
        if not hit:
            sys.exit(f"feuille {sheet!r} absente ; feuilles : {[n for n, _ in sheets]}")
        tgt = hit[0]
    else:
        tgt = sheets[0][1]
    tgt = tgt.lstrip("/")
    tgt = tgt if tgt.startswith("xl/") else "xl/" + tgt
    rows = []
    for row in ET.fromstring(z.read(tgt)).iter(f"{{{NS['m']}}}row"):
        cells = {}
        for c in row.findall("m:c", NS):
            t = c.get("t")
            v = c.find("m:v", NS)
            if t == "s" and v is not None:
                val = shared[int(v.text)]
            elif t == "inlineStr":
                val = "".join(x.text or "" for x in c.iter(f"{{{NS['m']}}}t"))
            else:
                val = v.text if v is not None else ""
            cells[col_index(c.get("r"))] = (val or "").strip()
        if cells:
            n = max(cells) + 1
            rows.append([cells.get(i, "") for i in range(n)])
    return rows


def read_vcf(path: Path):
    import gzip
    rows = [["chrom", "pos", "ref", "alt"]]
    op = gzip.open if path.name.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            for alt in f[4].split(","):
                if re.fullmatch(r"[ACGTNacgtn]+", alt):      # pas les SV symboliques
                    rows.append([f[0], f[1], f[3], alt])
    return rows


def mane_index(gtf: Path):
    """chrom -> [(start, end, tss, transcrit)] des transcrits MANE Select"""
    import gzip
    idx = {}
    op = gzip.open if gtf.name.endswith(".gz") else open
    with op(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#" or "MANE_Select" not in line:
                continue
            f = line.split("\t")
            if f[2] != "transcript":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m:
                continue
            c = f[0][3:] if f[0].startswith("chr") else f[0]
            s, e = int(f[3]), int(f[4])
            idx.setdefault(c, []).append((s, e, s if f[6] == "+" else e, m.group(1)))
    return idx


def pick_mane(idx, chrom, pos, max_dist):
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    best = None
    for s, e, tss, t in idx.get(c, []):
        d = 0 if s <= pos <= e else abs(pos - tss)
        if d <= max_dist and (best is None or d < best[0]):
            best = (d, t)
    return best[1] if best else "."


def read_table(path: Path, sheet):
    if path.suffix.lower() == ".xlsx":
        return read_xlsx(path, sheet)
    if path.name.lower().endswith((".vcf", ".vcf.gz")):
        return read_vcf(path)
    with open(path, newline="", encoding="utf-8-sig") as fh:
        text = fh.read()
    sep = "\t" if text.count("\t") >= text.count(",") else ","
    return [r for r in csv.reader(text.splitlines(), delimiter=sep)]


def g_notation(c, p, r, a):
    c = c[3:] if c.lower().startswith("chr") else c
    p = int(float(p))
    r, a = r.upper(), a.upper()
    if len(r) == 1 and len(a) == 1:
        return f"chr{c}:g.{p}{r}>{a}"
    # préfixe commun (représentation VCF)
    k = 0
    while k < min(len(r), len(a)) and r[k] == a[k]:
        k += 1
    r2, a2, s = r[k:], a[k:], p + k
    if not a2:
        return f"chr{c}:g.{s}del" if len(r2) == 1 else f"chr{c}:g.{s}_{s + len(r2) - 1}del"
    if not r2:
        return f"chr{c}:g.{s - 1}_{s}ins{a2}"
    e = s + len(r2) - 1
    return f"chr{c}:g.{s}delins{a2}" if s == e else f"chr{c}:g.{s}_{e}delins{a2}"


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("entree", type=Path)
    p.add_argument("sortie", type=Path)
    p.add_argument("--controle-positions", action="store_true",
                   help="entree = SORTIE.positions.tsv, sortie = FINAL.tsv à compléter")
    p.add_argument("--feuille")
    p.add_argument("--gtf", type=Path)
    p.add_argument("--max-dist", type=int, default=1_000_000)
    a = p.parse_args()
    if a.controle_positions:
        controle_positions(a.entree, a.sortie)
        return
    rows = [r for r in read_table(a.entree, a.feuille) if any(x.strip() for x in r)]
    if not rows:
        sys.exit(f"{a.entree} : vide")
    width = max(len(r) for r in rows)
    rows = [r + [""] * (width - len(r)) for r in rows]
    probe = rows[:200]

    def share(j, rx):
        vals = [r[j].strip() for r in probe if r[j].strip()]
        return sum(bool(rx.match(v.replace(" ", ""))) for v in vals) / max(len(vals), 1)

    tx = max(range(width), key=lambda j: share(j, TX))
    tx = tx if share(tx, TX) > 0.5 else None
    var = max(range(width), key=lambda j: share(j, HGVS))
    var = var if share(var, HGVS) > 0.5 else None
    head = [x.strip().lower() for x in rows[0]]
    names = {"chrom": ("chrom", "chr", "chromosome", "#chrom", "contig"),
             "pos": ("pos", "position", "start", "pos_hg38", "position_hg38"),
             "ref": ("ref", "reference"), "alt": ("alt", "alternate", "alternative")}
    vcfcols = {k: next((i for i, h in enumerate(head) if h in v), None)
               for k, v in names.items()}
    has_head = (var is None or not HGVS.match(rows[0][var].replace(" ", ""))) and \
        (tx is None or not TX.match(rows[0][tx]))
    body = rows[1:] if has_head else rows
    out, skipped, given = [], 0, []
    if var is not None:
        how = f"colonnes transcrit {tx + 1 if tx is not None else '-'} et variant {var + 1}"
        # c./n. sans transcrit : chromosome et position (contenu ou nom de colonne)
        chc = psc = idx = None
        if tx is None:
            chc = vcfcols["chrom"] if vcfcols["chrom"] is not None else next(
                (j for j in range(width) if j != var and share(j, CHROM) > 0.8), None)
            psc = vcfcols["pos"] if vcfcols["pos"] is not None else next(
                (j for j in range(width) if j not in (var, chc) and share(j, POS) > 0.8), None)
            if chc is not None and psc is not None and a.gtf and a.gtf.exists():
                idx = mane_index(a.gtf)
                how += (f", chromosome {chc + 1} et position {psc + 1} : transcrit MANE Select "
                        "à cette position (GTF)")
        for r in body:
            v = r[var].replace(" ", "")
            if not v:
                skipped += 1
                continue
            if not G_HGVS.match(v) and ":" in v:
                v = v.split(":", 1)[1]                  # NM_...:c.xxx -> c.xxx
            t = r[tx].strip() if tx is not None and r[tx].strip() else "."
            if idx is not None and not G_HGVS.match(v) and CHROM.match(r[chc].strip()) \
                    and POS.match(r[psc].strip()):
                pos = int(float(r[psc]))
                t = pick_mane(idx, r[chc].strip(), pos, a.max_dist)
                given.append((t, v, r[chc].strip(), pos))
            out.append((t, v))
    elif all(i is not None for i in vcfcols.values()):
        idx = mane_index(a.gtf) if tx is None and a.gtf and a.gtf.exists() else None
        how = ("colonnes " + ", ".join(f"{k}={i + 1}" for k, i in vcfcols.items())
               + (f", transcrit {tx + 1}" if tx is not None else ", sans transcrit"))
        for r in body:
            try:
                v = g_notation(r[vcfcols["chrom"]], r[vcfcols["pos"]],
                               r[vcfcols["ref"]], r[vcfcols["alt"]])
            except (ValueError, IndexError):
                skipped += 1
                continue
            t = r[tx].strip() if tx is not None and r[tx].strip() else "."
            if t == "." and idx is not None:
                t = pick_mane(idx, r[vcfcols["chrom"]], int(float(r[vcfcols["pos"]])),
                              a.max_dist)
            out.append((t, v))
        if idx is not None:
            how += " ; transcrit MANE Select le plus proche (GTF)"
    else:
        sys.exit("aucune colonne de variant reconnue (notation c./n./g. ou chrom/pos/ref/alt) ; "
                 f"en-tête lu : {rows[0][:15]}")
    with open(a.sortie, "w") as fh:
        for t, v in out:
            fh.write(f"{t}\t{v}\n")
    if given:
        pf = a.sortie.with_suffix(".positions.tsv")
        with open(pf, "w") as fh:
            for t, v, c, pos in given:
                fh.write(f"{t}\t{v}\t{c}\t{pos}\n")
        print(f"  {len(given)} positions saisies gardées pour le contrôle -> {pf}")
    print(f"{a.entree.name} : {how} ; {len(out)} variants -> {a.sortie}"
          + (f" ({skipped} lignes vides ou illisibles ignorées)" if skipped else ""))
    n_dot = sum(t == "." for t, _v in out)
    if n_dot:
        print(f"  {n_dot} variants sans transcrit : la sélection du transcrit ne se fera pas "
              "pour eux (donner --gtf pour prendre le MANE Select voisin)")


if __name__ == "__main__":
    main()
