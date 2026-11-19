#!/usr/bin/env python3
"""Chemin 2 - preparation (machine connectee, p. ex. Curta).

Deux etapes, chacune hors donnees patient :

1. coordonnees : lit un GTF GENCODE (hg38) et ecrit, pour les ARN de la table de
   localisation, les coordonnees du gene et l'union des exons des transcrits.
       python3 preparer_references.py coordonnees --gtf gencode.gtf --genes genes.txt --out arn_coordonnees_hg38.tsv

2. sequences : extrait la sequence genomique (brin +) des regions a annoter
   (gene entier si <= 20 kb, sinon exons +/- 150 nt) depuis un FASTA hg38 local
   indexe (.fai) ou, a defaut, par requetes HTTP Range sur une URL de FASTA indexe.
       python3 preparer_references.py sequences --coord arn_coordonnees_hg38.tsv --fasta GRCh38.fa --out ref_regions_hg38.fa
       python3 preparer_references.py sequences --coord arn_coordonnees_hg38.tsv --fasta-url URL --fai fichier.fai --out ref_regions_hg38.fa

Aucun calcul sur variant ici. Aucune dependance hors bibliotheque standard (curl pour le mode URL).
"""
import argparse
import re
import subprocess
import sys

MAX_SPAN = 20000
FLANC_EXON = 150
FLANC_GENE = 100


def attrs(champ):
    return dict(re.findall(r'(\S+) "([^"]*)"', champ))


def cmd_coordonnees(a):
    genes = [l.strip() for l in open(a.genes) if l.strip() and not l.startswith("#")]
    voulus = set(genes)
    gene_rows, exons = {}, {}
    ouvrir = open
    if a.gtf.endswith(".gz"):
        import gzip
        ouvrir = gzip.open
    with ouvrir(a.gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("gene", "exon"):
                continue
            at = attrs(f[8])
            g = at.get("gene_name")
            if g not in voulus:
                continue
            gid = at.get("gene_id", "")
            cle = (g, gid)
            if f[2] == "gene":
                gene_rows[cle] = (f[0], int(f[3]), int(f[4]), f[6], at.get("gene_type", ""))
            else:
                exons.setdefault(cle, set()).add((int(f[3]), int(f[4])))
    with open(a.out, "w") as out:
        out.write("gene\tgene_id\tchr\tdebut\tfin\tbrin\tgene_type\texons_union\tsource\n")
        for g in genes:
            for (gn, gid), (c, s, e, st, gt) in sorted(gene_rows.items()):
                if gn != g:
                    continue
                ex = sorted(exons.get((gn, gid), []))
                fus = []
                for x, y in ex:
                    if fus and x <= fus[-1][1] + 1:
                        fus[-1][1] = max(fus[-1][1], y)
                    else:
                        fus.append([x, y])
                out.write(f"{gn}\t{gid}\t{c}\t{s}\t{e}\t{st}\t{gt}\t"
                          + ",".join(f"{x}-{y}" for x, y in fus) + f"\t{a.source}\n")
            if not any(k[0] == g for k in gene_rows):
                print(f"AVERTISSEMENT : {g} absent du GTF", file=sys.stderr)


def lire_fai(chemin):
    fai = {}
    for l in open(chemin):
        p = l.split("\t")
        fai[p[0]] = (int(p[1]), int(p[2]), int(p[3]), int(p[4]))
    return fai


def extraire(fai, lecteur, c, s, e):
    _, off, lb, lby = fai[c]

    def bo(p):
        return off + (p // lb) * lby + (p % lb)
    a, b = bo(s - 1), bo(e - 1)
    return lecteur(a, b).replace("\n", "").replace("\r", "").upper()


def regions(row):
    s, e = int(row["debut"]), int(row["fin"])
    if e - s + 1 <= MAX_SPAN or not row["exons_union"]:
        return [(s - FLANC_GENE, e + FLANC_GENE)]
    out = []
    for x in row["exons_union"].split(","):
        u, v = map(int, x.split("-"))
        u, v = u - FLANC_EXON, v + FLANC_EXON
        if out and u <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], v))
        else:
            out.append((u, v))
    return out


def cmd_sequences(a):
    if a.fasta:
        fai = lire_fai(a.fasta + ".fai")
        fh = open(a.fasta, "rb")

        def lecteur(x, y):
            fh.seek(x)
            return fh.read(y - x + 1).decode()
        provenance = a.fasta
    else:
        fai = lire_fai(a.fai)

        def lecteur(x, y):
            r = subprocess.run(["curl", "-sf", "-r", f"{x}-{y}", a.fasta_url],
                               capture_output=True, text=True, check=True)
            return r.stdout
        provenance = a.fasta_url
    rows = []
    with open(a.coord) as fh2:
        tete = fh2.readline().rstrip("\n").split("\t")
        for l in fh2:
            rows.append(dict(zip(tete, l.rstrip("\n").split("\t"))))
    with open(a.out, "w") as out:
        out.write(f"; provenance={provenance}\n")
        vus = set()
        for r in rows:
            for s, e in regions(r):
                cle = (r["chr"], s, e)
                if cle in vus:
                    continue
                vus.add(cle)
                sq = extraire(fai, lecteur, r["chr"], s, e)
                if len(sq) != e - s + 1:
                    sys.exit(f"longueur inattendue pour {cle}")
                out.write(f">{r['chr']}:{s}-{e} {r['gene']}\n")
                for i in range(0, len(sq), 80):
                    out.write(sq[i:i + 80] + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    c = sp.add_parser("coordonnees")
    c.add_argument("--gtf", required=True)
    c.add_argument("--genes", required=True)
    c.add_argument("--out", required=True)
    c.add_argument("--source", default="GENCODE")
    s = sp.add_parser("sequences")
    s.add_argument("--coord", required=True)
    s.add_argument("--fasta")
    s.add_argument("--fasta-url")
    s.add_argument("--fai")
    s.add_argument("--out", required=True)
    a = p.parse_args()
    if a.cmd == "coordonnees":
        cmd_coordonnees(a)
    else:
        if not a.fasta and not (a.fasta_url and a.fai):
            sys.exit("--fasta, ou --fasta-url avec --fai")
        cmd_sequences(a)


if __name__ == "__main__":
    main()
