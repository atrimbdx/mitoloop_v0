#!/usr/bin/env python3
"""
preparer_locus.py - prepares the polymer simulation (simpol, path A: loop
extrusion stopped by oriented CTCF sites) of one MitoLoop locus, and the
in-silico deletions it is asked about.

Path A is the orthogonal check of AlphaGenome 2: its only input is the
position and orientation of the CTCF sites; it never sees the contact map
nor the rE2G elements. Path B of simpol (binders fitted on the map) is
circular here and is not used.

Subcommands
  preparer   locus window, anchors, measured Hi-C matrix and CTCF beds:
               loci.json      the locus, read by simpol/config.py through
                              SIMPOL_LOCI: window (TSS, loop summits and
                              elements of the gene, +/- --flanc), anchors
                              (promoter TSS +/- 2.5 kb, each loop summit
                              +/- --ancre-demi, each element on an anchor
                              of the tissue), hub = promoter + elements
               matrice.tsv    the measured contacts of the window (startA,
                              endA, startB, endB, N), read from a .hic with
                              hicstraw; used only to choose the parameter
                              set of path A, never as an input of it
               ctcf_SAUVAGE.bed      oriented CTCF sites (simpol
                                     ctcf_encode.py on --ctcf-pics)
               ctcf_SANSELEMENT.bed  the same without the sites inside the
                                     elements on anchor
               ctcf_SANSANCRE.bed    without the sites within --ancre-demi
                                     of the gene's distal loop summits
               plan_ajustement.tsv   path A, condition SAUVAGE, every
                                     parameter set, --replicats-ajustement
             A deletion that removes no CTCF site leaves path A with the
             same input: its differential is then null by construction, and
             the script says so (the element is not architectural under the
             CTCF-extrusion hypothesis).
  plan       plan_deletion.tsv: the chosen parameter set (--jeu, or the best
             path-A set of comparaison_SAUVAGE.tsv), every condition,
             --replicats replicates

Usage (mitoloop_v0, mitodark environment for hicstraw)
  python3 polymere/preparer_locus.py preparer --gene ABAT --type neurone \
      --cibles ~/reannot_git/results/ag2_ABAT/cibles.tsv \
      --loops results/loops_tissus/neurone/loops.tsv --hic <URL or .hic> \
      --ctcf-pics <ENCODE CTCF peaks, URL or bed.gz> --genome data/GRCh38.fa \
      --out polymere/abat
  python3 polymere/preparer_locus.py plan --dir polymere/abat
"""

import argparse
import csv
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

CONDITIONS = ("SAUVAGE", "SANSELEMENT", "SANSANCRE")
N_JEUX_A = 18          # simpol config.GRILLE_A : 3 x 3 x 2


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def matrice_hic(hic, chrom, debut, fin, res, norm, out):
    import hicstraw
    h = hicstraw.HiCFile(hic)
    noms = {c.name for c in h.getChromosomes()}
    c = chrom if chrom in noms else (chrom[3:] if chrom.startswith("chr") else "chr" + chrom)
    if c not in noms:
        sys.exit(f"{chrom} absent de {hic}")
    for nm in (norm, "NONE"):
        try:
            mzd = h.getMatrixZoomData(c, c, "observed", nm, "BP", res)
            recs = mzd.getRecords(debut, fin - 1, debut, fin - 1)
        except Exception as exc:              # normalisation absente : comptes bruts
            print(f"  {nm} : {exc}", file=sys.stderr)
            continue
        if recs:
            if nm != norm:
                print(f"  {norm} indisponible, comptes bruts (NONE)", file=sys.stderr)
            break
    else:
        sys.exit("matrice vide : vérifier la fenêtre et la résolution")
    n = 0
    with open(out, "w") as fo:
        fo.write("startA\tendA\tstartB\tendB\tN\n")
        for r in recs:
            if r.counts and r.counts == r.counts:          # NaN écarté
                fo.write(f"{r.binX}\t{r.binX + res}\t{r.binY}\t{r.binY + res}\t{r.counts:.4f}\n")
                n += 1
    return n


def lire_bed(path):
    with open(path) as fh:
        return [ln.rstrip("\n").split("\t") for ln in fh if ln.strip() and not ln.startswith("#")]


def ecrire_bed(rows, path):
    with open(path, "w") as fo:
        for r in rows:
            fo.write("\t".join(r) + "\n")


def cmd_preparer(a):
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    cib = [r for r in read_tsv(a.cibles) if r["gene"] == a.gene]
    promo = next((r for r in cib if r["classe"] == "promoteur"), None)
    if promo is None:
        sys.exit(f"pas de promoteur de {a.gene} dans {a.cibles}")
    chrom, tss = promo["chrom"], int(promo["tss"])
    elems = [r for r in cib if r["classe"] == "element_sur_ancre" and r["type"] == a.type]
    loops = [r for r in read_tsv(a.loops) if r["gene"] == a.gene and r["chrom"] == chrom]
    summits = sorted({int(float(r["summit"])) for r in loops})
    if not summits:
        sys.exit(f"aucune boucle de {a.gene} dans {a.loops}")
    pos = [tss] + summits + [int(e["debut"]) for e in elems] + [int(e["fin"]) for e in elems]
    debut = (min(pos) - a.flanc) // a.bin * a.bin
    fin = -(-(max(pos) + a.flanc) // a.bin) * a.bin
    ancres = [[f"promoteur_{a.gene}", tss - 2500, tss + 2500]]
    for s in summits:
        ancres.append([f"ancre_{s // 1000}kb", s - a.ancre_demi, s + a.ancre_demi])
    for i, e in enumerate(elems, 1):
        ancres.append([f"element_{i}", int(e["debut"]), int(e["fin"])])
    hub = [ancres[0][0]] + ([f"element_{i}" for i in range(1, len(elems) + 1)] or [ancres[1][0]])
    locus = {a.nom: dict(chrom=chrom, debut=debut, fin=fin, taille_bin=a.bin, marge=a.marge,
                         ancres=ancres, hub=hub, hub_fibro=[], gene=a.gene, type=a.type)}
    (out / "loci.json").write_text(json.dumps(locus, indent=2))
    print(f"locus {a.nom} : {chrom}:{debut:,}-{fin:,} (+/- {a.marge:,} simulés), billes de {a.bin} pb, "
          f"{(fin - debut + 2 * a.marge) // a.bin} billes")
    for nm, s, e in ancres:
        print(f"  ancre {nm:22s} {s:,}-{e:,}")
    print(f"  hub : {', '.join(hub)}")

    # matrice mesurée (choix du jeu de paramètres seulement)
    if a.hic:
        n = matrice_hic(a.hic, chrom, max(0, debut - a.marge), fin + a.marge, a.res, a.norm,
                        out / "matrice.tsv")
        print(f"matrice : {n:,} cases ({a.res} pb) -> {out / 'matrice.tsv'}")

    # sites CTCF orientés, puis les conditions de délétion
    pics = a.ctcf_pics
    if not pics and a.ctcf_biosample:
        # premier fichier de pics GRCh38 publié pour ce terme ENCODE (simpol ctcf_encode.py)
        sys.path.insert(0, str(Path(__file__).resolve().parent / "simpol"))
        from ctcf_encode import chercher_encode
        acc, pics = chercher_encode(a.ctcf_biosample)
        print(f"CTCF : ENCODE {acc} ({a.ctcf_biosample}) {pics}")
        (out / "ctcf_source.txt").write_text(f"{acc}\t{a.ctcf_biosample}\t{pics}\n")
    if not pics:
        sys.exit("--ctcf-pics ou --ctcf-biosample requis")
    if pics.startswith(("http://", "https://")):
        loc = out / Path(pics.split("?")[0]).name
        if not loc.exists():
            print(f"téléchargement {pics}", file=sys.stderr)
            urllib.request.urlretrieve(pics, loc)
        pics = str(loc)
    simpol = Path(__file__).resolve().parent / "simpol"
    sauvage = out / "ctcf_SAUVAGE.bed"
    subprocess.run([sys.executable, str(simpol / "ctcf_encode.py"), "--pics", pics,
                    "--genome", a.genome, "--chrom", chrom, "--debut", str(debut - a.marge),
                    "--fin", str(fin + a.marge), "--sortie", str(sauvage)], check=True)
    sites = lire_bed(sauvage)

    def centre(r):
        return (int(r[1]) + int(r[2])) / 2
    zones = {
        "SANSELEMENT": [(int(e["debut"]), int(e["fin"])) for e in elems],
        "SANSANCRE": [(s - a.ancre_demi, s + a.ancre_demi) for s in summits
                      if abs(s - tss) > a.ancre_demi],
    }
    bilan = []
    for cond, zs in zones.items():
        garde = [r for r in sites if not any(z0 <= centre(r) <= z1 for z0, z1 in zs)]
        ecrire_bed(garde, out / f"ctcf_{cond}.bed")
        retires = len(sites) - len(garde)
        bilan.append((cond, retires))
        msg = f"  {cond:12s} {retires} site(s) CTCF retiré(s) sur {len(sites)}"
        if not retires:
            msg += ("  -> même entrée que SAUVAGE : la voie A ne peut prédire aucun changement "
                    "(différentiel nul par construction)")
        print(msg)
    with open(out / "conditions.tsv", "w") as fo:
        fo.write("condition\tsites_retires\n")
        fo.write("SAUVAGE\t0\n")
        for c, n in bilan:
            fo.write(f"{c}\t{n}\n")
    with open(out / "plan_ajustement.tsv", "w") as fo:
        fo.write("condition\tvoie\tjeu\treplicat\n")
        for j in range(N_JEUX_A):
            for r in range(a.replicats_ajustement):
                fo.write(f"SAUVAGE\tA\t{j}\t{r}\n")
    print(f"plan_ajustement.tsv : {N_JEUX_A * a.replicats_ajustement} trajectoires "
          f"(SAUVAGE, {N_JEUX_A} jeux x {a.replicats_ajustement} réplicats)")


def cmd_plan(a):
    d = Path(a.dir)
    jeu = a.jeu
    if jeu is None:
        comp = d / "comparaison_SAUVAGE.tsv"
        if not comp.exists():
            sys.exit(f"{comp} absent : lancer comparer.py sur SAUVAGE, ou donner --jeu")
        rows = [r for r in read_tsv(comp) if r["voie"] == "A"]
        jeu = int(rows[0]["jeu"])                 # trié par accord des profils, décroissant
        print(f"jeu retenu : {jeu} (profils r={float(rows[0]['profils_pearson_moyen']):.3f}, "
              f"{rows[0]['parametres']})")
    retires = {r["condition"]: int(r["sites_retires"]) for r in read_tsv(d / "conditions.tsv")}
    conds = [c for c in CONDITIONS if c == "SAUVAGE" or retires.get(c, 0) > 0]
    for c in CONDITIONS:
        if c not in conds:
            print(f"  {c} : aucun site retiré, pas simulé (nul par construction)")
    n = 0
    with open(d / "plan_deletion.tsv", "w") as fo:
        fo.write("condition\tvoie\tjeu\treplicat\n")
        for c in conds:
            for r in range(a.replicats):
                fo.write(f"{c}\tA\t{jeu}\t{r}\n")
                n += 1
    (d / "jeu_retenu.txt").write_text(f"{jeu}\n")
    print(f"plan_deletion.tsv : {n} trajectoires ({', '.join(conds)} ; jeu {jeu} ; "
          f"{a.replicats} réplicats)")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd")
    s = sub.add_parser("preparer")
    s.add_argument("--gene", required=True)
    s.add_argument("--type", required=True, help="type cellulaire des éléments (ex. neurone)")
    s.add_argument("--cibles", required=True, help="cibles.tsv de ag2_balayage.py locus")
    s.add_argument("--loops", required=True, help="loops.tsv du tissu")
    s.add_argument("--hic", default="", help=".hic ou URL : matrice mesurée (choix du jeu)")
    s.add_argument("--res", type=int, default=10000)
    s.add_argument("--norm", default="SCALE")
    s.add_argument("--ctcf-pics", default="", help="pics CTCF ENCODE (URL ou bed[.gz])")
    s.add_argument("--ctcf-biosample", default="",
                   help="terme ENCODE : le fichier de pics CTCF est cherché sur le portail")
    s.add_argument("--genome", required=True)
    s.add_argument("--nom", default="", help="nom du locus (défaut : gène en minuscules)")
    s.add_argument("--bin", type=int, default=5000)
    s.add_argument("--marge", type=int, default=300000, help="flanc simulé de part et d'autre")
    s.add_argument("--flanc", type=int, default=150000, help="fenêtre analysée autour des ancres")
    s.add_argument("--ancre-demi", type=int, default=5000)
    s.add_argument("--replicats-ajustement", type=int, default=4)
    s.add_argument("--out", required=True)
    q = sub.add_parser("plan")
    q.add_argument("--dir", required=True)
    q.add_argument("--jeu", type=int, default=None)
    q.add_argument("--replicats", type=int, default=20)
    a = p.parse_args()
    if a.cmd == "preparer":
        a.nom = a.nom or a.gene.lower()
        cmd_preparer(a)
    elif a.cmd == "plan":
        cmd_plan(a)
    else:
        p.error("sous-commande attendue : preparer ou plan")


if __name__ == "__main__":
    main()
