#!/usr/bin/env python3
"""Préparation des tables du chemin 1 (machine connectée : Curta). Aucune donnée patient.

Sous-commandes (toutes écrivent des TSV, éventuellement gzippés) :
  mito         ensemble mito = MitoCarta3.0 (+ PanelApp), symboles HGNC à jour
  u12          met à jour les symboles HGNC de la table d'introns U12 ; table gène -> introns U12 ;
               avec --intronic : liste génome entier (intronIC) et part des gènes mito
  gtf          refFlat -> GTF minimal pour intronIC
  utr3         UTR 3' (transcrit codant à l'UTR 3' le plus long, par gène) avec séquence
  mirna        précurseurs miARN GENCODE + bras matures miRBase localisés dans la séquence génomique
  epissage     sites 5' (-3..+6), début d'intron (+1..+12) et fin d'intron (50 nt) de tous les introns codants
  elements     éléments ARN non codants (classes du chemin 1) avec séquence et chevauchements antisens
  partenaires  partenaires protéiques des ARN Pol III / ARNt / ARNr (groupes HGNC + ajouts manuels sourcés)
  snodb        (squelette) export snoDB -> snorna / cible_gene
  boucles      (squelette) boucles mitoloop -> ancres reliées à un gène mito
  coexpression (squelette) GTEx : r de Spearman lncRNA / module OXPHOS par tissu

Voir README.md pour l'ordre des commandes et les sources.
"""
import argparse
import csv
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from commun import (Fasta, Hgnc, adn, hypergeom_sf, introns, lire_liste_genes, lire_refflat, ouvrir,
                    revcomp, utr3_blocs)


def ecrire(path, header, rows):
    with (gzip.open(path, "wt") if path.endswith(".gz") else open(path, "w")) as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        n = 0
        for r in rows:
            w.writerow(r)
            n += 1
    print("%s : %d lignes" % (path, n), file=sys.stderr)


def seq_blocs(fa, chrom, strand, blocs):
    s = "".join(fa.seq(chrom, a, b) for a, b in sorted(blocs))
    return revcomp(s) if strand == "-" else s


# ---------------------------------------------------------------------- mito
def cmd_mito(a):
    h = Hgnc(a.hgnc)
    genes = {}
    for src, path in (("MitoCarta3.0", a.mitocarta), ("PanelApp", a.panelapp)):
        if not path:
            continue
        for g in lire_liste_genes(path, h):
            genes.setdefault(g, []).append(src)
    ecrire(a.out, ["gene", "sources"], ((g, "+".join(s)) for g, s in sorted(genes.items())))


# ----------------------------------------------------------------------- U12
def cmd_u12(a):
    h = Hgnc(a.hgnc)
    mito = lire_liste_genes(a.mito, h) if a.mito else set()
    rows, par_gene = [], {}
    with open(a.table) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        cols = rd.fieldnames
        for r in rd:
            g, st = h.maj(r["gene"])
            rows.append([g, r["gene"], st] + [r[c] for c in cols[1:]])
            d = par_gene.setdefault(g, {"introns": {}, "tx": set(), "orig": r["gene"], "statut": st})
            key = (r["chrom"], int(r["debut0"]), int(r["fin"]), r["brin"])
            d["introns"][key] = (r["sous_type"], r["debut_intron"], r.get("BPS_TCCTTRA", ""))
            d["tx"].add(r["transcrit"])
    ecrire(os.path.join(a.outdir, "introns_U12_genes_mito_hg38_hgnc.tsv"),
           ["gene_hgnc", "gene_table_source", "statut_hgnc"] + cols[1:], rows)

    def lignes():
        for g, d in sorted(par_gene.items()):
            intr = sorted(d["introns"])
            yield [g, d["orig"], d["statut"], "oui" if g in mito else ("non" if mito else "NA"), len(intr),
                   len(d["tx"]), ";".join("%s:%d-%d(%s)" % (c, s + 1, e, b) for c, s, e, b in intr),
                   ";".join(sorted({d["introns"][k][0] for k in intr})),
                   ";".join(d["introns"][k][1] for k in intr),
                   ";".join(d["introns"][k][2] or "." for k in intr)]
    ecrire(os.path.join(a.outdir, "gene_mito_introns_U12.tsv"),
           ["gene_hgnc", "gene_table_source", "statut_hgnc", "dans_ensemble_mito", "n_introns_U12",
            "n_transcrits", "introns_hg38_1based", "sous_types", "debut_intron_5ss", "BPS_TCCTTRA"], lignes())

    if a.intronic:
        cmd_u12_genome(a, h, mito)


def cmd_u12_genome(a, h, mito):
    """Lit la sortie intronIC v3 (<nom>.meta.iic : colonnes name, grandparent = gène, type_id = u12/u2 ;
    <nom>.bed.iic à côté pour les coordonnées). Les introns dupliqués entre transcrits sont déjà fusionnés
    par intronIC (dupe_map). Univers = gènes ayant au moins un intron classé."""
    bed = {}
    pbed = a.intronic.replace(".meta.iic", ".bed.iic")
    if os.path.exists(pbed):
        with ouvrir(pbed) as fh:
            for line in fh:
                p = line.rstrip("\n").split("\t")
                bed[p[3]] = (p[0], int(p[1]), int(p[2]), p[5])
    # introns identiques partagés entre transcrits (dupe_map) : un seul représentant, parfois d'un gène
    # de lecture continue ; on rattache l'intron à tous les gènes de ses transcrits (via --gtf).
    tx_gene, dupes = {}, {}
    if a.gtf:
        with ouvrir(a.gtf) as fh:
            for line in fh:
                m = re.search(r'gene_id "([^"]+)"; transcript_id "([^"]+)"', line)
                if m:
                    tx_gene[m.group(2)] = m.group(1)
        pd = a.intronic.replace(".meta.iic", ".dupe_map.iic")
        if os.path.exists(pd):
            with ouvrir(pd) as fh:
                for line in fh:
                    k, v = line.split()[:2]
                    dupes.setdefault(k, set()).add(v.split(":")[0])
    u12, tous = {}, set()
    with ouvrir(a.intronic) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes = {h.maj(r["grandparent"])[0]}
            m = re.match(r"[^@]+@([^_]+)_", r["name"])
            c = bed.get(r["name"])
            if m and c:
                cle = "%s:intron_%d_%d" % (m.group(1), c[1] + 1, c[2])
                genes |= {h.maj(tx_gene[t])[0] for t in dupes.get(cle, ()) if t in tx_gene}
            tous |= genes
            if r["type_id"] == "u12":
                for g in genes:
                    u12.setdefault(g, []).append((r["name"], r["dnts"], r["rel_score"]))
    N, K, n = len(tous), len(mito & tous), len(u12)
    k = len(set(u12) & mito)
    p = hypergeom_sf(k, N, K, n)

    def lignes(seulement_mito):
        for g, v in sorted(u12.items()):
            if seulement_mito and g not in mito:
                continue
            coords = []
            for name, dnts, sc in v:
                c = bed.get(name)
                coords.append("%s:%d-%d(%s)%s,score=%s" % (c[0], c[1] + 1, c[2], c[3], dnts, sc) if c
                              else "%s,%s,score=%s" % (name, dnts, sc))
            yield [g, "oui" if g in mito else "non", len(v), ";".join(coords)]
    hd = ["gene_hgnc", "dans_ensemble_mito", "n_introns_U12", "introns_hg38_1based_dnts_score"]
    ecrire(os.path.join(a.outdir, "genes_U12_genome_entier.tsv"), hd, lignes(False))
    ecrire(os.path.join(a.outdir, "gene_mito_introns_U12_intronIC.tsv"), hd, lignes(True))
    with open(os.path.join(a.outdir, "resume_U12_genome_entier.tsv"), "w") as fh:
        fh.write("genes_avec_introns\tgenes_mito_avec_introns\tgenes_U12\tgenes_U12_mito\t"
                 "part_mito_U12\tpart_mito_fond\tP_hypergeom_enrichissement\n")
        fh.write("%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.3g\n" % (N, K, n, k, k / n if n else 0, K / N if N else 0, p))
    print("U12 génome entier : %d gènes, dont %d mito (%.1f %% ; fond %.1f %% ; P = %.3g)"
          % (n, k, 100 * k / max(n, 1), 100 * K / max(N, 1), p), file=sys.stderr)


# ----------------------------------------------------------------------- GTF
def cmd_gtf(a):
    h = Hgnc(a.hgnc) if a.hgnc else None
    with open(a.out, "w") as fh:
        for tx in lire_refflat(a.refflat):
            if tx["cds_start"] >= tx["cds_end"] and not a.tous:
                continue
            g = h.maj(tx["gene"])[0] if h else tx["gene"]
            att = 'gene_id "%s"; transcript_id "%s";' % (g, tx["tx"])
            fh.write("%s\trefFlat\ttranscript\t%d\t%d\t.\t%s\t.\t%s\n" % (tx["chrom"], tx["start"] + 1, tx["end"],
                                                                      tx["strand"], att))
            for s, e in tx["exons"]:
                fh.write("%s\trefFlat\texon\t%d\t%d\t.\t%s\t.\t%s\n" % (tx["chrom"], s + 1, e, tx["strand"], att))
                cs, ce = max(s, tx["cds_start"]), min(e, tx["cds_end"])
                if cs < ce:
                    fh.write("%s\trefFlat\tCDS\t%d\t%d\t.\t%s\t.\t%s\n" % (tx["chrom"], cs + 1, ce, tx["strand"], att))


# ---------------------------------------------------------------------- UTR3
def cmd_utr3(a):
    h, fa = Hgnc(a.hgnc), Fasta(a.fasta)
    best = {}
    for tx in lire_refflat(a.refflat):
        if ("_" in tx["chrom"] or not tx["chrom"].startswith("chr")) or tx["chrom"] == "chrM" or tx["tx"].endswith("_PAR_Y"):
            continue
        g = h.maj(tx["gene"])[0]
        if h.locus_type(g) != "gene with protein product":
            continue
        bl = utr3_blocs(tx)
        L = sum(e - s for s, e in bl)
        if L >= a.min_len and (g not in best or L > best[g][0]):
            best[g] = (L, tx, bl)

    def lignes():
        for g, (L, tx, bl) in sorted(best.items()):
            yield [g, tx["tx"], tx["chrom"], tx["strand"], L,
                   ";".join("%d-%d" % b for b in sorted(bl)), seq_blocs(fa, tx["chrom"], tx["strand"], bl)]
    ecrire(a.out, ["gene", "transcrit", "chrom", "brin", "longueur", "blocs_0based", "sequence"], lignes())


# --------------------------------------------------------------------- miARN
def lire_fa(path, prefixe="hsa-"):
    d, name = {}, None
    with ouvrir(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                d[name] = []
            elif name:
                d[name].append(line.strip())
    return {k: adn("".join(v)) for k, v in d.items() if k.startswith(prefixe)}


def cmd_mirna(a):
    h, fa = Hgnc(a.hgnc), Fasta(a.fasta)
    hair, mat = lire_fa(a.hairpin), lire_fa(a.mature)
    vus = set()

    def lignes():
        for tx in lire_refflat(a.refflat):
            g = h.maj(tx["gene"])[0]
            if h.locus_type(g) != "RNA, micro" or ("_" in tx["chrom"] or not tx["chrom"].startswith("chr")):
                continue
            key = (tx["chrom"], tx["start"], tx["end"], tx["strand"])
            if key in vus:
                continue
            vus.add(key)
            F = a.flanc
            s0, e0 = tx["start"] - F, tx["end"] + F
            g_seq = fa.seq(tx["chrom"], s0, e0)
            seq = revcomp(g_seq) if tx["strand"] == "-" else g_seq
            core = seq[F:len(seq) - F]
            hp = sorted(n for n, s in hair.items() if core in s or s in seq)
            arms = []
            for n, m in mat.items():
                i = seq.find(m)
                if i >= 0:
                    arms.append((i, n, m))
            arms.sort()
            yield [g, tx["tx"], tx["chrom"], tx["strand"], s0, e0, ",".join(hp) or ".",
                   ";".join("%s:%d:%s" % (n, i, m) for i, n, m in arms) or ".", seq]
    ecrire(a.out, ["gene", "transcrit", "chrom", "brin", "debut0", "fin", "miRBase_hairpin",
                   "matures_nom:offset0:seq", "sequence_transcrit_flancs"], lignes())


# ------------------------------------------------------------------ épissage
def cmd_epissage(a):
    h, fa = Hgnc(a.hgnc), Fasta(a.fasta)
    vus = set()

    def lignes():
        for tx in lire_refflat(a.refflat):
            if tx["cds_start"] >= tx["cds_end"] or ("_" in tx["chrom"] or not tx["chrom"].startswith("chr")) or tx["chrom"] == "chrM":
                continue
            g = h.maj(tx["gene"])[0]
            if h.locus_type(g) != "gene with protein product":
                continue
            for s, e in introns(tx):
                key = (tx["chrom"], s, e, tx["strand"])
                if key in vus or e - s < 40:
                    continue
                vus.add(key)
                if tx["strand"] == "+":
                    ss5 = fa.seq(tx["chrom"], s - 3, s + 6)
                    deb = fa.seq(tx["chrom"], s, s + 12)
                    fin = fa.seq(tx["chrom"], e - 50, e)
                else:
                    ss5 = revcomp(fa.seq(tx["chrom"], e - 6, e + 3))
                    deb = revcomp(fa.seq(tx["chrom"], e - 12, e))
                    fin = revcomp(fa.seq(tx["chrom"], s, s + 50))
                yield [g, tx["chrom"], tx["strand"], s, e, ss5, deb, fin]
    ecrire(a.out, ["gene", "chrom", "brin", "debut0", "fin", "ss5_m3_p6", "intron_p1_p12", "intron_fin50"],
           lignes())


# ------------------------------------------------------------------ éléments
CLASSES = {"RNA, micro": "miARN", "RNA, small nuclear": "snRNA", "RNA, small nucleolar": "snoRNA",
           "RNA, long non-coding": "lncRNA", "RNA, ribosomal": "ARNr", "RNA, transfer": "ARNt",
           "RNA, vault": "ARN_voute", "RNA, Y": "ARN_Y", "RNA, misc": "ARN_misc"}
SOUS_SNRNA = [(r"^RNU4ATAC", "U4atac"), (r"^RNU6ATAC", "U6atac"), (r"^RNU11$", "U11"), (r"^RNU12$", "U12"),
              (r"^(RNU1-|RNVU1-)", "U1"), (r"^RNU2-", "U2"), (r"^RNU4-", "U4"), (r"^RNU5", "U5"),
              (r"^RNU6-", "U6"), (r"^RNU7-", "U7")]
POL3 = re.compile(r"^(RMRP|RPPH1|RN7SK|TERC|RN7SL\d|RNY\d|VTRNA|RNA5S|RNA45S|RNA5-8S|RNA18S|RNA28S)")


def sous_classe(g, cl):
    if cl == "snRNA":
        for rx, s in SOUS_SNRNA:
            if re.match(rx, g):
                return s
        return "autre"
    if POL3.match(g):
        return "PolIII_" + re.match(r"[A-Z]+\d?[A-Z]*", g).group(0)
    return cl


def cmd_elements(a):
    h, fa = Hgnc(a.hgnc), Fasta(a.fasta)
    codants, elems = [], {}
    for tx in lire_refflat(a.refflat):
        if ("_" in tx["chrom"] or not tx["chrom"].startswith("chr")) or tx["chrom"] == "chrM":
            continue
        g = h.maj(tx["gene"])[0]
        lt = h.locus_type(g)
        if lt == "gene with protein product" and tx["cds_start"] < tx["cds_end"]:
            codants.append((tx["chrom"], tx["start"], tx["end"], tx["strand"], g))
        cl = CLASSES.get(lt)
        if not cl:
            continue
        # un élément par locus : transcrits du même symbole, même brin, distants de moins de 10 kb
        lst = elems.setdefault((g, tx["chrom"]), [])
        e = next((x for x in lst if x["strand"] == tx["strand"] and tx["start"] < x["end"] + 10000
                  and x["start"] - 10000 < tx["end"]), None)
        if e is None:
            e = {"g": g, "cl": cl, "chrom": tx["chrom"], "strand": tx["strand"], "start": tx["start"],
                 "end": tx["end"], "tx": set()}
            lst.append(e)
        e["start"], e["end"] = min(e["start"], tx["start"]), max(e["end"], tx["end"])
        e["tx"].add(tx["tx"])
    # chevauchements antisens (lncRNA <-> gène codant, brin opposé)
    par_chr = {}
    for c in codants:
        par_chr.setdefault(c[0], []).append(c)

    def lignes():
        for e in sorted((x for lst in elems.values() for x in lst), key=lambda x: (x["chrom"], x["start"])):
            g, chrom = e["g"], e["chrom"]
            L = e["end"] - e["start"]
            seq = fa.seq(chrom, e["start"], e["end"]) if L <= a.max_seq else ""
            anti = sorted({c[4] for c in par_chr.get(chrom, []) if c[1] < e["end"] and e["start"] < c[2]
                           and c[3] != e["strand"]}) if e["cl"] == "lncRNA" else []
            yield [g, e["cl"], sous_classe(g, e["cl"]), chrom, e["strand"], e["start"], e["end"],
                   ",".join(sorted(e["tx"])), ",".join(anti) or ".", seq or "."]
    ecrire(a.out, ["gene", "classe", "sous_classe", "chrom", "brin", "debut0", "fin", "transcrits",
                   "antisens_chevauchant_codant", "sequence_brin_plus"], lignes())


# --------------------------------------------------------------- partenaires
# Groupes HGNC (fait vérifiable : appartenance au groupe dans hgnc_complete_set) par classe d'ARN.
GROUPES = {
    "RMRP": ["RNase MRP complex subunits"],
    "RPPH1": ["RNase MRP complex subunits"],
    "TERC": ["Telomerase holoenzyme subunits", "H/ACA ribonucleoprotein complex"],
    "RN7SK": ["7SK snRNP complex"],
    "RN7SL": ["Signal recognition particle", "Signal recognition particle GTPases"],
    "ARNt": ["@tRNA"],  # tous les groupes HGNC dont le nom contient « tRNA »
}
# Ajouts manuels : uniquement des faits vérifiés (référence avec DOI dans README.md).
MANUELS = [
    ("RMRP", "PNPT1", "import mitochondrial de RMRP dépendant de PNPase", "Wang et al. 2010 Cell doi:10.1016/j.cell.2010.06.035"),
    ("RPPH1", "PNPT1", "import mitochondrial de l'ARN de la RNase P dépendant de PNPase", "Wang et al. 2010 Cell doi:10.1016/j.cell.2010.06.035"),
    ("RNA5S", "PNPT1", "import mitochondrial de l'ARNr 5S dépendant de PNPase", "Wang et al. 2010 Cell doi:10.1016/j.cell.2010.06.035"),
]


def cmd_partenaires(a):
    h = Hgnc(a.hgnc)
    mito = {}
    with ouvrir(a.mito) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mito[r["gene"]] = r.get("sources", "mito")
    groupes = {}
    for s, r in h.info.items():
        for gname in (r.get("gene_group") or "").strip('"').split("|"):
            if gname and h.locus_type(s) == "gene with protein product":
                groupes.setdefault(gname, set()).add(s)
    rows = []
    for arn, gl in GROUPES.items():
        noms = sorted(n for n in groupes if "tRNA" in n) if gl == ["@tRNA"] else gl
        for gname in noms:
            for p in sorted(groupes.get(gname, ())):
                rows.append([arn, p, mito.get(p, "non"), cible_nucleaire(h, p, mito),
                             "groupe HGNC : " + gname, "HGNC hgnc_complete_set"])
    for arn, p, role, ref in MANUELS:
        rows.append([arn, p, mito.get(p, "non"), "oui (" + role + ")", role, ref])
    ecrire(a.out, ["arn", "partenaire", "ensemble_mito", "agit_sur_arn_nucleaire", "lien", "source"], rows)


def cible_nucleaire(h, p, mito):
    """Le partenaire agit-il sur l'ARN nucléaire (et pas seulement sur son homologue mitochondrial) ?
    Règle prudente : nom HGNC contenant « mitochondrial » -> non ; sinon, présent dans l'ensemble mito
    -> a_verifier (double localisation possible, à confirmer dans la littérature) ; sinon oui."""
    nom = h.info.get(p, {}).get("name", "")
    if "mitochondrial" in nom:
        return "non (paralogue mitochondrial : " + nom + ")"
    if p in mito:
        return "a_verifier (gene de l'ensemble mito sans mention mitochondrial dans le nom HGNC)"
    return "oui"


# ------------------------------------------------- squelettes snoRNA / lncRNA
def cmd_snodb(a):
    """Export snoDB (ou autre table d'interactions snoRNA-cible) -> snorna / cible_gene, symboles HGNC.
    Les noms de colonnes de l'export ne sont pas figés ici : les passer en option après inspection."""
    h = Hgnc(a.hgnc)
    vus = set()
    with ouvrir(a.table) as fh:
        for r in csv.DictReader(fh, delimiter=a.sep):
            s, c = r[a.col_snorna].strip(), r[a.col_cible].strip()
            if s and c:
                vus.add((h.maj(s)[0], h.maj(c)[0]))
    ecrire(a.out, ["snorna", "cible_gene"], sorted(vus))


def cmd_boucles(a):
    """Boucles mitoloop (ex. /scratch/.../MITODARK/loops/loops.tsv) -> gene_mito / chrom / debut0 / fin
    de l'ancre distale. Colonnes à préciser en option (format de loops.tsv non figé ici)."""
    h = Hgnc(a.hgnc)
    with ouvrir(a.table) as fh:
        rows = [[h.maj(r[a.col_gene])[0], r[a.col_chrom], int(float(r[a.col_debut])), int(float(r[a.col_fin]))]
                for r in csv.DictReader(fh, delimiter="\t")]
    ecrire(a.out, ["gene_mito", "chrom", "debut0", "fin"], rows)


def cmd_coexpression(a):
    """Co-expression GTEx lncRNA / module OXPHOS (MitoCarta3.0, MitoPathway « OXPHOS subunits »).
    Par tissu : score du module = médiane des z-scores log2(TPM+1) des gènes du module ; r de Spearman
    entre chaque lncRNA (TPM médian >= --min-tpm) et ce score. Écrit lncrna / module / tissu / r / n."""
    import numpy as np
    h = Hgnc(a.hgnc)
    module = set()
    with ouvrir(a.mitocarta) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if "OXPHOS subunits" in (r.get("MitoCarta3.0_MitoPathways") or ""):
                module.add(h.maj(r["Symbol"])[0])
    tissu_de = {}
    with ouvrir(a.attributs) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            tissu_de[r["SAMPID"]] = r["SMTSD"]
    tissus = set(a.tissus.split(","))
    with ouvrir(a.gct) as fh:
        fh.readline(); fh.readline()
        cols = fh.readline().rstrip("\n").split("\t")
        idx = {t: [i for i, c in enumerate(cols) if tissu_de.get(c) == t] for t in tissus}
        lnc, mod = {}, {}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            g, lt = h.maj(p[1])[0], None
            lt = h.locus_type(g)
            if g in module or lt == "RNA, long non-coding":
                v = np.array(p[2:], dtype=float)
                (mod if g in module else lnc)[g] = v
    def rang(x):
        return np.argsort(np.argsort(x)).astype(float)
    rows = []
    for t, ii in idx.items():
        if len(ii) < a.min_n:
            continue
        ii = [i - 2 for i in ii]
        M = np.array([np.log2(v[ii] + 1) for v in mod.values()])
        Z = (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9)
        score = rang(np.median(Z, 0))
        for g, v in lnc.items():
            x = v[ii]
            if np.median(x) < a.min_tpm:
                continue
            r = np.corrcoef(rang(x), score)[0, 1]
            rows.append([g, "OXPHOS_subunits", t, "%.3f" % r, len(ii)])
    ecrire(a.out, ["lncrna", "module", "tissu", "r", "n_echantillons"], rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    p = sp.add_parser("mito"); p.add_argument("--hgnc", required=True); p.add_argument("--mitocarta", required=True)
    p.add_argument("--panelapp"); p.add_argument("--out", required=True)
    p = sp.add_parser("u12"); p.add_argument("--hgnc", required=True); p.add_argument("--table", required=True)
    p.add_argument("--mito"); p.add_argument("--gtf", help="GTF donné à intronIC (rattache les introns dupliqués à leurs gènes)")
    p.add_argument("--intronic", help="<nom>.meta.iic de intronIC (seuil par défaut 90 %%)")
    p.add_argument("--outdir", required=True)
    p = sp.add_parser("gtf"); p.add_argument("--refflat", required=True); p.add_argument("--hgnc")
    p.add_argument("--tous", action="store_true"); p.add_argument("--out", required=True)
    p = sp.add_parser("utr3"); p.add_argument("--refflat", required=True); p.add_argument("--fasta", required=True)
    p.add_argument("--hgnc", required=True); p.add_argument("--min-len", type=int, default=1)
    p.add_argument("--out", required=True)
    p = sp.add_parser("mirna"); p.add_argument("--refflat", required=True); p.add_argument("--fasta", required=True)
    p.add_argument("--hgnc", required=True); p.add_argument("--hairpin", required=True)
    p.add_argument("--mature", required=True); p.add_argument("--flanc", type=int, default=15)
    p.add_argument("--out", required=True)
    p = sp.add_parser("epissage"); p.add_argument("--refflat", required=True); p.add_argument("--fasta", required=True)
    p.add_argument("--hgnc", required=True); p.add_argument("--out", required=True)
    p = sp.add_parser("elements"); p.add_argument("--refflat", required=True); p.add_argument("--fasta", required=True)
    p.add_argument("--hgnc", required=True); p.add_argument("--max-seq", type=int, default=2000)
    p.add_argument("--out", required=True)
    p = sp.add_parser("partenaires"); p.add_argument("--hgnc", required=True); p.add_argument("--mito", required=True)
    p.add_argument("--out", required=True)
    p = sp.add_parser("snodb"); p.add_argument("--hgnc", required=True); p.add_argument("--table", required=True)
    p.add_argument("--col-snorna", required=True); p.add_argument("--col-cible", required=True)
    p.add_argument("--sep", default="\t"); p.add_argument("--out", required=True)
    p = sp.add_parser("boucles"); p.add_argument("--hgnc", required=True); p.add_argument("--table", required=True)
    for c in ("gene", "chrom", "debut", "fin"):
        p.add_argument("--col-" + c, required=True)
    p.add_argument("--out", required=True)
    p = sp.add_parser("coexpression"); p.add_argument("--hgnc", required=True); p.add_argument("--gct", required=True)
    p.add_argument("--attributs", required=True); p.add_argument("--mitocarta", required=True)
    p.add_argument("--tissus", default="Muscle - Skeletal,Heart - Left Ventricle,Brain - Cortex,Kidney - Cortex")
    p.add_argument("--min-tpm", type=float, default=1.0); p.add_argument("--min-n", type=int, default=30)
    p.add_argument("--out", required=True)
    a = ap.parse_args()
    globals()["cmd_" + a.cmd](a)


if __name__ == "__main__":
    main()
