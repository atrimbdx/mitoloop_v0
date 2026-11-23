#!/usr/bin/env python3
"""Chemin 3, préparation (machine connectée, Curta).

Construit la table sORF -> niveau mito, un GTF personnalisé pour VEP (--gtf) et un
BED (--custom), à partir :
  - des sORF traduits : catalogue GENCODE des Ribo-seq ORFs (BED12) ou fichiers
    d'entrée de la phase I (dépôt jorruior/gencode-riboseqORFs), OpenProt (BED12) ;
  - de la table manuelle des micropeptides mitochondriaux établis ;
  - des prédictions TargetP 2.0, MitoFates, TMbed ou DeepTMHMM si elles existent ;
  - de MitoCarta3.0 et de listes de peptides détectés dans des protéomes
    mitochondriaux publiés.

Les sORF de l'ADNmt (chrM : humanine, MOTS-c, SHLP) sont exclus.
Seuils et règles : voir sorf_lib.py et README.md (fixés a priori).
"""
import argparse
import csv
import hashlib
import re
import sys
from collections import defaultdict

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
import sorf_lib as L  # noqa: E402


def lire_ensemble(chemin, colonnes):
    """Fichier texte (une valeur par ligne) ou TSV avec une des colonnes données."""
    vals = set()
    with L.ouvrir(chemin) as f:
        lignes = [l.rstrip("\n") for l in f if l.strip() and not l.startswith("#")]
    if not lignes:
        return vals
    tete = lignes[0].split("\t")
    idx = [i for i, c in enumerate(tete) if c in colonnes]
    if idx:
        for l in lignes[1:]:
            c = l.split("\t")
            for i in idx:
                if i < len(c) and c[i].strip():
                    vals.add(c[i].strip().rstrip("*"))
    else:
        vals = {l.split("\t")[0].strip().rstrip("*") for l in lignes}
    return vals


def lire_targetp(chemin):
    """Sortie résumée de TargetP 2.0 (« *_summary.targetp2 ») : ID, Prediction, ..."""
    r = {}
    with L.ouvrir(chemin) as f:
        for l in f:
            if l.startswith("#") or not l.strip():
                continue
            c = l.rstrip("\n").split("\t")
            if len(c) >= 2:
                r[c[0].split()[0]] = c[1].strip()
    return r


def lire_mitofates(chemin):
    """Sortie de MitoFates : colonnes « Sequence ID », « Probability of presequence »,
    « Prediction ». On garde la décision par défaut de l'outil."""
    r = {}
    with L.ouvrir(chemin) as f:
        lec = csv.reader(f, delimiter="\t")
        tete = next(lec)
        try:
            ip = [i for i, c in enumerate(tete) if c.strip().lower() == "prediction"][0]
        except IndexError:
            sys.exit("MitoFates : colonne « Prediction » introuvable")
        for c in lec:
            if c:
                r[c[0].split()[0]] = "oui" if c[ip].strip().lower().startswith("possessing") else "non"
    return r


def lire_tm_3lignes(chemin):
    """Format 3 lignes (TMbed, DeepTMHMM) : >id [| classe], séquence, étiquettes.
    Hélice TM : étiquettes H/h (TMbed) ou M (DeepTMHMM), ou classe « TM »/« SP+TM »."""
    r, lignes = {}, [l.rstrip("\n") for l in L.ouvrir(chemin) if l.strip()]
    i = 0
    while i < len(lignes):
        if lignes[i].startswith(">"):
            tete = lignes[i][1:]
            oid = tete.split()[0].split("|")[0].strip()
            classe = tete.split("|")[1].strip() if "|" in tete else ""
            etiq = lignes[i + 2] if i + 2 < len(lignes) else ""
            tm = classe in ("TM", "SP+TM") or bool(re.search("[HhM]", etiq))
            r[oid] = "oui" if tm else "non"
            i += 3
        else:
            i += 1
    return r


def lire_mitocarta(chemin):
    """MitoCarta3.0 exportée en TSV (feuille « A Human MitoCarta3.0 ») : colonne
    Symbol, et Synonyms si présente."""
    s = set()
    with L.ouvrir(chemin) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("Symbol"):
                s.add(r["Symbol"].strip())
            for x in (r.get("Synonyms") or "").replace("|", ",").split(","):
                if x.strip() and x.strip() != "-":
                    s.add(x.strip())
    return s


def gtf_cds_canoniques(gtf, symboles):
    """CDS (+ stop) du transcrit Ensembl_canonical des gènes demandés (GTF GENCODE)."""
    par_tx = defaultdict(lambda: {"exons": [], "gene": "", "chrom": "", "strand": ""})
    with L.ouvrir(gtf) as f:
        for l in f:
            if l.startswith("#"):
                continue
            c = l.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] not in ("CDS", "stop_codon"):
                continue
            g = re.search(r'gene_name "([^"]+)"', c[8])
            if not g or g.group(1) not in symboles or "Ensembl_canonical" not in c[8]:
                continue
            tx = re.search(r'transcript_id "([^"]+)"', c[8]).group(1)
            d = par_tx[tx]
            d.update(gene=g.group(1), chrom=L.norm_chrom(c[0]), strand=c[6])
            d["exons"].append((int(c[3]) - 1, int(c[4])))
    out = []
    for tx, d in par_tx.items():
        ex = fusion(sorted(d["exons"]))
        out.append({"orf_id": "%s|%s" % (d["gene"], tx), "source": "gencode_cds_canonique",
                    "etude": "", "gene_hote": d["gene"], "chrom": d["chrom"],
                    "strand": d["strand"], "exons": ex, "peptide": ""})
    return out


def hote_par_chevauchement(orfs, gtf):
    """sORF sans gène hôte (OpenProt) : gène du GTF qui contient son 1er exon, même brin d'abord, sinon brin
    opposé (noté « antisens:GENE ») ; le plus court s'il y en a plusieurs."""
    sans = [o for o in orfs if not o.get("gene_hote")]
    if not sans:
        return 0
    genes = defaultdict(list)
    with L.ouvrir(gtf) as f:
        for l in f:
            if l.startswith("#"):
                continue
            c = l.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            g = re.search(r'gene_name "([^"]+)"', c[8])
            if g:
                genes[L.norm_chrom(c[0])].append((int(c[3]) - 1, int(c[4]), c[6], g.group(1)))
    n = 0
    for o in sans:
        a, b = o["exons"][0]
        cands = [x for x in genes.get(o["chrom"], ()) if x[0] <= a and b <= x[1]]
        meme = sorted((x for x in cands if x[2] == o["strand"]), key=lambda x: x[1] - x[0])
        opp = sorted((x for x in cands if x[2] != o["strand"]), key=lambda x: x[1] - x[0])
        if meme:
            o["gene_hote"] = meme[0][3]
        elif opp:
            o["gene_hote"] = "antisens:" + opp[0][3]
        else:
            continue
        n += 1
    return n


def fusion(ex):
    out = []
    for a, b in ex:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def lire_manuel(chemin):
    with open(chemin) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def apparier_manuel(orfs, manuel, log):
    """Règles d'appariement de la colonne « appariement » :
    symbole:A|B   -> sORF dont le gène hôte est A ou B (CDS GENCODE ou Ribo-seq) ;
    hote:X;longueur:N -> sORF du gène hôte X, de longueur N aa si N est donné.
    Si plusieurs candidats restent, aucun niveau manuel n'est donné (prudence)."""
    for m in manuel:
        regles = dict(x.split(":", 1) for x in m["appariement"].split(";") if ":" in x)
        cand = orfs
        if "symbole" in regles:
            s = set(regles["symbole"].split("|"))
            cand = [o for o in cand if o["gene_hote"] in s]
            # On préfère la CDS GENCODE si elle existe.
            gc = [o for o in cand if o["source"] == "gencode_cds_canonique"]
            cand = gc or cand
        if "hote" in regles:
            cand = [o for o in cand if o["gene_hote"] == regles["hote"]]
        if "longueur" in regles:
            n = int(regles["longueur"])
            cand = [o for o in cand if o.get("longueur_aa") == n]
        if "type" in regles:
            cand = [o for o in cand if o.get("orf_type", regles["type"]) == regles["type"]]
        if len(cand) == 1:
            cand[0]["manuel_peptide"] = m["peptide"]
            cand[0]["manuel_niveau"] = m["niveau_mito"]
            log.append("manuel %s -> %s" % (m["peptide"], cand[0]["orf_id"]))
        else:
            log.append("manuel %s : %d candidat(s), non apparié%s" % (
                m["peptide"], len(cand),
                (" (" + ", ".join(o["orf_id"] for o in cand[:5]) + ")") if cand else ""))
            for o in cand:
                o.setdefault("justif_extra", []).append(
                    "candidat ambigu pour %s (table manuelle)" % m["peptide"])


def dedoublonner(orfs):
    vus = {}
    for o in orfs:
        k = (o["chrom"], o["strand"], tuple(o["exons"]))
        if k in vus:
            p = vus[k]
            if o.get("etude") and o["etude"] not in p["etude"].split(";"):
                p["etude"] = ";".join(x for x in (p["etude"], o["etude"]) if x)
            if o["source"] == "gencode_cds_canonique":
                p["source"], p["orf_id"] = o["source"], o["orf_id"]
        else:
            vus[k] = o
    return list(vus.values())


# ---------------------------------------------------------------------------
# GTF et BED pour VEP
# ---------------------------------------------------------------------------

def segments_tx(ex, strand, i, j):
    """Intervalles génomiques de la portion [i, j) du transcrit (ordre 5'->3')."""
    seq = ex if strand == "+" else list(reversed(ex))
    out, pos = [], 0
    for a, b in seq:
        n = b - a
        lo, hi = max(i, pos), min(j, pos + n)
        if lo < hi:
            if strand == "+":
                out.append((a + lo - pos, a + hi - pos))
            else:
                out.append((b - (hi - pos), b - (lo - pos)))
        pos += n
    return sorted(out)


def ecrire_gtf_bed(orfs, gtf, bed, prefixe):
    lg = []
    lb = []
    for o in orfs:
        ex, s = o["exons"], o["strand"]
        n = sum(b - a for a, b in ex)
        a_stop = o.get("cds_nt", "")[-3:].upper() in ("TAA", "TAG", "TGA") or o["peptide"].endswith("*")
        fin_cds = n - 3 if a_stop else n
        chrom = prefixe + o["chrom"]
        gid = re.sub(r"[^A-Za-z0-9_.:|-]", "_", o["orf_id"])
        attr = ('gene_id "%s"; transcript_id "%s.t"; gene_name "%s"; gene_biotype "protein_coding"; '
                'transcript_biotype "protein_coding"; tag "chemin3_sorf"; niveau_mito "%s";'
                % (gid, gid, o["gene_hote"] or gid, o["niveau_mito"]))
        lo, hi = ex[0][0], ex[-1][1]
        lg.append((chrom, lo, ["gene", lo + 1, hi, ".", s, ".", attr]))
        lg.append((chrom, lo, ["transcript", lo + 1, hi, ".", s, ".", attr]))
        for k, (a, b) in enumerate(ex):
            lg.append((chrom, a, ["exon", a + 1, b, ".", s, ".", attr + ' exon_number "%d";' % (k + 1)]))
        # Phase GTF : nombre de nt à sauter pour atteindre le prochain codon.
        cumul = 0
        segs = segments_tx(ex, s, 0, fin_cds)
        for a, b in (segs if s == "+" else list(reversed(segs))):
            phase = (3 - cumul % 3) % 3
            lg.append((chrom, a, ["CDS", a + 1, b, ".", s, str(phase), attr]))
            cumul += b - a
        for a, b in segments_tx(ex, s, 0, 3):
            lg.append((chrom, a, ["start_codon", a + 1, b, ".", s, "0", attr]))
        if a_stop:
            for a, b in segments_tx(ex, s, n - 3, n):
                lg.append((chrom, a, ["stop_codon", a + 1, b, ".", s, "0", attr]))
        for a, b in ex:
            lb.append((chrom, a, b, "%s|%s|%s" % (gid, o["niveau_mito"], o["gene_hote"]), s))
    lg.sort(key=lambda x: (x[0], x[1]))
    with open(gtf, "w") as f:
        for chrom, _, c in lg:
            f.write("\t".join([chrom, "chemin3"] + [str(x) for x in c]) + "\n")
    lb.sort()
    with open(bed, "w") as f:
        for chrom, a, b, n, s in lb:
            f.write("%s\t%d\t%d\t%s\t0\t%s\n" % (chrom, a, b, n, s))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phaseI-bed")
    ap.add_argument("--phaseI-fa")
    ap.add_argument("--bed12", action="append", default=[], help="BED12 de sORF (GENCODE Ribo-seq ORFs, OpenProt)")
    ap.add_argument("--bed12-fa", action="append", default=[], help="FASTA protéique associé (même ordre)")
    ap.add_argument("--bed12-source", action="append", default=[])
    ap.add_argument("--bed12-garder", action="append", default=[],
                    help="regex sur le nom (même ordre que --bed12 ; « - » = tout) : OpenProt '^IP_' = altprots seuls")
    ap.add_argument("--genome", help="FASTA hg38 (indexé .fai de préférence)")
    ap.add_argument("--gencode-gtf", help="GTF GENCODE pour les micropeptides déjà annotés codants")
    ap.add_argument("--manuel", default=L.ici("..", "micropeptides_mito_etablis.tsv"))
    ap.add_argument("--mitocarta", help="MitoCarta3.0 humaine exportée en TSV")
    ap.add_argument("--proteome-mito", action="append", default=[],
                    help="peptides ou orf_id détectés dans un protéome mitochondrial publié")
    ap.add_argument("--targetp")
    ap.add_argument("--mitofates")
    ap.add_argument("--tm", action="append", default=[], help="sortie 3 lignes TMbed ou DeepTMHMM")
    ap.add_argument("--fasta-peptides", help="écrit les peptides (pour lancer TargetP, MitoFates, TMbed)")
    ap.add_argument("--prefixe-chr", default="chr", help="préfixe des chromosomes dans le GTF/BED")
    ap.add_argument("--sortie", required=True, help="préfixe des sorties")
    a = ap.parse_args()

    orfs = []
    if a.phaseI_bed:
        orfs += L.lire_phaseI(a.phaseI_bed, a.phaseI_fa)
    for i, b in enumerate(a.bed12):
        fa = a.bed12_fa[i] if i < len(a.bed12_fa) else None
        src = a.bed12_source[i] if i < len(a.bed12_source) else "bed12"
        g = a.bed12_garder[i] if i < len(a.bed12_garder) and a.bed12_garder[i] != "-" else None
        orfs += L.lire_bed12(b, fa, src, g)
    manuel = lire_manuel(a.manuel)
    if a.gencode_gtf:
        symb = set()
        for m in manuel:
            for x in m["appariement"].split(";"):
                if x.startswith("symbole:"):
                    symb |= set(x[8:].split("|"))
        orfs += gtf_cds_canoniques(a.gencode_gtf, symb)

    log = []
    n0 = len(orfs)
    orfs = [o for o in orfs if o["exons"] and o["chrom"] not in ("M", "MT")]
    log.append("sORF lus : %d ; après exclusion chrM et vides : %d" % (n0, len(orfs)))
    orfs = dedoublonner(orfs)
    log.append("après dédoublonnage par coordonnées : %d" % len(orfs))
    if a.gencode_gtf:
        n = hote_par_chevauchement(orfs, a.gencode_gtf)
        log.append("gène hôte attribué par chevauchement avec le GTF (même brin, sinon brin opposé) : %d" % n)

    genome = L.Genome(a.genome) if a.genome else None
    discord = 0
    for o in orfs:
        if genome:
            nt = L.cds_nt(o, genome)
            o["cds_nt"], o["codon_init"] = nt, nt[:3]
            tr = L.traduire(nt)
            if not o["peptide"]:
                o["peptide"] = tr
            elif tr.rstrip("*")[1:] != o["peptide"].rstrip("*")[1:]:
                discord += 1  # le 1er acide aminé peut différer (départ non AUG lu M)
        o.update(L.caracteristiques(o["peptide"]))
    if genome:
        log.append("peptides discordants avec la traduction du génome : %d" % discord)

    apparier_manuel(orfs, manuel, log)

    tp = lire_targetp(a.targetp) if a.targetp else {}
    mf = lire_mitofates(a.mitofates) if a.mitofates else {}
    tm = {}
    for t in a.tm:
        tm.update(lire_tm_3lignes(t))
    prot = set()
    for p in a.proteome_mito:
        prot |= lire_ensemble(p, ("peptide", "sequence", "orf_id"))
    mc = lire_mitocarta(a.mitocarta) if a.mitocarta else set()
    # Clé stable d'un lancement à l'autre (empreinte du peptide) pour relier les
    # sorties des outils de prédiction lancés sur --fasta-peptides.
    cle = {o["orf_id"]: "p" + hashlib.md5(o["peptide"].rstrip("*").encode()).hexdigest()[:12]
           for o in orfs}
    for o in orfs:
        k = cle[o["orf_id"]]
        o["targetp"] = tp.get(k, "absent") if tp else "non_calcule"
        o["mitofates"] = mf.get(k, "absent") if mf else "non_calcule"
        o["tm_outil"] = tm.get(k, "absent") if tm else "non_calcule"
        o["proteome_mito"] = ("oui" if (o["orf_id"] in prot or o["peptide"].rstrip("*") in prot) else "non") if prot else "non_calcule"
        o["hote_mitocarta"] = ("oui" if o["gene_hote"] in mc else "non") if mc else "non_calcule"
        niv, j = L.niveau_mito(o)
        j += o.get("justif_extra", [])
        if o["hote_mitocarta"] == "oui" and o["source"] != "gencode_cds_canonique":
            j.append("gène hôte dans MitoCarta3.0 (un variant peut aussi agir sur la CDS principale : voie régulatrice)")
        o["niveau_mito"], o["justification_mito"] = niv, " ; ".join(j)

    if a.fasta_peptides:
        with open(a.fasta_peptides, "w") as f:
            ecrits = set()
            for o in orfs:
                p, k = o["peptide"].rstrip("*"), cle[o["orf_id"]]
                if len(p) >= 10 and k not in ecrits:
                    ecrits.add(k)
                    f.write(">%s\n%s\n" % (k, p))
        with open(a.fasta_peptides + ".ids.tsv", "w") as f:
            f.write("cle\torf_id\n")
            for k, v in cle.items():
                f.write("%s\t%s\n" % (v, k))

    L.ecrire_table(orfs, a.sortie + ".table_sorf.tsv.gz")
    ecrire_gtf_bed(orfs, a.sortie + ".gtf", a.sortie + ".bed", a.prefixe_chr)
    cpt = defaultdict(int)
    for o in orfs:
        cpt[o["niveau_mito"]] += 1
    log.append("niveaux mito : " + ", ".join("%s=%d" % (n, cpt[n]) for n in reversed(L.ORDRE_NIVEAU)))
    with open(a.sortie + ".log.txt", "w") as f:
        f.write("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
