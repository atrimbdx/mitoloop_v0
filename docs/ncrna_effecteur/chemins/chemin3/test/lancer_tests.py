#!/usr/bin/env python3
"""Test du chemin 3 sur un JEU CONSTRUIT À LA MAIN (synthétique).

ATTENTION : le génome « chrT1 », les sORF TOYA à TOYD, les prédictions et les
variants ci-dessous sont inventés pour vérifier le code. Ce ne sont ni des données
réelles ni des résultats. Aucun variant de patient n'est utilisé.

Lancer : python3 lancer_tests.py   (écrit dans jeu_main/, puis vérifie les sorties)
"""
import hashlib
import os
import random
import subprocess
import sys

ICI = os.path.dirname(os.path.abspath(__file__))
J = os.path.join(ICI, "jeu_main")
S = os.path.join(ICI, "..", "scripts")
sys.path.insert(0, S)
import sorf_lib as L  # noqa: E402

os.makedirs(J, exist_ok=True)
CODON = {}
for k, v in L.CODE.items():
    CODON.setdefault(v, k)


def coder(p, stop="TAA"):
    return "".join(CODON[a] for a in p) + stop


random.seed(1)
g = [random.choice("ACGT") for _ in range(3000)]


def poser(a, s):
    g[a:a + len(s)] = list(s)


# TOYA : brin +, 1 exon, peptide de 45 aa (table manuelle synthétique -> fort)
pA = "MRRLLRKALSRLAGQPTEVNSWYHCFIKDAMPGQSTNVWYHCFIK"
cA = coder(pA)
A0 = 100
poser(A0, cA)
# TOYB : brin -, 2 exons (60 + 63 nt), hélice hydrophobe -> hypothese
pB = "MSTN" + "L" * 22 + "KQDEGSTNPQRSTV"
cB = coder(pB, "TGA")
rc = L.revcomp(cB)
poser(1000, rc[:60])
poser(1200, rc[60:])
# TOYC : brin +, peptide dans une liste de protéome mito synthétique -> modere
pC = "MSDEKQNPTSGAWHYVCFIRLMSDEKQNPT"
cC = coder(pC, "TAG")
poser(2000, cC)
# TOYD : brin +, TargetP synthétique = mTP -> faible
pD = "MSDEKQNPTSGAWHYVCFIRLMSDEKQNPTGAWHY"
cD = coder(pD)
poser(2500, cD)

with open(os.path.join(J, "genome_synthetique.fa"), "w") as f:
    f.write(">chrT1 JEU_CONSTRUIT_A_LA_MAIN\n" + "".join(g) + "\n")


def bed12(nom, ex, st, hote):
    a, b = ex[0][0], ex[-1][1]
    return "\t".join(map(str, ["chrT1", a, b, nom, 0, st, a, b, 0, len(ex),
                               ",".join(str(y - x) for x, y in ex),
                               ",".join(str(x - a) for x, y in ex), hote])) + "\n"


with open(os.path.join(J, "sorf_synthetiques.bed"), "w") as f:
    f.write("# JEU CONSTRUIT A LA MAIN (synthetique)\n")
    f.write(bed12("TOYA_orf", [(A0, A0 + len(cA))], "+", "TOYA"))
    f.write(bed12("TOYB_orf", [(1000, 1060), (1200, 1263)], "-", "TOYB"))
    f.write(bed12("TOYC_orf", [(2000, 2000 + len(cC))], "+", "TOYC"))
    f.write(bed12("TOYD_orf", [(2500, 2500 + len(cD))], "+", "TOYD"))
with open(os.path.join(J, "manuel_synthetique.tsv"), "w") as f:
    f.write("peptide\tappariement\tniveau_mito\n")
    f.write("TOYA-peptide (synthetique)\thote:TOYA;longueur:45\tfort\n")
with open(os.path.join(J, "proteome_mito_synthetique.txt"), "w") as f:
    f.write("peptide\n" + pC + "\n")
kD = "p" + hashlib.md5(pD.encode()).hexdigest()[:12]
with open(os.path.join(J, "targetp_synthetique.targetp2"), "w") as f:
    f.write("# JEU CONSTRUIT A LA MAIN\n# ID\tPrediction\tnoTP\tSP\tmTP\n%s\tmTP\t0.1\t0.0\t0.9\n" % kD)

kB = "p" + hashlib.md5(pB.encode()).hexdigest()[:12]
kC = "p" + hashlib.md5(pC.encode()).hexdigest()[:12]
with open(os.path.join(J, "tm_synthetique.3line"), "w") as f:
    f.write(">%s\n%s\n%s\n" % (kB, pB, "...." + "H" * 22 + "." * 14))
    f.write(">%s\n%s\n%s\n" % (kD, pD, "." * len(pD)))
with open(os.path.join(J, "mitofates_synthetique.tsv"), "w") as f:
    f.write("Sequence ID\tProbability of presequence\tPrediction\n")
    f.write("%s\t0.9\tPossessing mitochondrial presequence\n%s\t0.1\tNo mitochondrial presequence\n" % (kC, kB))

# Variants (pos 1-based) et conséquence / niveau attendus
V = []


def v(pos, ref, alt, csq, niv, note):
    V.append((pos, ref, alt, csq, niv, note))


v(A0 + 1, "A", "G", "start_lost", "fort", "TOYA ATG>GTG")
s5 = A0 + 12  # codon 5 (L, TTA) de TOYA, 1re base
v(s5 + 1, cA[12], "T" if cA[12] != "T" else "G", None, None, "TOYA codon 5")
v(A0 + 18 + 1, cA[18], "T", None, None, "TOYA codon 7 AAA>TAA")
v(A0 + 3 * 44 + 1, cA[132], "T", None, None, "TOYA codon 45")
v(A0 + 3 * 45 + 3, "A", "C", "stop_lost", "modere", "TOYA TAA>TAC")
v(A0 + 30, cA[29] + cA[30], cA[29], "frameshift", "fort", "TOYA délétion 1 nt")
v(A0 + 21, "N", "A", None, "aucun", "référence discordante")
v(1030, rc[29], "A" if rc[29] != "A" else "C", None, None, "TOYB exon 1 (fin du sORF, brin -)")
v(1061, g[1060], "A" if g[1060] != "A" else "C", "splice_site_orf", "hypothese", "TOYB site d'épissage")
v(1250, g[1249], g[1249] + "G", "frameshift", "hypothese", "TOYB insertion 1 nt")
v(2000 + 3 * 5 + 3, cC[17], "A" if cC[17] != "A" else "C", None, None, "TOYC codon 6")
v(2500 + 2, "T", "C", "start_lost", "faible", "TOYD ATG>ACG")
v(500, g[499], "A" if g[499] != "A" else "C", None, "aucun", "intergénique")
with open(os.path.join(J, "variants_synthetiques.tsv"), "w") as f:
    f.write("chr\tpos\tref\talt\tnote\n")
    for pos, ref, alt, _, _, note in V:
        f.write("chrT1\t%d\t%s\t%s\t%s\n" % (pos, ref, alt, note))
    f.write("chrM\t3243\tA\tG\tchrM hors périmètre\n")

py = sys.executable
subprocess.run([py, os.path.join(S, "preparer_sorf.py"), "--bed12", os.path.join(J, "sorf_synthetiques.bed"),
                "--bed12-source", "jeu_main", "--genome", os.path.join(J, "genome_synthetique.fa"),
                "--manuel", os.path.join(J, "manuel_synthetique.tsv"),
                "--proteome-mito", os.path.join(J, "proteome_mito_synthetique.txt"),
                "--targetp", os.path.join(J, "targetp_synthetique.targetp2"),
                "--mitofates", os.path.join(J, "mitofates_synthetique.tsv"),
                "--tm", os.path.join(J, "tm_synthetique.3line"),
                "--fasta-peptides", os.path.join(J, "peptides.fa"),
                "--sortie", os.path.join(J, "sortie")], check=True)
subprocess.run([py, os.path.join(S, "annoter_chemin3.py"), "--variants", os.path.join(J, "variants_synthetiques.tsv"),
                "--table", os.path.join(J, "sortie.table_sorf.tsv.gz"),
                "--sortie", os.path.join(J, "variants_annotes.tsv")], check=True)

tab = {o["orf_id"]: o for o in L.lire_table(os.path.join(J, "sortie.table_sorf.tsv.gz"))}
erreurs = []
for k, n in (("TOYA_orf", "fort"), ("TOYB_orf", "hypothese"), ("TOYC_orf", "modere"), ("TOYD_orf", "faible")):
    if tab[k]["niveau_mito"] != n:
        erreurs.append("%s : niveau %s, attendu %s" % (k, tab[k]["niveau_mito"], n))
for k, p in (("TOYA_orf", pA), ("TOYB_orf", pB), ("TOYC_orf", pC), ("TOYD_orf", pD)):
    if L.traduire(tab[k]["cds_nt"]).rstrip("*") != p:
        erreurs.append("%s : traduction incorrecte" % k)

lignes = open(os.path.join(J, "variants_annotes.tsv")).read().strip().split("\n")[1:]
for (pos, ref, alt, csq, niv, note), l in zip(V, lignes):
    res = l.split("\t")[-1]
    score, niveau, just = res.split("|", 2)
    print("%-38s %-6s %-10s %s" % (note, score, niveau, just[:110]))
    if csq and not just.startswith(csq):
        erreurs.append("%s : conséquence %s attendue, obtenu %s" % (note, csq, just[:40]))
    if niv and niveau != niv:
        erreurs.append("%s : niveau %s attendu, obtenu %s" % (note, niv, niveau))
attendus_texte = {
    "TOYA codon 5": ("missense", "L5V", "modere"),
    "TOYA codon 7 AAA>TAA": ("stop_gained", "", "fort"),
    "TOYA codon 45": ("stop_gained", "", "modere"),
    "TOYB exon 1 (fin du sORF, brin -)": ("missense", "S32A", "hypothese"),
    "TOYC codon 6": ("missense", "Q6H", "faible"),
}
for (pos, ref, alt, csq, niv, note), l in zip(V, lignes):
    if note in attendus_texte:
        c, motif, n = attendus_texte[note]
        res = l.split("\t")[-1]
        if not res.split("|")[2].startswith(c) or motif not in res or res.split("|")[1] != n:
            erreurs.append("%s : attendu %s %s %s, obtenu %s" % (note, c, motif, n, res[:80]))
if tab["TOYB_orf"]["tm_outil"] != "oui" or "outil dédié" not in tab["TOYB_orf"]["justification_mito"]:
    erreurs.append("lecture TM 3 lignes incorrecte")
if tab["TOYC_orf"]["mitofates"] != "oui" or tab["TOYB_orf"]["mitofates"] != "non":
    erreurs.append("lecture MitoFates incorrecte")
if "chrM hors périmètre" not in lignes[-1]:
    erreurs.append("ligne chrM absente")
elif "|aucun|" not in lignes[-1]:
    erreurs.append("chrM non exclu")

# Contrôle du GTF : la CDS exclut le stop et la somme des CDS vaut 3 x longueur.
cds_len = {}
for l in open(os.path.join(J, "sortie.gtf")):
    c = l.split("\t")
    if c[2] == "CDS":
        gid = c[8].split('"')[1]
        cds_len[gid] = cds_len.get(gid, 0) + int(c[4]) - int(c[3]) + 1
for k, p in (("TOYA_orf", pA), ("TOYB_orf", pB)):
    if cds_len.get(k) != 3 * len(p):
        erreurs.append("GTF %s : CDS %s nt, attendu %d" % (k, cds_len.get(k), 3 * len(p)))

print("\nJEU CONSTRUIT À LA MAIN : %s" % ("OK" if not erreurs else "ÉCHECS"))
for e in erreurs:
    print("  - " + e)
sys.exit(1 if erreurs else 0)
