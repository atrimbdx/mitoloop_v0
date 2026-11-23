"""Fonctions communes du chemin 3 (micropeptides) : lecture des sORF, caractéristiques
des peptides, règles de niveau mito et de conséquence, toutes fixées a priori.

Aucune dépendance hors bibliothèque standard (la VM est hors ligne).
Coordonnées internes : 0-based, demi-ouvertes [début, fin), comme le BED.
"""
import re
import csv
import gzip
import os

# ---------------------------------------------------------------------------
# Constantes a priori (ne jamais ajuster sur des variants de patients ni sur des
# sorties AlphaGenome / SpliceAI / SpliceAI2 / AlphaFold 3)
# ---------------------------------------------------------------------------

# Hydropathie de Kyte et Doolittle (1982, doi:10.1016/0022-2836(82)90515-0).
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
# Fenêtre de 19 résidus et moyenne >= 1,6 : critère proposé par Kyte et Doolittle
# pour repérer un segment transmembranaire. Repli si aucun outil dédié n'a tourné.
KD_FENETRE = 19
KD_SEUIL_TM = 1.6

# Indice heuristique de présequence (repli, niveau « hypothese » seulement) :
# les présequences sont riches en résidus basiques et pauvres en résidus acides.
# Seuils choisis a priori, non calibrés : >= 3 R/K et <= 1 D/E dans les 20 premiers
# résidus, peptide d'au moins 30 aa (sinon la présequence ne laisserait presque rien).
MTS_NTERM = 20
MTS_MIN_BASIQUES = 3
MTS_MAX_ACIDES = 1
MTS_LONGUEUR_MIN = 30

# Poids des niveaux mito du peptide.
POIDS_MITO = {"fort": 1.0, "modere": 0.7, "faible": 0.4, "hypothese": 0.15, "aucun": 0.0}
ORDRE_NIVEAU = ["aucun", "hypothese", "faible", "modere", "fort"]

# Poids des conséquences sur le sORF.
POIDS_CSQ = {
    "start_lost": 1.0,
    "stop_gained": 1.0,
    "frameshift": 1.0,
    "splice_site_orf": 0.8,
    "stop_lost": 0.6,
    "inframe_indel": 0.6,
    "missense": 0.5,
    "stop_retained": 0.0,
    "synonymous": 0.0,
}
# Troncature retirant < 10 % du peptide : poids réduit, par analogie avec l'arbre
# PVS1 (Abou Tayoun et al. 2018, doi:10.1002/humu.23626).
POIDS_TRONCATURE_DISTALE = 0.5
FRACTION_DISTALE = 0.10

# Seuils du niveau final à partir du score = poids_csq x poids_mito.
SEUILS_SCORE = [(0.8, "fort"), (0.45, "modere"), (0.2, "faible"), (1e-9, "hypothese")]

CODE = {}
_b = "TCAG"
_aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_i = 0
for x in _b:
    for y in _b:
        for z in _b:
            CODE[x + y + z] = _aa[_i]
            _i += 1

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def traduire(nt):
    return "".join(CODE.get(nt[i:i + 3].upper(), "X") for i in range(0, len(nt) - 2, 3))


def ouvrir(chemin, mode="rt"):
    if chemin.endswith(".gz"):
        return gzip.open(chemin, mode)
    return open(chemin, mode)


def norm_chrom(c):
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


# ---------------------------------------------------------------------------
# Lecture des sORF
# ---------------------------------------------------------------------------

def lire_fasta(chemin):
    seqs, cle = {}, None
    with ouvrir(chemin) as f:
        for l in f:
            l = l.strip()
            if not l:
                continue
            if l.startswith(">"):
                cle = l[1:].split()[0].split("|")[0]   # OpenProt : >IP_123|... ; GENCODE : nom--étude
                seqs[cle] = []
            elif cle is not None:
                seqs[cle].append(l)
    return {k: "".join(v) for k, v in seqs.items()}


def lire_phaseI(bed, fasta, source="riboseq_phaseI_entree"):
    """Format du dépôt jorruior/gencode-riboseqORFs (phaseI/) : BED6 par exon,
    coordonnées 1-based inclusives (stop compris), nom en colonne 4, étude en 5 ;
    FASTA protéique avec en-tête nom--étude."""
    pep = lire_fasta(fasta)
    orfs = {}
    with ouvrir(bed) as f:
        for l in f:
            if not l.strip() or l.startswith(("#", "track")):
                continue
            c = l.rstrip("\n").split("\t")
            oid = c[3] + "--" + c[4]
            o = orfs.setdefault(oid, {"orf_id": oid, "source": source, "chrom": norm_chrom(c[0]),
                                      "strand": c[5], "exons": [], "etude": c[4]})
            o["exons"].append((int(c[1]) - 1, int(c[2])))
    for o in orfs.values():
        o["exons"].sort()
        o["peptide"] = pep.get(o["orf_id"], "")
        o["gene_hote"] = gene_depuis_nom(o["orf_id"])
    return list(orfs.values())


def gene_depuis_nom(nom):
    """Phase II du catalogue GENCODE (c1riboseqorf2598:LINC00467:lncRNA) : gène hôte en 2e champ ;
    phase I : le nom commence par le symbole hôte (GENE_..., GENE:...). Accessions OpenProt (IP_, II_)
    sans gène : hôte attribué ensuite par chevauchement avec le GTF (preparer_sorf --gencode-gtf)."""
    if re.match(r"^(IP|II)_\d+", nom):
        return ""
    m = re.match(r"^c[0-9XYM]+riboseqorf\d+:([^:]+):", nom)
    if m:
        return m.group(1)
    n = nom.split("--")[0]
    for sep in (":", "_"):
        if sep in n:
            return n.split(sep)[0]
    return ""


def lire_bed12(bed, fasta=None, source="bed12", garder=None):
    """BED12 standard (0-based) ; thickStart/thickEnd = CDS du sORF (stop compris
    si présent). Le nom (col. 4) sert d'identifiant ; si la colonne 13 existe, elle
    est prise comme gène hôte."""
    pep = lire_fasta(fasta) if fasta else {}
    orfs = []
    with ouvrir(bed) as f:
        for l in f:
            if not l.strip() or l.startswith(("#", "track", "browser")):
                continue
            c = l.rstrip("\n").split("\t")
            st, n = int(c[1]), c[3]
            if garder and not re.search(garder, n):
                continue
            ts, te = int(c[6]), int(c[7])
            tailles = [int(x) for x in c[10].strip(",").split(",")]
            debuts = [int(x) for x in c[11].strip(",").split(",")]
            ex = []
            for t, d in zip(tailles, debuts):
                a, b = max(st + d, ts), min(st + d + t, te)
                if a < b:
                    ex.append((a, b))
            orfs.append({"orf_id": n, "source": source, "chrom": norm_chrom(c[0]), "strand": c[5],
                         "exons": sorted(ex), "peptide": pep.get(n, ""),
                         "gene_hote": c[12] if len(c) > 12 else gene_depuis_nom(n), "etude": ""})
    return orfs


# ---------------------------------------------------------------------------
# Séquences
# ---------------------------------------------------------------------------

class Genome:
    """Accès FASTA : pyfaidx si présent (gros génome sur Curta), sinon chargement
    complet (petits fichiers de test)."""

    def __init__(self, chemin):
        self.fa = None
        try:
            import pyfaidx  # noqa
            self.fa = pyfaidx.Fasta(chemin, as_raw=True, sequence_always_upper=True)
            self.noms = {norm_chrom(k): k for k in self.fa.keys()}
        except ImportError:
            s = lire_fasta(chemin)
            self.seqs = {norm_chrom(k): v.upper() for k, v in s.items()}

    def get(self, chrom, a, b):
        chrom = norm_chrom(chrom)
        if self.fa is not None:
            return self.fa[self.noms[chrom]][a:b]
        return self.seqs[chrom][a:b]


def cds_nt(orf, genome):
    s = "".join(genome.get(orf["chrom"], a, b) for a, b in orf["exons"])
    return revcomp(s) if orf["strand"] == "-" else s


# ---------------------------------------------------------------------------
# Caractéristiques du peptide (a priori)
# ---------------------------------------------------------------------------

def kd_max(p, w=KD_FENETRE):
    p = p.rstrip("*")
    if len(p) < w:
        return None
    vals = [KD.get(a, 0.0) for a in p]
    s = sum(vals[:w])
    m = s
    for i in range(w, len(vals)):
        s += vals[i] - vals[i - w]
        m = max(m, s)
    return round(m / w, 3)


def indices_mts(p):
    p = p.rstrip("*")
    n = p[:MTS_NTERM]
    basiques = sum(a in "RK" for a in n)
    acides = sum(a in "DE" for a in n)
    ok = len(p) >= MTS_LONGUEUR_MIN and basiques >= MTS_MIN_BASIQUES and acides <= MTS_MAX_ACIDES
    return basiques, acides, ok


def caracteristiques(p):
    km = kd_max(p)
    b, a, mts = indices_mts(p)
    return {"longueur_aa": len(p.rstrip("*")), "kd_max19": "" if km is None else km,
            "tm_kd": "oui" if (km is not None and km >= KD_SEUIL_TM) else "non",
            "nterm20_RK": b, "nterm20_DE": a, "mts_heuristique": "oui" if mts else "non"}


def niveau_mito(o):
    """Règle a priori du niveau mito d'un sORF (voir README, section critères).
    o : dict avec manuel_niveau, proteome_mito, targetp, mitofates, tm_outil,
    tm_kd, mts_heuristique, hote_mitocarta."""
    j = []
    if o.get("manuel_niveau") in ("fort", "modere"):
        j.append("table manuelle : %s (%s)" % (o.get("manuel_peptide", ""), o["manuel_niveau"]))
        return o["manuel_niveau"], j
    if o.get("proteome_mito") == "oui":
        j.append("peptide détecté dans un protéome mitochondrial publié")
        return "modere", j
    outil_mts = o.get("targetp") == "mTP" or o.get("mitofates") == "oui"
    if outil_mts:
        j.append("présequence prédite (%s)" % ",".join(
            x for x, ok in (("TargetP", o.get("targetp") == "mTP"),
                            ("MitoFates", o.get("mitofates") == "oui")) if ok))
        return "faible", j
    # Les replis (heuristique N-terminale, Kyte-Doolittle) ne servent que si l'outil
    # dédié n'a pas tourné : un outil dédié l'emporte toujours sur le repli.
    sans_outil_mts = o.get("targetp", "non_calcule") == "non_calcule" and \
        o.get("mitofates", "non_calcule") == "non_calcule"
    sans_outil_tm = o.get("tm_outil", "non_calcule") == "non_calcule"
    if sans_outil_mts and o.get("mts_heuristique") == "oui":
        j.append("indice N-terminal basique (repli, sans TargetP ni MitoFates)")
    if o.get("tm_outil") == "oui":
        j.append("hélice transmembranaire prédite (outil dédié), non spécifique de la mitochondrie")
    elif sans_outil_tm and o.get("tm_kd") == "oui":
        j.append("hélice transmembranaire prédite (repli Kyte-Doolittle), non spécifique de la mitochondrie")
    if j:
        return "hypothese", j
    return "aucun", ["aucun indice mitochondrial"]


def niveau_depuis_score(s):
    for seuil, n in SEUILS_SCORE:
        if s >= seuil:
            return n
    return "aucun"


# ---------------------------------------------------------------------------
# Table sORF (TSV) partagée entre Curta et la VM
# ---------------------------------------------------------------------------

COLONNES_TABLE = ["orf_id", "source", "etude", "gene_hote", "chrom", "strand", "exons",
                  "longueur_aa", "peptide", "cds_nt", "codon_init", "kd_max19", "tm_kd",
                  "nterm20_RK", "nterm20_DE", "mts_heuristique", "targetp", "mitofates",
                  "tm_outil", "proteome_mito", "hote_mitocarta", "manuel_peptide",
                  "manuel_niveau", "niveau_mito", "justification_mito"]


def exons_txt(ex):
    return ",".join("%d-%d" % e for e in ex)


def exons_lire(t):
    return [tuple(int(x) for x in e.split("-")) for e in t.split(",") if e]


def ecrire_table(orfs, chemin):
    with ouvrir(chemin, "wt") as f:
        w = csv.DictWriter(f, COLONNES_TABLE, delimiter="\t", extrasaction="ignore",
                           lineterminator="\n")
        w.writeheader()
        for o in orfs:
            r = dict(o)
            r["exons"] = exons_txt(o["exons"])
            w.writerow(r)


def lire_table(chemin):
    out = []
    with ouvrir(chemin) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            r["exons"] = exons_lire(r["exons"])
            out.append(r)
    return out


def ici(*p):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), *p)
