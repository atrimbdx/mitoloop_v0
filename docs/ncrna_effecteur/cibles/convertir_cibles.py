#!/usr/bin/env python3
"""Construit cibles_ncrna.tsv (cibles des éléments non codants effecteurs) pour reannot.

Bibliothèque standard uniquement. Chaque sous-commande lit une source et écrit des lignes au
format d'interface (10 colonnes, voir COLONNES). La sous-commande `fusion` concatène, complète
les ENSG manquants (HGNC) et dédoublonne.

  python3 convertir_cibles.py snodb       snoDB_All_V2.0.tsv                     > t_snodb.tsv
  python3 convertir_cibles.py curees      cibles_curees.tsv                      > t_curees.tsv
  python3 convertir_cibles.py pol3        partenaires_pol3_arnt_arnr.tsv         > t_pol3.tsv
  python3 convertir_cibles.py pseudogenes --hgnc hgnc_complete_set.txt --gtf gencode.gtf > t_pseudo.tsv
  python3 convertir_cibles.py humanin     --hgnc hgnc_complete_set.txt           > t_humanin.tsv
  python3 convertir_cibles.py lncrna_cis  --gtf gencode.gtf [--distance 50000]   > t_lnc.tsv
  python3 convertir_cibles.py mirtarbase  hsa_MTI.csv --hairpin hairpin.fa --hgnc hgnc_complete_set.txt > t_mtb.tsv
  python3 convertir_cibles.py generique   FICHIER --type-id X --col-element A --col-cible B [...] > t_gen.tsv
  python3 convertir_cibles.py fusion      --hgnc hgnc_complete_set.txt --out cibles_ncrna.tsv t_*.tsv
"""
import argparse
import csv
import gzip
import re
import sys

COLONNES = ["type_id", "element", "element_id", "cible", "type_cible", "nature",
            "position_cible", "preuve", "source", "reference", "region_guide"]
# region_guide (facultative, demandée par reannot) : région antisens du snoRNA au format guide:A-B
# (positions 1-based dans l'ARN mature du snoRNA) ; "." quand aucune source vérifiée ne la donne.
TYPES_CIBLE = {"arnm", "snrna", "arnr", "arnt", "proteine_partenaire", "gene_parent",
               "gene_voisin_cis", "arnt_mito_origine"}
NATURES = {"modification_2OMe", "pseudouridylation", "repression_graine", "clivage",
           "sequestration", "liaison", "cis", "origine"}
PREUVES = {"experimentale_humaine", "experimentale_animale", "prediction", "base_de_donnees"}

csv.field_size_limit(10**9)


def ouvrir(chemin):
    if chemin == "-":
        return sys.stdin
    if chemin.endswith(".gz"):
        return gzip.open(chemin, "rt", encoding="utf-8", errors="replace")
    return open(chemin, encoding="utf-8", errors="replace")


def ecrire(lignes, sortie=sys.stdout, entete=True):
    if entete:
        sortie.write("\t".join(COLONNES) + "\n")
    for l in lignes:
        valider(l)
        sortie.write("\t".join(str(l[c]) for c in COLONNES) + "\n")


def valider(l):
    for c in COLONNES:
        v = l.get(c, "")
        if v == "" or "\t" in str(v) or "\n" in str(v):
            raise ValueError(f"colonne {c} vide ou invalide : {l}")
    if l["type_cible"] not in TYPES_CIBLE:
        raise ValueError(f"type_cible inconnu : {l['type_cible']}")
    if l["nature"] not in NATURES:
        raise ValueError(f"nature inconnue : {l['nature']}")
    if l["preuve"] not in PREUVES:
        raise ValueError(f"preuve inconnue : {l['preuve']}")


def ligne(**kw):
    d = {c: kw.get(c, ".") for c in COLONNES}
    for c in COLONNES:
        if d[c] in (None, ""):
            d[c] = "."
    return d


# ---------------------------------------------------------------- HGNC
def lire_hgnc(chemin):
    """Renvoie (symbole -> ligne dict, ancien/alias -> symbole approuvé)."""
    par_symbole, alias = {}, {}
    with ouvrir(chemin) as f:
        for d in csv.DictReader(f, delimiter="\t"):
            s = d["symbol"]
            par_symbole[s] = d
            for champ in ("prev_symbol", "alias_symbol"):
                for a in d.get(champ, "").strip('"').split("|"):
                    if a and a not in alias:
                        alias[a] = s
    return par_symbole, alias


# ---------------------------------------------------------------- snoDB
ARNR = {"18S": "RNA18SN1", "28S": "RNA28SN1", "5.8S": "RNA5-8SN1", "5S": "RNA5S1",
        "5S_rRNA": "RNA5S1"}
SNRNA = {"U1": "RNU1-1", "U2": "RNU2-1", "U4": "RNU4-1", "U5": "RNU5A-1", "U6": "RNU6-1",
         "U7": "RNU7-1", "U11": "RNU11", "U12": "RNU12", "U4atac": "RNU4ATAC",
         "U6atac": "RNU6ATAC", "U3": "SNORD3A"}
RE_ARNR = re.compile(r"^(18S|28S|5\.8S|5S)(?:-(\d+))?$")
RE_SNRNA = re.compile(r"^(U4atac|U6atac|U11|U12|U1|U2|U3|U4|U5|U6|U7)(?:\.(\d+))?(?:-(\d+))?$")
U3_U8_U13 = {"SNORD3A", "SNORD3B-1", "SNORD3B-2", "SNORD3C", "SNORD3D", "SNORD118", "SNORD13",
             "U3", "U8", "U13"}


def classe_scarna(seq):
    """Type de guide d'un scaRNA déduit des boîtes (heuristique documentée) :
    ACA dans les 6 derniers nt -> H/ACA ; sinon boîte C (RUGAUGA, 30 premiers nt) ou
    boîte D (CUGA, 12 derniers nt) -> C/D ; sinon inconnu."""
    s = seq.upper().replace("U", "T")
    if "ACA" in s[-6:]:
        return "H/ACA"
    if re.search("[AG]TGATGA", s[:30]) or "CTGA" in s[-12:]:
        return "C/D"
    return "inconnu"


def type_id_snodb(nom, box):
    if re.match(r"^SNORD11[3-6]", nom):
        return "snorna_orphelin_15q11_14q32"
    if nom in U3_U8_U13 or re.match(r"^SNORD3[A-D]", nom):
        return "snorna_u3_u8_u13"
    if box == "scaRNA" or nom.startswith("SCARNA"):
        return "scarna"
    if box in ("C/D", "SNORD-like"):
        return "snorna_cd"
    if box in ("H/ACA", "AluACA"):
        return "snorna_haca"
    if box == "telomerase RNA":
        return "terc"
    if nom.startswith("SNORD"):
        return "snorna_cd"
    if nom.startswith("SNORA"):
        return "snorna_haca"
    return None


def nature_modif(guide):
    return {"C/D": "modification_2OMe", "H/ACA": "pseudouridylation"}.get(guide)


def cmd_snodb(args):
    ref = "snoDB v2.0 (Bergeron 2023, doi:10.1093/nar/gkac835)"
    stats = {"sans_classe": 0, "ignores_snorna_mirna": 0}
    sortie = []
    # Type de guide des scaRNA : un symbole peut avoir plusieurs entrées (copies) ; on retient
    # la première classe déterminée par les boîtes pour toutes ses copies.
    guide_scarna = {}
    with ouvrir(args.fichier) as f:
        for d in csv.DictReader(f, delimiter="\t"):
            if d["box_type"].strip() == "scaRNA":
                g = classe_scarna(d["sequence"])
                if g != "inconnu":
                    guide_scarna.setdefault(d["gene_name"].strip() or d["snodb_id"], g)
    with ouvrir(args.fichier) as f:
        for d in csv.DictReader(f, delimiter="\t"):
            nom = d["gene_name"].strip()
            box = d["box_type"].strip()
            tid = type_id_snodb(nom, box)
            if tid is None:
                stats["sans_classe"] += 1
                continue
            guide = box if box in ("C/D", "H/ACA") else (
                "C/D" if box == "SNORD-like" else "H/ACA" if box == "AluACA" else None)
            note_guide = ""
            if box == "scaRNA":
                guide = guide_scarna.get(nom or d["snodb_id"], classe_scarna(d["sequence"]))
                note_guide = f" ; type de guide du scaRNA déduit des boîtes : {guide}"
            ensg = [e for e in d["ensembl_id"].split(";") if e] or ["."]
            element = nom or (ensg[0] if ensg[0] != "." else d["snodb_id"])
            cibles = []
            for tok in filter(None, d["rrna_targets"].split(";")):
                m = RE_ARNR.match(tok)
                if m:
                    cibles.append((ARNR[m.group(1)], "arnr", m.group(2), tok))
                elif tok in ARNR:
                    cibles.append((ARNR[tok], "arnr", None, tok))
                else:
                    cibles.append((tok, "arnr", None, tok))
            for tok in filter(None, d["snrna_targets"].split(";")):
                m = RE_SNRNA.match(tok)
                if m:
                    cibles.append((SNRNA[m.group(1)], "snrna", m.group(3), tok))
                else:
                    cibles.append((tok, "snrna", None, tok))
            for col, tc in (("protein_coding_targets", "arnm"), ("lncrna_targets", "arnm"),
                            ("pseudogene_targets", "arnm"), ("other_targets", "arnm"),
                            ("ncrna_targets", "snrna"), ("trna_targets", "arnt")):
                for tok in filter(None, d[col].split(";")):
                    t = tc
                    if col == "ncrna_targets" and tok.startswith("MT-T"):
                        t = "arnt"
                    cibles.append((tok, t, None, tok + " [" + col + "]"))
            for col in ("snorna_targets", "mirna_targets"):
                stats["ignores_snorna_mirna"] += len([t for t in d[col].split(";") if t])
            for cible, tc, pos, brut in cibles:
                if pos is not None:
                    nat = nature_modif(guide)
                    if nat is None:
                        nat = "liaison"
                    src = f"snoDB v2.0 site de modification guidé ({brut}){note_guide}"
                else:
                    nat = "liaison"
                    src = f"snoDB v2.0 interaction ARN-ARN sans position, haut débit ({brut})"
                for e in ensg:
                    sortie.append(ligne(type_id=tid, element=element, element_id=e, cible=cible,
                                        type_cible=tc, nature=nat, position_cible=pos or ".",
                                        preuve="base_de_donnees", source=src, reference=ref))
    ecrire(sortie)
    print(f"[snodb] {len(sortie)} lignes ; {stats}", file=sys.stderr)


# ---------------------------------------------------------------- curées
def cmd_curees(args):
    sortie = []
    with ouvrir(args.fichier) as f:
        for d in csv.DictReader(f, delimiter="\t"):
            sortie.append(ligne(**{c: d.get(c, ".") for c in COLONNES}))
    ecrire(sortie)
    print(f"[curees] {len(sortie)} lignes", file=sys.stderr)


# ---------------------------------------------------------------- Pol III (chemin 1)
POL3 = {"RMRP": ("rmrp_rpph1", ["RMRP"]), "RPPH1": ("rmrp_rpph1", ["RPPH1"]),
        "TERC": ("terc", ["TERC"]), "RN7SK": ("rn7sk", ["RN7SK"]),
        "RN7SL": ("srp_7sl", ["RN7SL1", "RN7SL2", "RN7SL3"]),
        "ARNt": ("trna_nucleaire", ["*"]), "RNA5S": ("rrna_5s", ["*"])}


def cmd_pol3(args):
    sortie, vus = [], set()
    with ouvrir(args.fichier) as f:
        for d in csv.DictReader(f, delimiter="\t"):
            if d["arn"] not in POL3:
                continue
            tid, elements = POL3[d["arn"]]
            src = d["source"]
            if src.startswith("HGNC"):
                preuve, reference = "base_de_donnees", "HGNC hgnc_complete_set (groupes de gènes ; table chemin1 partenaires_pol3_arnt_arnr.tsv)"
                source = "HGNC groupe de gènes : " + d.get("lien", "")
            else:
                preuve = "experimentale_humaine"
                m = re.search(r"doi:(\S+)", src)
                reference = m.group(1) if m else src
                source = src + " : " + d.get("lien", "")
            for e in elements:
                cle = (tid, e, d["partenaire"])
                if cle in vus:
                    continue
                vus.add(cle)
                sortie.append(ligne(type_id=tid, element=e, cible=d["partenaire"],
                                    type_cible="proteine_partenaire", nature="liaison",
                                    preuve=preuve, source=source, reference=reference))
    ecrire(sortie)
    print(f"[pol3] {len(sortie)} lignes", file=sys.stderr)


# ---------------------------------------------------------------- GTF
def lire_genes_gtf(chemin):
    genes = []
    with ouvrir(chemin) as f:
        for l in f:
            if l.startswith("#"):
                continue
            c = l.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "gene":
                continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', c[8]))
            genes.append({"chr": c[0], "deb": int(c[3]), "fin": int(c[4]), "brin": c[6],
                          "id": a.get("gene_id", ".").split(".")[0],
                          "type": a.get("gene_type", "."), "nom": a.get("gene_name", ".")})
    return genes


TYPES_PSEUDO = {"transcribed_processed_pseudogene", "transcribed_unprocessed_pseudogene",
                "transcribed_unitary_pseudogene", "translated_processed_pseudogene"}
TYPES_LNC = {"lncRNA", "lincRNA", "antisense", "antisense_RNA", "sense_intronic",
             "sense_overlapping", "processed_transcript", "3prime_overlapping_ncRNA",
             "bidirectional_promoter_lncRNA", "macro_lncRNA", "non_coding"}
RE_PSEUDO_NOM = re.compile(r"^(.*?)\s+pseudogene(?:\s+\d+)?$")


def parent_hgnc(d, par_symbole, par_nom):
    """Gène parent d'un pseudogène selon la nomenclature HGNC ('<parent> pseudogene N')."""
    m = RE_PSEUDO_NOM.match(d["name"])
    if not m:
        return None
    base = m.group(1)
    if base in par_symbole and par_symbole[base]["locus_group"] != "pseudogene":
        return base
    return par_nom.get(base.lower())


def cmd_pseudogenes(args):
    par_symbole, _ = lire_hgnc(args.hgnc)
    par_nom = {d["name"].lower(): s for s, d in par_symbole.items()
               if d["locus_group"] != "pseudogene" and d["status"] == "Approved"}
    par_ensg = {d.get("ensembl_gene_id", ""): d for d in par_symbole.values()
                if d.get("ensembl_gene_id")}
    garder = None
    if args.gtf:
        garder = {g["id"]: g for g in lire_genes_gtf(args.gtf) if g["type"] in TYPES_PSEUDO}
    sortie = []
    for s, d in par_symbole.items():
        if d["locus_group"] != "pseudogene":
            continue
        ensg = d.get("ensembl_gene_id", "") or "."
        if garder is not None and ensg not in garder:
            continue
        p = parent_hgnc(d, par_symbole, par_nom)
        if not p:
            continue
        sortie.append(ligne(type_id="pseudogene_eponge", element=s, element_id=ensg, cible=p,
                            type_cible="gene_parent", nature="origine", preuve="base_de_donnees",
                            source=f"HGNC nom approuvé « {d['name']} » ; biotype GENCODE "
                                   f"{garder[ensg]['type'] if garder else '.'}",
                            reference="HGNC hgnc_complete_set (copie locale) ; GENCODE v27"))
    ecrire(sortie)
    print(f"[pseudogenes] {len(sortie)} lignes", file=sys.stderr)


def cmd_humanin(args):
    par_symbole, _ = lire_hgnc(args.hgnc)
    sortie = []
    for s, d in sorted(par_symbole.items()):
        if re.match(r"^MTRNR2L\d+$", s) and d["name"].startswith("MT-RNR2 like"):
            sortie.append(ligne(type_id="micropeptide_numt_humanin", element=s,
                                element_id=d.get("ensembl_gene_id") or ".", cible="MT-RNR2",
                                type_cible="gene_parent", nature="origine",
                                preuve="base_de_donnees",
                                source=f"HGNC nom approuvé « {d['name']} »",
                                reference="HGNC hgnc_complete_set (copie locale)"))
    ecrire(sortie)
    print(f"[humanin] {len(sortie)} lignes", file=sys.stderr)


def cmd_lncrna_cis(args):
    genes = lire_genes_gtf(args.gtf)
    pc = {}
    for g in genes:
        if g["type"] == "protein_coding":
            pc.setdefault(g["chr"], []).append(g)
    for v in pc.values():
        v.sort(key=lambda g: g["deb"])
    D = args.distance
    sortie = []
    for g in genes:
        if g["type"] not in TYPES_LNC:
            continue
        voisins = pc.get(g["chr"], [])
        tss = g["deb"] if g["brin"] == "+" else g["fin"]
        relations = {}
        amont = aval = None
        for p in voisins:
            if p["fin"] < g["deb"] - D or p["deb"] > g["fin"] + D:
                continue
            chevauche = p["deb"] <= g["fin"] and p["fin"] >= g["deb"]
            ptss = p["deb"] if p["brin"] == "+" else p["fin"]
            if chevauche and p["brin"] != g["brin"]:
                relations.setdefault(p["nom"], "antisens chevauchant")
            elif chevauche:
                relations.setdefault(p["nom"], "même brin chevauchant")
            elif p["brin"] != g["brin"] and abs(ptss - tss) <= 1000:
                relations.setdefault(p["nom"], "bidirectionnel (TSS à ≤ 1 kb)")
            else:
                dist = g["deb"] - p["fin"] if p["fin"] < g["deb"] else p["deb"] - g["fin"]
                if p["fin"] < g["deb"]:
                    if amont is None or dist < amont[0]:
                        amont = (dist, p["nom"])
                else:
                    if aval is None or dist < aval[0]:
                        aval = (dist, p["nom"])
        for v in (amont, aval):
            if v and v[1] not in relations:
                relations[v[1]] = f"plus proche gène codant non chevauchant ({v[0]} pb)"
        for nom, rel in sorted(relations.items()):
            sortie.append(ligne(type_id="lncrna_general", element=g["nom"], element_id=g["id"],
                                cible=nom, type_cible="gene_voisin_cis", nature="cis",
                                preuve="prediction",
                                source=f"GENCODE v27, relation génomique : {rel} ; fenêtre {D} pb",
                                reference="GENCODE v27 (primary_assembly)"))
    ecrire(sortie)
    print(f"[lncrna_cis] {len(sortie)} lignes", file=sys.stderr)


# ---------------------------------------------------------------- miRTarBase
def symbole_mir(precurseur):
    """hsa-mir-181a-1 -> MIR181A1 ; hsa-mir-1-1 -> MIR1-1 ; hsa-let-7a-1 -> MIRLET7A1 (règle HGNC : le tiret
    avant le numéro de copie est gardé seulement après un chiffre)."""
    p = precurseur.lower()
    if not p.startswith("hsa-"):
        return None
    p = p[4:]
    if p.startswith("mir-"):
        pre, reste = "MIR", p[4:]
    elif p.startswith("let-"):
        pre, reste = "MIRLET", p[4:]
    else:
        return None
    reste = re.sub(r"(?<=[a-z])-(?=\d+$)", "", reste)
    return pre + reste.upper()


def precurseurs_de(mature, hairpins):
    """hsa-miR-181a-5p -> précurseurs hsa-mir-181a-1, hsa-mir-181a-2 (hairpin.fa de miRBase)."""
    b = re.sub(r"-(5p|3p)$", "", mature.strip()).rstrip("*").lower()
    return [h for h in hairpins if h == b or re.match(re.escape(b) + r"-\d+$", h)]


def est_humain(v):
    return v.strip() in ("hsa", "Homo sapiens")


def cmd_mirtarbase(args):
    """miRTarBase (CSV ou TSV, gz accepté) : paires miARN humain - gène humain, dédoublonnées (une ligne par
    expérience dans la source). preuve = experimentale_humaine si au moins un « Functional MTI » non « Weak »,
    sinon base_de_donnees. --fonctionnelles-seulement écarte les paires « Weak » (CLIP, séquençage seuls)."""
    hairpins = set()
    with ouvrir(args.hairpin) as f:
        for l in f:
            if l.startswith(">hsa-"):
                hairpins.add(l[1:].split()[0].lower())
    par_symbole, alias = lire_hgnc(args.hgnc)
    paires, lues = {}, 0
    with ouvrir(args.fichier) as f:
        dialecte = "\t" if args.fichier.endswith((".tsv", ".tsv.gz", ".txt")) else ","
        for d in csv.DictReader(f, delimiter=dialecte):
            lues += 1
            if not est_humain(d.get("Species (miRNA)", "hsa")) or not est_humain(d.get("Species (Target Gene)", "hsa")):
                continue
            cle = (d["miRNA"].strip(), d["Target Gene"].strip())
            p = paires.setdefault(cle, {"fort": False, "support": set(), "exp": set(), "pmid": set(), "id": set()})
            sup = d.get("Support Type", "").strip()
            p["fort"] |= "Functional MTI" in sup and "Weak" not in sup
            p["support"].add(sup)
            p["exp"].update(x for x in d.get("Experiments", "").split("//") if x)
            p["pmid"].add(d.get("References (PMID)", ".").strip())
            p["id"].add(d.get("miRTarBase ID", ".").strip())
    sortie, sans_hgnc, faibles = [], set(), 0
    cache = {}
    for (mature, cible), p in paires.items():
        if args.fonctionnelles_seulement and not p["fort"]:
            faibles += 1
            continue
        if mature not in cache:
            cache[mature] = [s for s in (symbole_mir(h) for h in precurseurs_de(mature, hairpins))
                             if s and s in par_symbole]
        elems = cache[mature]
        if not elems:
            sans_hgnc.add(mature)
            continue
        cible = cible if cible in par_symbole else alias.get(cible, cible)
        exp = sorted(p["exp"])
        src = (f"miRTarBase {','.join(sorted(p['id'])[:3])} {mature} ; {' / '.join(sorted(p['support']))} ; "
               f"{len(exp)} méthode(s) : {'//'.join(exp)[:80]}").replace("\t", " ")
        pm = sorted(p["pmid"])
        ref = f"miRTarBase {args.version} ; PMID {','.join(pm[:5])}" + (f" (+{len(pm) - 5})" if len(pm) > 5 else "")
        for e in elems:
            sortie.append(ligne(type_id="mirna_canonique", element=e, cible=cible, type_cible="arnm",
                                nature="repression_graine",
                                preuve="experimentale_humaine" if p["fort"] else "base_de_donnees",
                                source=src, reference=ref))
    ecrire(sortie)
    print(f"[mirtarbase] {lues} lignes lues ; {len(paires)} paires humaines ; {len(sortie)} lignes écrites ; "
          f"{faibles} paires Weak écartées ; {len(sans_hgnc)} miARN matures sans précurseur HGNC "
          f"(ex. {', '.join(sorted(sans_hgnc)[:5])})", file=sys.stderr)


# ---------------------------------------------------------------- générique
def cmd_generique(args):
    sortie = []
    with ouvrir(args.fichier) as f:
        for d in csv.DictReader(f, delimiter=args.sep):
            el = d[args.col_element].strip()
            ci = d[args.col_cible].strip()
            if not el or not ci:
                continue
            pos = d.get(args.col_position, ".") if args.col_position else "."
            src = args.source
            if args.col_source and d.get(args.col_source):
                src += " " + d[args.col_source]
            sortie.append(ligne(type_id=args.type_id, element=el,
                                element_id=d.get(args.col_element_id, ".") if args.col_element_id else ".",
                                cible=ci, type_cible=args.type_cible, nature=args.nature,
                                position_cible=pos or ".", preuve=args.preuve, source=src,
                                reference=args.reference))
    ecrire(sortie)
    print(f"[generique] {len(sortie)} lignes", file=sys.stderr)


# ---------------------------------------------------------------- fusion
def cmd_fusion(args):
    par_symbole, alias = lire_hgnc(args.hgnc) if args.hgnc else ({}, {})
    vues, toutes = set(), []
    for chemin in args.entrees:
        with ouvrir(chemin) as f:
            for d in csv.DictReader(f, delimiter="\t"):
                l = ligne(**{c: d.get(c, ".") for c in COLONNES})
                if l["element_id"] == "." and l["element"] in par_symbole:
                    l["element_id"] = par_symbole[l["element"]].get("ensembl_gene_id") or "."
                if l["cible"] not in par_symbole and l["cible"] in alias:
                    l["source"] += f" ; cible renommée {l['cible']}->{alias[l['cible']]} (HGNC)"
                    l["cible"] = alias[l["cible"]]
                cle = tuple(l[c] for c in COLONNES[:7])
                if cle in vues:
                    continue
                vues.add(cle)
                toutes.append(l)
    toutes.sort(key=lambda l: (l["type_id"], l["element"], l["cible"], l["position_cible"]))
    with open(args.out, "w", encoding="utf-8") as s:
        ecrire(toutes, s)
    print(f"[fusion] {len(toutes)} lignes -> {args.out}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest="cmd", required=True)
    for nom in ("snodb", "curees", "pol3"):
        p = sp.add_parser(nom)
        p.add_argument("fichier")
    p = sp.add_parser("pseudogenes")
    p.add_argument("--hgnc", required=True)
    p.add_argument("--gtf")
    p = sp.add_parser("humanin")
    p.add_argument("--hgnc", required=True)
    p = sp.add_parser("lncrna_cis")
    p.add_argument("--gtf", required=True)
    p.add_argument("--distance", type=int, default=50000)
    p = sp.add_parser("mirtarbase")
    p.add_argument("fichier")
    p.add_argument("--hairpin", required=True)
    p.add_argument("--hgnc", required=True)
    p.add_argument("--version", default="(version à noter)")
    p.add_argument("--fonctionnelles-seulement", action="store_true",
                   help="garder seulement les paires « Functional MTI » (non Weak)")
    p = sp.add_parser("generique")
    p.add_argument("fichier")
    p.add_argument("--sep", default="\t")
    p.add_argument("--type-id", required=True)
    p.add_argument("--col-element", required=True)
    p.add_argument("--col-element-id")
    p.add_argument("--col-cible", required=True)
    p.add_argument("--col-position")
    p.add_argument("--col-source")
    p.add_argument("--type-cible", required=True, choices=sorted(TYPES_CIBLE))
    p.add_argument("--nature", required=True, choices=sorted(NATURES))
    p.add_argument("--preuve", required=True, choices=sorted(PREUVES))
    p.add_argument("--source", required=True)
    p.add_argument("--reference", required=True)
    p = sp.add_parser("fusion")
    p.add_argument("--hgnc")
    p.add_argument("--out", required=True)
    p.add_argument("entrees", nargs="+")
    a = ap.parse_args()
    {"snodb": cmd_snodb, "curees": cmd_curees, "pol3": cmd_pol3, "pseudogenes": cmd_pseudogenes,
     "humanin": cmd_humanin, "lncrna_cis": cmd_lncrna_cis, "mirtarbase": cmd_mirtarbase,
     "generique": cmd_generique, "fusion": cmd_fusion}[a.cmd](a)


if __name__ == "__main__":
    main()
