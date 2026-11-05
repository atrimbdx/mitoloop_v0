"""
promoter_elements.py - elements of the core promoter and promoter type

Positioned TATA box, initiator and G-quadruplexes on the gene strand, tissue
specificity of expression (tau) and the promoter type (housekeeping, tissue,
developmental; Lenhard et al. 2012) used by promoter_features.py. The same
rules are applied to the variants by the reannot repository
(vm/ajout_mitoloop.py); keep both in step.
"""

from __future__ import annotations

import math
import re

COMP = str.maketrans("ACGTacgt", "TGCAtgca")

# éléments du promoteur minimal, lus sur le brin du gène
TATA = re.compile(r"(?=TATA[AT]A[AT][AG])")          # TATAWAWR (Bucher 1990)
INR = re.compile(r"(?=[CT][CT]A[ACGT][AT][CT][CT])")  # YYANWYY, A = +1
G4 = re.compile(r"(?=(G{3,}[ACGT]{1,7}?G{3,}[ACGT]{1,7}?G{3,}[ACGT]{1,7}?G{3,}))")
C4 = re.compile(r"(?=(C{3,}[ACGT]{1,7}?C{3,}[ACGT]{1,7}?C{3,}[ACGT]{1,7}?C{3,}))")


def oriented(seq, tss, strand, up=40, down=20):
    """séquence de -up à +down du TSS (indice tss dans seq, base +1), sur le
    brin du gène ; indice up = +1 ; '' si la fenêtre sort de seq"""
    if strand == "+":
        a, b = tss - up, tss + down
        return seq[a:b] if a >= 0 and b <= len(seq) else ""
    a, b = tss - down + 1, tss + up + 1
    return seq[a:b].translate(COMP)[::-1] if a >= 0 and b <= len(seq) else ""


def tata_box(o, up=40):
    """boîte TATA positionnée : début entre -35 et -24 (fin vers -20)"""
    return any(up - 35 <= m.start() <= up - 24 for m in TATA.finditer(o))


def inr(o, up=40):
    """initiateur YYANWYY dont le A est à +1 (à 1 pb près, TSS des bases)"""
    return any(up - 1 <= m.start() + 2 <= up + 1 for m in INR.finditer(o))


def g4_spans(seq):
    """G-quadruplex potentiels (G3+N1-7 x4, sur l'un ou l'autre brin) : (début, fin)"""
    out = []
    for rx in (G4, C4):
        out += [(m.start(), m.start() + len(m.group(1))) for m in rx.finditer(seq)]
    return out


def tau(values):
    """spécificité tissulaire (Yanai et al. 2005) sur log2(TPM + 1) : 0 =
    ubiquitaire, 1 = un seul tissu ; None sans expression"""
    x = [math.log2(v + 1) for v in values if v == v]
    m = max(x) if x else 0
    if len(x) < 2 or m <= 0:
        return None
    return sum(1 - v / m for v in x) / (len(x) - 1)


def promoter_type(cls, tata, tau_, polycomb, cage_width=None):
    """ménage, tissu, développement ou indéterminé (Lenhard et al. 2012) :
    cls classe CpG (HCP/ICP/LCP, Weber 2007) ; tata boîte TATA positionnée ;
    tau_ spécificité d'expression ; polycomb part des tissus où le TSS est
    réprimé par Polycomb (ChromHMM ReprPC) ou None ; cage_width largeur du pic
    CAGE (initiation précise si <= 10 pb)"""
    sharp = cage_width is not None and cage_width <= 10
    if cls in ("HCP", "ICP") and polycomb is not None and polycomb >= 0.5:
        return "developpement"
    if cls == "LCP" or (tau_ is not None and tau_ >= 0.8) or \
            ((tata or sharp) and cls != "HCP") or \
            ((tata or sharp) and tau_ is not None and tau_ >= 0.6):
        return "tissu"
    if cls == "HCP" and (tau_ is None or tau_ < 0.5):
        return "menage"
    return "indetermine"
