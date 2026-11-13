#!/usr/bin/env python3
"""
archiver_curta.py - Curta clean-up of the MITODARK / mitoloop work: keep
(1) the version behind the ANR text (AAPG 2027, MitoLoop) and (2) the current
development version, archive everything else.

Run inventaire_curta.sh first, then fill a keep list (garder.tsv, see
garder.exemple.tsv), one line per path:

  anr     PATH [NAME]   result behind the ANR text: COPIED, read-only, into
                        DEST/reference_<date>/ANR_AAPG2027/ (NAME: name of the
                        copy, default the file name)
  dev     PATH [NAME]   final result of the current version (final loops, rE2G
                        elements, overlap...): kept in place and COPIED into
                        DEST/reference_<date>/dev/resultats_finaux/
  garder  PATH          kept in place, not copied (data/ to re-run, ...)
  depot   PATH          git repository of the current version: its .git and
                        its tracked files are kept in place ; its untracked
                        files (old results...) are archived unless listed

Everything else under --racines is archived: one tar.gz per root in
DEST/ancien_<date>/, with the list of its members. Safety rules:
  - default is a dry run: the plan is printed and written (plan.tsv), nothing
    is touched; --appliquer creates the reference and the archives;
  - originals are removed only with --appliquer --supprimer, and only after
    the archive has been read back and every archived path found in it;
  - symbolic links inside kept paths and inside --proteger paths (default
    ~/reannot and ~/reannot_git) are followed: their targets are kept, so that
    no link in use is broken;
  - nothing outside --racines is ever touched; a .git directory is archived
    only with its whole repository.
The ANR version of the code is frozen as a git bundle of --anr-commit
(REPO=COMMIT, repeatable) next to the copied results, with the ANR PDF
(--anr-pdf) and a MANIFEST.tsv (sha256 of every copy).

Usage
  python3 archiver_curta.py --garder garder.tsv --racines ~/MITODARK ~/mitoloop_v0 \\
      --dest ~/mitoloop_archives --anr-pdf MitoLoop_final.pdf \\
      --anr-commit ~/mitoloop_v0=fb64088            # dry run
  ... --appliquer                                   # reference + archives
  ... --appliquer --supprimer                       # and remove the originals
Python 3 standard library ; tar and git.
"""

import argparse
import datetime
import hashlib
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path


def real(p):
    return os.path.realpath(os.path.expanduser(p))


def under(p, q):
    """p is q or inside q"""
    return p == q or p.startswith(q.rstrip("/") + "/")


def size_of(p):
    if os.path.islink(p) or not os.path.isdir(p):
        try:
            return os.lstat(p).st_size
        except OSError:
            return 0
    n = 0
    for dp, dn, fn in os.walk(p):
        for f in fn:
            try:
                n += os.lstat(os.path.join(dp, f)).st_size
            except OSError:
                pass
    return n


def hum(n):
    for k, u in ((4, "To"), (3, "Go"), (2, "Mo"), (1, "Ko")):
        if n >= 1024 ** k:
            return f"{n / 1024 ** k:.1f} {u}"
    return f"{n} o"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_keep(path):
    out = []
    for n, ln in enumerate(open(path), 1):
        ln = ln.rstrip("\n")
        if not ln.strip() or ln.lstrip().startswith("#"):
            continue
        f = ln.split("\t")
        if f[0] not in ("anr", "dev", "garder", "depot") or len(f) < 2:
            sys.exit(f"{path}:{n} : attendu 'anr|dev|garder|depot<TAB>chemin[<TAB>nom]'")
        if not os.path.lexists(os.path.expanduser(f[1])):
            sys.exit(f"{path}:{n} : {f[1]} n'existe pas")
        lp = os.path.abspath(os.path.expanduser(f[1]))       # le chemin lui-même (lien éventuel)
        out.append((f[0], real(f[1]), f[2] if len(f) > 2 and f[2] else os.path.basename(lp), lp))
    return out


def link_targets(paths, roots):
    """targets (inside the roots) of the symbolic links found under paths"""
    out = {}
    for top in paths:
        if os.path.islink(top):
            items = [top]
        elif os.path.isdir(top):
            items = []
            for dp, dn, fn in os.walk(top):
                items += [os.path.join(dp, x) for x in dn + fn if os.path.islink(os.path.join(dp, x))]
        else:
            items = []
        for l in items:
            t = real(l)
            if any(under(t, r) for r in roots):
                out.setdefault(t, l)
    return out


def plan(roots, keep):
    """units to archive: maximal paths under the roots that hold nothing kept"""
    units = []

    def walk(p):
        if any(under(p, k) for k in keep):
            return                                   # kept (or inside a kept path)
        if not any(under(k, p) for k in keep):
            units.append(p)                          # nothing kept below: archive whole
            return
        if os.path.isdir(p) and not os.path.islink(p):
            for x in sorted(os.listdir(p)):
                walk(os.path.join(p, x))
    for r in roots:
        walk(r)
    # a .git directory goes only with its whole repository
    bad = [u for u in units if os.path.basename(u) == ".git"]
    for u in bad:
        units.remove(u)
    return units, bad


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--garder", required=True)
    p.add_argument("--racines", nargs="+", required=True)
    p.add_argument("--dest", required=True)
    p.add_argument("--proteger", nargs="*", default=["~/reannot", "~/reannot_git"])
    p.add_argument("--anr-pdf", default="")
    p.add_argument("--anr-commit", action="append", default=[], help="DEPOT=COMMIT, répétable")
    p.add_argument("--appliquer", action="store_true")
    p.add_argument("--supprimer", action="store_true")
    a = p.parse_args()
    if a.supprimer and not a.appliquer:
        sys.exit("--supprimer demande --appliquer")
    today = datetime.date.today().strftime("%Y%m%d")
    roots = [real(r) for r in a.racines if os.path.isdir(os.path.expanduser(r))]
    dest = real(a.dest)
    prot = [real(x) for x in a.proteger if os.path.exists(os.path.expanduser(x))]
    if any(under(dest, r) for r in roots):
        sys.exit(f"--dest {dest} est dans une racine : choisir un dossier à part")
    for r in roots:
        if any(under(r, q) or under(q, r) for q in prot):
            sys.exit(f"racine {r} et dossier protégé se recouvrent : retirer l'un des deux")
    keep_rows = read_keep(a.garder)
    for _c, k, _n, _l in keep_rows:
        if not any(under(k, r) for r in roots):
            print(f"note : {k} est hors des racines (gardé de toute façon)")
    keep = {k for c, k, _n, _l in keep_rows if c != "depot"} | \
        {l for c, _k, _n, l in keep_rows if c != "depot"}
    for c, k, _n, _l in keep_rows:                 # dépôt : .git et fichiers suivis
        if c == "depot":
            res = subprocess.run(["git", "-C", k, "ls-files", "-z"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
            if res.returncode:
                sys.exit(f"{k} n'est pas un dépôt git")
            keep |= {os.path.join(k, ".git")} | {os.path.join(k, x) for x in res.stdout.split("\0") if x}
    # targets of links in kept paths, protected paths and the roots themselves
    lt = link_targets(sorted(keep) + prot, roots)
    extra = {t: l for t, l in lt.items() if not any(under(t, k) for k in keep)}
    keep |= set(extra)
    units, gits = plan(roots, sorted(keep))
    sizes = {u: size_of(u) for u in units}

    print("\n== gardé")
    for c, k, n, _l in keep_rows:
        print(f"   {c:7s} {hum(size_of(k)):>10s}  {k}" + (f"  -> {n}" if c in ("anr", "dev") else ""))
    for t, l in sorted(extra.items()):
        print(f"   lien    {hum(size_of(t)):>10s}  {t}   (cible de {l})")
    print(f"\n== à archiver : {len(units)} éléments, {hum(sum(sizes.values()))}")
    for u in units:
        print(f"   {hum(sizes[u]):>10s}  {u}")
    for g in gits:
        print(f"   (non archivé : {g}, son dépôt garde des fichiers ; à traiter à la main)")
    os.makedirs(dest, exist_ok=True)
    with open(os.path.join(dest, f"plan_{today}.tsv"), "w") as fo:
        fo.write("action\tchemin\ttaille\n")
        for c, k, _n, _l in keep_rows:
            fo.write(f"garder_{c}\t{k}\t{size_of(k)}\n")
        for t in sorted(extra):
            fo.write(f"garder_lien\t{t}\t{size_of(t)}\n")
        for u in units:
            fo.write(f"archiver\t{u}\t{sizes[u]}\n")
    print(f"\nplan : {dest}/plan_{today}.tsv")
    if not a.appliquer:
        print("SIMULATION : rien n'a été fait. Relancer avec --appliquer.")
        return

    # 1. reference : ANR (copied, read-only) and final results of the dev version
    ref = Path(dest) / f"reference_{today}"
    anr, dev = ref / "ANR_AAPG2027", ref / "dev" / "resultats_finaux"
    man = []
    for c, k, n, _l in keep_rows:
        if c in ("garder", "depot"):
            continue
        d = (anr / "resultats") if c == "anr" else dev
        d.mkdir(parents=True, exist_ok=True)
        tgt = d / n
        tgt.parent.mkdir(parents=True, exist_ok=True)
        if os.path.isdir(k):
            if tgt.exists():
                shutil.rmtree(str(tgt))
            shutil.copytree(k, str(tgt), symlinks=False)       # Python 3.6 : pas de dirs_exist_ok
            files = [x for x in tgt.rglob("*") if x.is_file()]
        else:
            shutil.copy2(k, tgt)
            files = [tgt]
        for f in files:
            man.append((c, k, str(f), f.stat().st_size, sha256(f)))
    if a.anr_pdf:
        anr.mkdir(parents=True, exist_ok=True)
        t = anr / os.path.basename(a.anr_pdf)
        shutil.copy2(os.path.expanduser(a.anr_pdf), t)
        man.append(("anr", real(a.anr_pdf), str(t), t.stat().st_size, sha256(t)))
    for spec in a.anr_commit:
        repo, commit = spec.split("=", 1)
        repo = real(repo)
        anr.mkdir(parents=True, exist_ok=True)
        subprocess.run(["git", "-C", repo, "tag", "-f", "anr-aapg2027", commit], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        b = anr / f"{os.path.basename(repo)}_anr-aapg2027.bundle"
        subprocess.run(["git", "-C", repo, "bundle", "create", str(b), "anr-aapg2027"], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        full = subprocess.run(["git", "-C", repo, "rev-parse", commit], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True).stdout.strip()
        man.append(("anr", f"{repo}@{full}", str(b), b.stat().st_size, sha256(b)))
    if man:
        with open(ref / "MANIFEST.tsv", "w") as fo:
            fo.write("categorie\tsource\tcopie\ttaille\tsha256\n")
            for r in man:
                fo.write("\t".join(map(str, r)) + "\n")
        (ref / "LISEZMOI.txt").write_text(
            "ANR_AAPG2027/ : résultats et code (git bundle, tag anr-aapg2027) derrière le texte\n"
            "déposé à l'ANR (AAPG 2027, MitoLoop), en lecture seule.\n"
            "dev/resultats_finaux/ : listes finales de la version en développement (boucles,\n"
            "éléments rE2G, recouvrement), copiées ; les originaux restent en place.\n"
            "MANIFEST.tsv : source et sha256 de chaque copie. Restaurer le code :\n"
            "  git clone <bundle> mitoloop_anr && cd mitoloop_anr && git checkout anr-aapg2027\n")
        if anr.exists():
            for dp, dn, fn in os.walk(anr):
                for x in fn + dn:
                    q = os.path.join(dp, x)
                    os.chmod(q, os.stat(q).st_mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
        print(f"référence : {ref} ({len(man)} fichiers, MANIFEST.tsv)")

    # 2. archives, one per root, read back before any removal
    arch = Path(dest) / f"ancien_{today}"
    arch.mkdir(parents=True, exist_ok=True)
    for r in roots:
        us = [u for u in units if under(u, r)]
        if not us:
            continue
        name = r.strip("/").replace("/", "_")
        lst = arch / f"{name}.liste.txt"
        tgz = arch / f"{name}.tar.gz"
        lst.write_text("".join(os.path.relpath(u, "/") + "\n" for u in us))
        print(f"archive : {tgz} ({len(us)} éléments, {hum(sum(sizes[u] for u in us))})")
        subprocess.run(["tar", "-czf", str(tgz), "-C", "/", "--files-from", str(lst)], check=True)
        members = set(subprocess.run(["tar", "-tzf", str(tgz)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                                     check=True).stdout.splitlines())
        members = {m.rstrip("/") for m in members}
        missing = [u for u in us if os.path.relpath(u, "/") not in members]
        if missing:
            print(f"   ARCHIVE INCOMPLÈTE : {len(missing)} éléments absents ({missing[:3]}) ; rien supprimé")
            continue
        with open(arch / f"{name}.contenu.txt", "w") as fo:
            fo.write("\n".join(sorted(members)) + "\n")
        print(f"   relue : {len(members):,} entrées, tous les éléments présents")
        if a.supprimer:
            for u in us:
                if os.path.isdir(u) and not os.path.islink(u):
                    shutil.rmtree(u)
                else:
                    os.remove(u)
            print(f"   originaux supprimés ({len(us)} éléments)")
    print("terminé.")


if __name__ == "__main__":
    main()
