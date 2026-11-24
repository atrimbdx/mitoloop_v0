#!/usr/bin/env bash
# inventaire_curta.sh - read-only inventory of the MITODARK / mitoloop work on
#                       Curta, before any clean-up.
#
# Purpose      Lists every file under the given roots (or, by default, every
#              directory of $HOME, /gpfs/home/$USER and /scratch/$USER whose
#              name contains mitodark or mitoloop, down to depth 3), and
#              describes the key result files so that each one can be matched
#              to the numbers of the ANR text (AAPG 2027, MitoLoop):
#                1,203 promoters ; 2,193 loops at 907 promoters (reciprocal,
#                >= 3 tissues) ; 3,135 replicable rE2G elements (80 % of the
#                64 skeletal-muscle samples) ; ~20 % of muscle loops overlapping
#                an element (208 / 1,057).
#              Nothing is modified, moved or deleted.
# Outputs      inventaire_<date>/
#                resume.txt     sizes per root and subdirectory, file types,
#                               largest files, git repositories, symbolic
#                               links leaving the roots   <- to send back
#                cles.tsv       key files: path, size, date, lines, header,
#                               distinct genes / promoters / terms  <- to send back
#                fichiers.tsv.gz  every file: path, size, date
#                liens.tsv      symbolic links: link -> target
# Usage        bash inventaire_curta.sh [ROOT ...]
#              (e.g. bash inventaire_curta.sh ~/MITODARK ~/mitoloop_v0 ~/reannot_git)

set -uo pipefail
OUT="inventaire_$(date +%Y%m%d_%H%M)"
mkdir -p "$OUT"

roots=("$@")
if [ ${#roots[@]} -eq 0 ]; then
    for base in "$HOME" "/gpfs/home/${USER:-x}" "/scratch/${USER:-x}"; do
        [ -d "$base" ] || continue
        while IFS= read -r d; do roots+=("$d"); done < <(
            find "$base" -maxdepth 3 -type d \( -iname '*mitodark*' -o -iname '*mitoloop*' \) \
                 -not -path '*/.git/*' -prune 2>/dev/null)
    done
    [ -d "$HOME/reannot_git" ] && roots+=("$HOME/reannot_git")
fi
# unique real paths, nested roots removed
mapfile -t roots < <(for r in "${roots[@]}"; do readlink -f "$r"; done | sort -u |
                     awk '{ for (i in seen) if (index($0, seen[i] "/") == 1) next; seen[NR] = $0; print }')
[ ${#roots[@]} -gt 0 ] || { echo "aucune racine trouvée : donner les dossiers en argument"; exit 1; }
printf '%s\n' "${roots[@]}" > "$OUT/racines.txt"
echo "racines :"; sed 's/^/  /' "$OUT/racines.txt"

# every file (no link followed)
for r in "${roots[@]}"; do
    find "$r" -path '*/.git' -prune -o -type f -printf '%p\t%s\t%TY-%Tm-%Td\n' 2>/dev/null
done | gzip > "$OUT/fichiers.tsv.gz"
for r in "${roots[@]}"; do
    find "$r" -path '*/.git' -prune -o -type l -printf '%p\t%l\n' 2>/dev/null
done > "$OUT/liens.tsv"

python3 - "$OUT" "${roots[@]}" <<'PY'
import collections, csv, gzip, io, os, re, subprocess, sys
out, roots = sys.argv[1], sys.argv[2:]
rows = [l.rstrip("\n").split("\t") for l in gzip.open(f"{out}/fichiers.tsv.gz", "rt")]
rows = [(p, int(s), d) for p, s, d in rows]
hum = lambda n: next(f"{n / 1024 ** k:.1f} {u}" for k, u in ((4, "To"), (3, "Go"), (2, "Mo"), (1, "Ko"), (0, "o")) if n >= 1024 ** k or k == 0)
res = []
res.append(f"Inventaire {out} ({len(rows):,} fichiers, {hum(sum(s for _, s, _ in rows))})\n")
for r in roots:
    sub = [(p, s, d) for p, s, d in rows if p.startswith(r + "/")]
    res.append(f"\n== {r} : {len(sub):,} fichiers, {hum(sum(s for _, s, _ in sub))}, "
               f"dates {min((d for _, _, d in sub), default='-')} .. {max((d for _, _, d in sub), default='-')}")
    by = collections.defaultdict(lambda: [0, 0, "0000"])
    for p, s, d in sub:
        rel = os.path.relpath(p, r).split("/")
        keys = [rel[0] + "/"] if len(rel) > 1 else ["."]
        if len(rel) > 2:
            keys.append("   " + "/".join(rel[:2]) + "/")
        for key in keys:
            v = by[key]; v[0] += 1; v[1] += s; v[2] = max(v[2], d)
    for k, (n, s, d) in sorted(by.items(), key=lambda x: x[0].strip()):
        res.append(f"   {k:55s} {n:8,} fichiers {hum(s):>10s}  dernier {d}")
ext = collections.Counter()
size_ext = collections.Counter()
for p, s, _ in rows:
    e = re.sub(r"^.*?(\.[A-Za-z0-9]+(\.gz|\.bgz)?)$", r"\1", os.path.basename(p)) if "." in os.path.basename(p) else "(sans)"
    ext[e] += 1; size_ext[e] += s
res.append("\n== types de fichiers (nombre, taille)")
for e, n in ext.most_common(25):
    res.append(f"   {e:15s} {n:8,} {hum(size_ext[e]):>10s}")
res.append("\n== 30 plus gros fichiers")
for p, s, d in sorted(rows, key=lambda x: -x[1])[:30]:
    res.append(f"   {hum(s):>10s}  {d}  {p}")
res.append("\n== dépôts git")
for r in roots:
    for dp, dn, fn in os.walk(r):
        if ".git" in dn:
            g = dp
            q = lambda *a: subprocess.run(["git", "-C", g, *a], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True).stdout.strip()
            res.append(f"   {g}\n      remote {q('remote', 'get-url', 'origin') or '-'} ; branche {q('rev-parse', '--abbrev-ref', 'HEAD')} ; "
                       f"commit {q('log', '-1', '--format=%h %cs %s')[:90]} ; modifiés non commités {len(q('status', '--porcelain').splitlines())}")
            dn.remove(".git")
        dn[:] = [x for x in dn if x not in ("data", ".venv-regulome", "__pycache__", "cache")]
res.append("\n== liens symboliques qui sortent des racines (à ne pas casser)")
links = [l.rstrip("\n").split("\t") for l in open(f"{out}/liens.tsv")]
n_out = 0
for l, t in links:
    tt = os.path.normpath(os.path.join(os.path.dirname(l), t))
    if not any(tt == r or tt.startswith(r + "/") for r in roots):
        n_out += 1
        if n_out <= 40:
            res.append(f"   {l} → {t}")
res.append(f"   ({n_out} liens sortants, {len(links)} liens en tout ; liste complète : liens.tsv)")

# key files
PAT = re.compile(r"(loops?|boucles?|re2g|element|overlap|peaks_vs|assoc|libraries|genes\.txt|l1|regions|"
                 r"v4c|peakc|anr|mitoloop|mitodark|aapg|figure)", re.I)
TAB = (".tsv", ".tsv.gz", ".txt", ".csv", ".bed", ".bed.gz", ".csv.gz")
DOC = (".pdf", ".docx", ".doc", ".pptx", ".svg", ".png", ".xlsx")
cles = []
for p, s, d in rows:
    b = os.path.basename(p).lower()
    if "/data/" in p or "/cache/" in p or "/.venv" in p:
        continue
    if not PAT.search(b) or not b.endswith(TAB + DOC):
        continue
    info = {"chemin": p, "taille": s, "date": d, "lignes": ".", "entete": ".", "genes": ".",
            "promoteurs_ou_loci": ".", "termes": "."}
    if b.endswith(TAB) and s < 2_000_000_000:
        try:
            fh = gzip.open(p, "rt", errors="replace") if b.endswith(".gz") else open(p, errors="replace")
            head = fh.readline().rstrip("\n")
            sep = "\t" if head.count("\t") >= head.count(",") else ","
            cols = [c.strip().lstrip("#").lower() for c in head.split(sep)]
            sets = collections.defaultdict(set)
            n = 0
            for ln in fh:
                n += 1
                if n > 5_000_000:
                    break
                f = ln.rstrip("\n").split(sep)
                for i, c in enumerate(cols[:len(f)]):
                    if c in ("gene", "gene_name", "symbol", "target_gene", "targetgene"):
                        sets["genes"].add(f[i])
                    elif c in ("promoter", "locus", "viewpoint", "bait", "tss_id"):
                        sets["promoteurs_ou_loci"].add(f[i])
                    elif c in ("term", "terms", "tissue", "encode_term", "biosample", "called"):
                        sets["termes"].update(x for x in re.split(r"[|;,]", f[i]) if x)
            fh.close()
            info.update(lignes=n, entete=head[:200].replace("\t", " | "),
                        **{k: len(v) for k, v in sets.items()})
        except Exception as ex:          # fichier illisible : on le signale seulement
            info["entete"] = f"(illisible : {ex})"[:200]
    cles.append(info)
cols = ["chemin", "taille", "date", "lignes", "genes", "promoteurs_ou_loci", "termes", "entete"]
with open(f"{out}/cles.tsv", "w") as fo:
    fo.write("\t".join(cols) + "\n")
    for c in sorted(cles, key=lambda x: x["chemin"]):
        fo.write("\t".join(str(c[k]) for k in cols) + "\n")
res.append(f"\n== fichiers clés : {len(cles)} (cles.tsv)")
open(f"{out}/resume.txt", "w").write("\n".join(res) + "\n")
print("\n".join(res[:3]))
PY
echo
echo "Inventaire : $OUT/  (aucun fichier modifié)"
echo "À m'envoyer : $OUT/resume.txt et $OUT/cles.tsv"
