#!/usr/bin/env bash
# run_mcia.sh - the mitochondrial regulome on public data, on a cluster
#               (MCIA Curta)
#
# Purpose      One entry point for every step, so that a run gives a single
#              archive to send back: setup (Python environment), fetch
#              (public inputs; needs internet, i.e. the login node), run
#              (regulome.py pathways and modes, control genes, promoter
#              features, comparison; L1 when the Jiang 2020 per-sample
#              tables are in data/), bundle (tar.gz of logs and results),
#              or all of them in a row.
# Inputs       data/ (filled by fetch); optional, put by hand:
#              data/jiang2020_protein.tsv and config/protein_tissues.tsv
#              (protein-mRNA correlation in modes), data/aneva_vg.tsv (Vg),
#              data/gtex_jiang_reads.gct.gz + data/jiang2020_samples.tsv
#              (L1 score), data/genes_extra.txt (genes added to the set)
# Outputs      .venv-regulome/, data/, genes.txt, controls.txt,
#              results/regulome, results/promoters, results/l1,
#              mitoloop_regulome_<date>.tar.gz (everything to send back)
# Usage        bash run_mcia.sh setup|fetch|run|bundle|all
#              PYTHON=/path/to/python3 to use another interpreter for setup;
#              REUSE=dir to link files already downloaded (default
#              ../mitoshow/data); NO_GENOME=1 to skip the GRCh38 FASTA.
# Method       see README.md, section "Mitochondrial regulome on MCIA"

set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
STEP="${1:-all}"
mkdir -p results logs
LOG="logs/run_$(date +%Y%m%d_%H%M%S)_${STEP}.log"
exec > >(tee -a "$LOG") 2>&1

# input files (overridable, for tests or files kept elsewhere)
MITOCARTA="${MITOCARTA:-data/Human.MitoCarta3.0.xls}"
GTF="${GTF:-data/gencode.gtf.gz}"
TPM="${TPM:-data/gtex_median_tpm.gct.gz}"
EQTL="${EQTL:-data/GTEx_eQTL_independent.tar}"
GAF="${GAF:-data/goa_human.gaf.gz}"
OBO="${OBO:-data/go-basic.obo}"

activate() {
    if [ -f .venv-regulome/bin/activate ]; then
        # shellcheck disable=SC1091
        source .venv-regulome/bin/activate
    fi
}

env_info() {
    echo "== $(date -Is) step=$STEP host=$(hostname) user=${USER:-?} job=${SLURM_JOB_ID:-none}"
    echo "== git $(git rev-parse --short HEAD 2>/dev/null || echo '?') $(git status --porcelain 2>/dev/null | wc -l) local changes"
    echo "== python $(command -v python3) $(python3 --version 2>&1)"
    python3 - <<'PY' 2>&1
import importlib
for m in ("numpy", "pandas", "scipy", "xlrd", "pysam"):
    try:
        print(f"== {m} {importlib.import_module(m).__version__}")
    except Exception as e:
        print(f"== {m} missing ({type(e).__name__})")
PY
}

# a step, timed, with its exit status in the log
step() {
    local name="$1"; shift
    echo "== [$name] $*"
    local t0=$SECONDS
    "$@"
    local rc=$?
    echo "== [$name] status $rc, $((SECONDS - t0)) s"
    return $rc
}

do_setup() {
    # the Python of a MitoShow setup when there is one (the system python3
    # of the login node may be too old)
    local def=python3
    [ -x ../mitoshow/.venv/bin/python3 ] && def=../mitoshow/.venv/bin/python3
    PY="${PYTHON:-$def}"
    echo "== setup with $PY ($($PY --version 2>&1))"
    "$PY" -m venv .venv-regulome || { echo "venv failed: load a Python >= 3.9 module first (module avail python)"; return 1; }
    activate
    python3 -m pip install --upgrade pip >/dev/null
    python3 -m pip install numpy pandas scipy xlrd requests pysam
}

do_fetch() {
    local reuse=()
    for r in ${REUSE:-../mitoshow/data}; do [ -d "$r" ] && reuse+=(--reuse "$r"); done
    local opt=(--datadir data --genes genes.txt --validation)
    [ -n "${NO_GENOME:-}" ] && opt+=(--no-genome)
    step fetch_regulome python3 scripts/fetch_regulome.py "${opt[@]}" ${reuse[@]+"${reuse[@]}"}
    local r1=$?
    step gtex_fetch python3 scripts/gtex_eqtl.py fetch --datadir data
    local r2=$?
    [ -s data/genes_extra.txt ] && { cat genes.txt data/genes_extra.txt | grep -v '^#' | sort -u > genes.tmp && mv genes.tmp genes.txt; echo "== genes.txt + data/genes_extra.txt: $(wc -l < genes.txt) genes"; }
    return $(( r1 || r2 ))
}

do_run() {
    local rc=0 miss=()
    for x in genes.txt "$MITOCARTA" "$GTF" "$TPM" "$EQTL" "$GAF" "$OBO"; do [ -s "$x" ] || miss+=("$x"); done
    if [ ${#miss[@]} -gt 0 ]; then
        echo "== missing inputs: ${miss[*]}"
        echo "   run 'bash run_mcia.sh fetch' on the login node first (see data/MANIFEST.tsv for failed downloads)"
        return 1
    fi
    # 1. functional layer of each gene (MitoPathways; GO inference elsewhere)
    step pathways python3 scripts/regulome.py pathways --mitocarta "$MITOCARTA" ${MITOPATHWAYS:+--mitopathways "$MITOPATHWAYS"} \
        --genes genes.txt --gaf "$GAF" --obo "$OBO" \
        --gtf "$GTF" --out results/regulome || rc=1
    # 2. control genes matched on expression, variability, length, density
    step controls python3 scripts/make_controls.py --gtf "$GTF" \
        --tpm "$TPM" --genes genes.txt --out controls.txt \
        --pairs results/control_pairs.csv --report results/controls_check.tsv || rc=1
    # 3. modes of regulation: variability of the mRNA (tissues, eQTL, Vg,
    #    protein-mRNA correlation of tissue medians when the protein table is there)
    local modes=(--gtf "$GTF" --genes results/regulome/pathways.tsv
                 --tpm "$TPM" --eqtl "$EQTL"
                 --tissue Muscle_Skeletal --out results/regulome)
    [ -s controls.txt ] && modes+=(--controls controls.txt)
    if [ -s data/jiang2020_protein.tsv ]; then
        modes+=(--protein data/jiang2020_protein.tsv)
        [ -s config/protein_tissues.tsv ] && modes+=(--tissue-map config/protein_tissues.tsv)
    fi
    [ -s data/aneva_vg.tsv ] && modes+=(--vg data/aneva_vg.tsv)
    step modes python3 scripts/regulome.py modes "${modes[@]}" || rc=1
    # 4. promoters of the set against the controls: sequence, CpG islands,
    #    motifs, bidirectional promoters
    local prom=(--gtf "$GTF" --genes genes.txt --tpm "$TPM"
                --focus muscle --out results/promoters)
    [ -s controls.txt ] && prom+=(--controls controls.txt)
    [ -s data/GRCh38.fa ] && prom+=(--fasta data/GRCh38.fa)
    [ -s data/cpgIslandExt.txt.gz ] && prom+=(--cpg-islands data/cpgIslandExt.txt.gz)
    [ -s config/methylation.tsv ] && prom+=(--meth-table config/methylation.tsv)
    step promoters python3 scripts/promoter_features.py "${prom[@]}" || rc=1
    # 5. L1, when the matched RNA and protein of the GTEx tissue proteome are there
    local score=(--modes results/regulome/modes.tsv --layer variability)
    if [ -s data/jiang2020_samples.tsv ] && [ -s data/gtex_jiang_reads.gct.gz ] && [ -s data/jiang2020_protein.tsv ]; then
        step l1 python3 scripts/l1_score.py score --rna data/gtex_jiang_reads.gct.gz \
            --protein data/jiang2020_protein.tsv --samples data/jiang2020_samples.tsv \
            --gtf "$GTF" --genes results/regulome/pathways.tsv \
            --pathways results/regulome/pathways.tsv --out results/l1 || rc=1
        if [ -s results/l1/l1.tsv ]; then
            score=(--scores results/l1/l1.tsv)
            local val=(--l1 results/l1/l1.tsv --gtf "$GTF" --out results/l1/validation)
            [ -s data/kopajtich/raw_counts.tsv ] && [ -s data/kopajtich/pairs.tsv ] && \
                val+=(--patients-rna data/kopajtich/raw_counts.tsv
                      --patients-protein data/kopajtich/proteomics_not_normalized.tsv
                      --patients-samples data/kopajtich/pairs.tsv)
            step l1_validate python3 scripts/l1_score.py validate "${val[@]}" || rc=1
        fi
    else
        echo "== L1: data/jiang2020_protein.tsv, data/jiang2020_samples.tsv and data/gtex_jiang_reads.gct.gz not all present: comparison along the variability score of modes"
    fi
    # 6. promoter architecture along the score
    local cmp=("${score[@]}" --expression results/regulome/modes.tsv
               --pathways results/regulome/pathways.tsv --out results/regulome)
    [ -s results/promoters/promoter_features.tsv ] && cmp+=(--promoters results/promoters/promoter_features.tsv)
    [ -s results/promoters/methylation_by_group.tsv ] && cmp+=(--methylation results/promoters/methylation_by_group.tsv)
    step compare python3 scripts/regulome.py compare "${cmp[@]}" || rc=1
    return $rc
}

do_bundle() {
    local f="mitoloop_regulome_$(date +%Y%m%d_%H%M).tar.gz"
    cp "$LOG" logs/current_run.log
    local items=(logs)
    for x in genes.txt controls.txt data/MANIFEST.tsv data/panelapp_genes_all.tsv results; do
        [ -e "$x" ] && items+=("$x")
    done
    # results only: the inputs stay on the cluster
    tar czf "$f" --exclude='*.part' "${items[@]}"
    echo "== bundle: $(pwd)/$f ($(du -h "$f" | cut -f1)) - send this file back"
}

activate
env_info
# the scripts need Python >= 3.8
if [ "$STEP" != setup ] && ! python3 -c 'import sys; sys.exit(sys.version_info < (3, 8))'; then
    echo "== python3 is $(python3 --version 2>&1): too old. Run 'bash run_mcia.sh setup' first"
    echo "   (it creates .venv-regulome; if it fails: module avail python, module load <python 3.x>)"
    exit 1
fi
case "$STEP" in
    setup) do_setup ;;
    fetch) do_fetch ;;
    run) do_run ;;
    bundle) do_bundle ;;
    all) do_fetch; do_run; r=$?; do_bundle; (exit $r) ;;
    *) echo "unknown step $STEP (setup|fetch|run|bundle|all)"; exit 2 ;;
esac
rc=$?
echo "== $(date -Is) end of $STEP (status $rc)"
exit "$rc"
