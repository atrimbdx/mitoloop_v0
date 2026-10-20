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
#              data/aneva_vg.tsv (Vg),
#              data/*mmc3.xlsx + data/*mmc4.xlsx (Tables S2 and S3 of Jiang
#              et al. 2020: L1 score, converted by jiang2020_tables.py),
#              data/genes_extra.txt (genes added to the set)
# Outputs      .venv-regulome/, data/, genes.txt, controls.txt,
#              results/regulome, results/promoters, results/l1,
#              mitoloop_regulome_<date>.tar.gz (everything to send back)
# Usage        bash run_mcia.sh setup|selftest|fetch|run|reannot|bundle|all
#              (reannot: the MitoLoop package of the reannot pipeline, in
#              transfert_mitoloop/, see reannot/README_MitoLoop.md)
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
for m in ("numpy", "pandas", "scipy", "xlrd", "openpyxl", "pysam"):
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
    python3 -m pip install numpy pandas scipy xlrd openpyxl requests pysam
}

do_fetch() {
    local reuse=()
    for r in ${REUSE:-../mitoshow/data}; do [ -d "$r" ] && reuse+=(--reuse "$r"); done
    # GTEx and GENCODE first: the gene set is written under the GTF symbols
    step gtex_fetch python3 scripts/gtex_eqtl.py fetch --datadir data
    local r2=$?
    local opt=(--datadir data --genes genes.txt --validation --promoter-annotations --gtf "$GTF")
    [ -n "${NO_GENOME:-}" ] && opt+=(--no-genome)
    step fetch_regulome python3 scripts/fetch_regulome.py "${opt[@]}" ${reuse[@]+"${reuse[@]}"}
    local r1=$?
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
    # 0. Tables S2 (mmc3) and S3 (mmc4) of Jiang et al. 2020, when put in
    #    data/: per-sample tables for L1, per-tissue protein for modes
    do_jiang || rc=1
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
    # experimental annotations (config/promoter_annotations.tsv), those present
    if [ -s config/promoter_annotations.tsv ]; then
        while IFS=$'\t' read -r name file mode win _src; do
            case "$name" in ''|'#'*|name) continue ;; esac
            [ -s "$file" ] || continue
            [ "$win" = - ] && win=""
            prom+=(--annotation "$name=$file:$mode${win:+:$win}")
        done < config/promoter_annotations.tsv
    fi
    step promoters python3 scripts/promoter_features.py "${prom[@]}" || rc=1
    # 5. L1, when the matched RNA and protein of the GTEx tissue proteome are
    #    there (Tables S2 and S3 of Jiang et al. 2020, converted above)
    local score=(--modes results/regulome/modes.tsv --layer variability)
    if [ -s data/jiang2020_samples.tsv ] && [ -s data/jiang2020_rna_samples.tsv ] && [ -s data/jiang2020_protein_samples.tsv ]; then
        local l1=(--rna data/jiang2020_rna_samples.tsv --protein data/jiang2020_protein_samples.tsv
                  --samples data/jiang2020_samples.tsv --gtf "$GTF"
                  --genes results/regulome/pathways.tsv
                  --pathways results/regulome/pathways.tsv --out results/l1)
        l1+=(--rna-units tpm)
        [ -s controls.txt ] && l1+=(--controls controls.txt)
        # ischemic time, when every sample has one (GTEx attributes)
        awk -F'\t' 'NR == 1 {for (i = 1; i <= NF; i++) if ($i == "SMTSISCH") c = i; next}
                    !c || $c == "" {bad = 1} END {exit !(c && !bad && NR > 1)}' \
            data/jiang2020_samples.tsv && l1+=(--covariates SMTSISCH)
        step l1 python3 scripts/l1_score.py score "${l1[@]}" || rc=1
        # sensitivity: without the mitochondrial content adjustment, which
        # only the gene set receives
        step l1_raw python3 scripts/l1_score.py score "${l1[@]/results\/l1/results/l1_raw}" \
            --no-content-adjustment || rc=1
        if [ -s results/l1/l1.tsv ]; then
            # genes ranked on L1 at equal mRNA variation (score_resid)
            score=(--scores results/l1/l1.tsv --score-column score_resid)
            local val=(--l1 results/l1/l1.tsv --gtf "$GTF" --out results/l1/validation) ext=()
            [ -s results/l1_raw/l1.tsv ] && val+=(--replicate results/l1_raw/l1.tsv)
            # protein buffering of copy number in tumours (Goncalves et al.
            # 2017): buffered proteins should have a low L1
            [ -s data/attenuation.csv ] && ext+=(--external attenuation=data/attenuation.csv:attenuation:-)
            # patient fibroblasts (Kopajtich et al.): does a collapse of the
            # mRNA reach the protein, more for genes of high L1?
            if [ -s data/kopajtich/raw_data/raw_counts.tsv ] && \
                    [ ! data/kopajtich/patients_samples.tsv -nt data/kopajtich/raw_data/raw_counts.tsv ]; then
                step kopajtich python3 scripts/kopajtich_tables.py --dir data/kopajtich/raw_data \
                    --gtf "$GTF" --out data/kopajtich --report results/kopajtich || rc=1
            fi
            [ -s data/kopajtich/patients_samples.tsv ] && \
                ext+=(--patients-rna data/kopajtich/patients_rna.tsv
                      --patients-protein data/kopajtich/patients_protein.tsv
                      --patients-samples data/kopajtich/patients_samples.tsv)
            step l1_validate python3 scripts/l1_score.py validate "${val[@]}" ${ext[@]+"${ext[@]}"} || rc=1
            # the same validation for the L1 without content adjustment, the
            # one used to annotate variants (reannot)
            if [ -s results/l1_raw/l1.tsv ]; then
                step l1_raw_validate python3 scripts/l1_score.py validate --l1 results/l1_raw/l1.tsv \
                    --gtf "$GTF" --out results/l1_raw/validation ${ext[@]+"${ext[@]}"} || rc=1
            fi
        fi
    else
        echo "== L1: no Jiang 2020 tables (Tables S2 and S3, *mmc3.xlsx and *mmc4.xlsx, in data/): comparison along the variability score of modes"
    fi
    # 6. promoter architecture along the score
    local cmp=("${score[@]}" --expression results/regulome/modes.tsv
               --pathways results/regulome/pathways.tsv --out results/regulome)
    [ -s results/promoters/promoter_features.tsv ] && cmp+=(--promoters results/promoters/promoter_features.tsv)
    [ -s results/promoters/methylation_by_group.tsv ] && cmp+=(--methylation results/promoters/methylation_by_group.tsv)
    step compare python3 scripts/regulome.py compare "${cmp[@]}" || rc=1
    return $rc
}

# Tables S2 and S3 of Jiang et al. 2020 (Cell 183:269), put by hand in data/
# (or JIANG_DIR), to the three tables of l1_score.py; again when the
# spreadsheets are newer than the tables
do_jiang() {
    local dir="${JIANG_DIR:-data}" s2 s3
    s2=$(ls "$dir"/*mmc3.xlsx 2>/dev/null | head -1)
    s3=$(ls "$dir"/*mmc4.xlsx 2>/dev/null | head -1)
    [ -n "$s2" ] && [ -n "$s3" ] || return 0
    if [ -s data/jiang2020_samples.tsv ] && [ data/jiang2020_samples.tsv -nt "$s2" ] \
            && [ data/jiang2020_samples.tsv -nt "$s3" ]; then
        return 0
    fi
    python3 -c 'import openpyxl' 2>/dev/null || python3 -m pip install openpyxl \
        || { echo "== openpyxl missing: bash run_mcia.sh setup on the login node"; return 1; }
    step jiang2020 python3 scripts/jiang2020_tables.py --s2 "$s2" --s3 "$s3" --dir "$dir" \
        --gtf "$GTF" --tpm "$TPM" --attributes data/GTEx_v8_SampleAttributesDS.txt \
        --download-attributes --out data --report results/jiang2020
}

# the MitoLoop package of the reannot pipeline (VEP, SpliceAI, MORFEE,
# AlphaGenome on an offline VM): gene table with the L1 class, regulatory
# regions of the gene set, annotations cut to these regions, JASPAR motifs,
# the VM scripts; then a check of ajout_mitoloop.py on variants drawn from
# the regions. Written to transfert_mitoloop/, to copy to the VM.
do_reannot() {
    local T=transfert_mitoloop P=data/promoters_exp rc=0
    [ -s results/l1_raw/l1.tsv ] || { echo "== reannot: results/l1_raw/l1.tsv missing (run first)"; return 1; }
    mkdir -p results/reannot
    local reg=(--gtf "$GTF" --genes genes.txt --out results/reannot/regions.bed)
    [ -s results/re2g/re2g_elements.tsv ] && reg+=(--elements results/re2g/re2g_elements.tsv --min-prevalence 5)
    [ -s results/loops/loops.tsv ] && reg+=(--loops results/loops/loops.tsv)
    step reannot_regions python3 scripts/vcf_annotate.py regions "${reg[@]}" || return 1
    rm -rf "$T"; mkdir -p "$T"
    local pk=(--genes genes.txt --l1 results/l1_raw/l1.tsv
              --fibroblast results/l1/validation/patient_genes.tsv
              --regions results/reannot/regions.bed --out "$T"
              --epd "$P/Hs_EPDnew_006_hg38.bed"
              --cage "$P/hg38_fair+new_CAGE_peaks_phase1and2.bed.gz"
              --chromhmm "muscle=$P/E108_15_coreMarks_hg38lift_mnemonics.bed.gz"
              --chromhmm "brain=$P/E073_15_coreMarks_hg38lift_mnemonics.bed.gz"
              --jaspar "$P/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt")
    local f tf
    for f in "$P"/remap2022_*_nr_macs2_hg38_v1_0.bed.gz; do
        [ -s "$f" ] || continue
        tf=$(basename "$f"); tf=${tf#remap2022_}; tf=${tf%%_*}
        pk+=(--remap "$tf=$f")
    done
    step reannot_package python3 reannot/mitoloop_paquet.py "${pk[@]}" || return 1
    cp reannot/ajout_mitoloop.py reannot/07_maj_mitoloop_vm.sh reannot/README_MitoLoop.md "$T/"
    echo "mitoloop_v0 $(git rev-parse --short HEAD 2>/dev/null || echo '?') $(date +%Y-%m-%d)" > "$T/VERSION"
    ( cd "$T" && sha256sum $(find . -type f ! -name SHA256SUMS | sed 's|^\./||' | sort) > SHA256SUMS )
    # check: 300 SNV drawn from the promoters, as a reannot table
    if [ -s data/GRCh38.fa ]; then
        python3 - "$T" <<'PY'
import gzip, random, sys
import pysam
fa = pysam.FastaFile("data/GRCh38.fa")
rows = [l.split("\t") for l in gzip.open(f"{sys.argv[1]}/regions.bed.gz", "rt")]
rows = [r for r in rows if r[4].strip() == "promoter" and r[0] in fa.references]
random.seed(1)
comp = {"A": "G", "G": "A", "C": "T", "T": "C"}
with open("results/reannot/check.final.tsv", "w") as fh:
    fh.write("transcrit\tvariant_entree\tchrom\tpos\tref\talt\tgene\n")
    for r in random.sample(rows, min(300, len(rows))):
        p = random.randint(int(r[1]) + 1, int(r[2]))
        b = fa.fetch(r[0], p - 1, p).upper()
        if b in comp:
            fh.write(f".\t.\t{r[0]}\t{p}\t{b}\t{comp[b]}\t{r[3]}\n")
PY
        step reannot_check python3 "$T/ajout_mitoloop.py" --data "$T" --fasta data/GRCh38.fa \
            results/reannot/check.final.tsv || rc=1
        [ -s results/reannot/check.final.mitoloop.tsv ] && awk -F'\t' 'NR > 1 { n[$(NF-1)]++ } END { for (k in n) print "   " n[k] " variants with " k " argument(s)" }' \
            results/reannot/check.final.mitoloop.tsv | sort -k4n
    fi
    # validation: do the arguments point to the causal variants of the
    # fine-mapped GTEx eQTL of the mitochondrial genes?
    if [ -s data/GTEx_v8_finemapping_DAPG.txt.gz ] && [ -s data/GRCh38.fa ]; then
        step reannot_validation python3 reannot/validation_eqtl.py --data "$T" \
            --dapg data/GTEx_v8_finemapping_DAPG.txt.gz --fasta data/GRCh38.fa \
            --out results/reannot/validation_eqtl.tsv || rc=1
    fi
    cp "$T/MANIFEST.tsv" results/reannot/MANIFEST.tsv
    echo "== reannot: $(pwd)/$T ($(du -sh "$T" | cut -f1)) - copy to the VM, then bash 07_maj_mitoloop_vm.sh"
    return $rc
}

# the tests of the scripts used here, on synthetic data: checks the Python
# environment of the cluster in a few seconds, before the real run
do_selftest() {
    local rc=0
    for t in tests/test_regulome.sh tests/test_l1.sh tests/test_promoters.sh tests/test_eqtl_controls.sh tests/test_jiang.sh tests/test_kopajtich.sh tests/test_reannot.sh; do
        step selftest bash "$t" || rc=1
    done
    return $rc
}

do_bundle() {
    local f="mitoloop_regulome_$(date +%Y%m%d_%H%M).tar.gz"
    cp "$LOG" logs/current_run.log
    local items=(logs)
    for x in genes.txt controls.txt data/MANIFEST.tsv data/panelapp_genes_all.tsv data/genes_renamed.tsv data/kopajtich/CONTENTS.tsv results; do
        [ -e "$x" ] && items+=("$x")
    done
    # the converted Jiang 2020 tables (derived, ~15 MB compressed): further
    # analyses without another run
    for x in data/jiang2020_samples.tsv data/jiang2020_protein.tsv data/jiang2020_rna_samples.tsv data/jiang2020_protein_samples.tsv; do
        [ -s "$x" ] && items+=("$x")
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
    selftest) do_selftest ;;
    fetch) do_fetch ;;
    run) do_run ;;
    bundle) do_bundle ;;
    all) do_fetch; do_run; r=$?; do_reannot || r=1; do_bundle; (exit $r) ;;
    reannot) do_reannot ;;
    *) echo "unknown step $STEP (setup|selftest|fetch|run|reannot|bundle|all)"; exit 2 ;;
esac
rc=$?
echo "== $(date -Is) end of $STEP (status $rc)"
exit "$rc"
