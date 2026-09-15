#!/usr/bin/env bash
# Usage:
#   ./run_all_tools.sh <r1.fastq[.gz]> <r2.fastq[.gz]> <output_dir> [threads]
#
# ONLY_TOOL=fargene|rgi|deeparg runs just one tool instead of all three.
# RGI_LOCALDB_DIR and DEEPARG_HF_DIR must be set (see pixi.toml's
# prepare-rgi-db / download-deeparg-db tasks).
set -euo pipefail

if [ "$#" -lt 3 ] || [ "$#" -gt 4 ]; then
    echo "Usage: $0 <r1.fastq[.gz]> <r2.fastq[.gz]> <output_dir> [threads]" >&2
    exit 1
fi

R1="$1"
R2="$2"
OUTDIR="$3"
THREADS="${4:-4}"
ONLY_TOOL="${ONLY_TOOL:-all}"

if [ ! -f "$R1" ]; then
    echo "error: R1 not found: $R1" >&2
    exit 1
fi
if [ ! -f "$R2" ]; then
    echo "error: R2 not found: $R2" >&2
    exit 1
fi
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [ "$THREADS" -lt 1 ]; then
    echo "error: threads must be a positive integer, got: $THREADS" >&2
    exit 1
fi
case "$ONLY_TOOL" in
    all|fargene|rgi|deeparg) ;;
    *) echo "ONLY_TOOL must be all|fargene|rgi|deeparg, got: $ONLY_TOOL" >&2; exit 1 ;;
esac

RGI_LOCALDB_DIR="${RGI_LOCALDB_DIR:-/EDIT/ME/rgi_work_dir}"
DEEPARG_HF_DIR="${DEEPARG_HF_DIR:-/EDIT/ME/deeparg_hf_bundle}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FARGENE_BIN="$SCRIPT_DIR/.pixi/envs/fargene/bin/fargene"
RGI_BIN="$SCRIPT_DIR/.pixi/envs/rgi/bin/rgi"
DEEPARG_BIN="$SCRIPT_DIR/.pixi/envs/deeparg/bin/deeparg"

check_bin() {
    local bin="$1" tool="$2"
    if [ ! -x "$bin" ]; then
        echo "error: $bin not found -- run 'pixi install' in $SCRIPT_DIR (env: $tool)" >&2
        return 1
    fi
}

# Paths registered here are removed once, on exit, whether we succeed, fail,
# or one tool's failure leaves another's leftovers behind (set -e means a
# failing command can jump straight past any cleanup written after it).
CLEANUP_PATHS=()
cleanup() {
    local p
    for p in "${CLEANUP_PATHS[@]:-}"; do
        [ -n "$p" ] && rm -rf "$p"
    done
}
trap cleanup EXIT

mkdir -p "$OUTDIR"
R1="$(cd "$(dirname "$R1")" && pwd)/$(basename "$R1")"
R2="$(cd "$(dirname "$R2")" && pwd)/$(basename "$R2")"

is_gzip() {
    [ "$(head -c2 "$1" | od -An -tx1 | tr -d ' \n')" = "1f8b" ]
}

decompress_or_copy() {
    local src="$1" dst="$2"
    if is_gzip "$src"; then
        gunzip -c "$src" > "$dst"
    else
        cp "$src" "$dst"
    fi
}

run_fargene() {
    check_bin "$FARGENE_BIN" fargene || return 1
    declare -A models=(
        [class_a]=class_a               [class_b1_b2]=class_b_1_2        [class_b3]=class_b_3
        [class_c]=class_c               [class_d1]=class_d_1             [class_d2]=class_d_2
        [mph]=mph                       [erm_1]=erm_type_a               [erm_2]=erm_type_f
        [tet_enzyme]=tet_enzyme         [tet_rpg]=tet_rpg
        [aph2b]=aminoglycoside_model_g  [aph3p]=aminoglycoside_model_h   [aph6]=aminoglycoside_model_i
        [aac2p]=aminoglycoside_model_a  [aac3_1]=aminoglycoside_model_b  [aac3_2]=aminoglycoside_model_c
        [aac6p_1]=aminoglycoside_model_d [aac6p_2]=aminoglycoside_model_e [aac6p_3]=aminoglycoside_model_f
        [tet_efflux]=tet_efflux         [qnr]=qnr
    )
    mkdir -p "$OUTDIR/fargene"

    # fargene can't read .gz directly -- decompress once into a scratch dir
    # under outdir, reuse across all 22 classes, delete when done (even on
    # failure, via the trap).
    local tmp_dir="$OUTDIR/fargene/.decompressed"
    mkdir -p "$tmp_dir"
    CLEANUP_PATHS+=("$tmp_dir")
    local r1_plain="$tmp_dir/r1.fastq"
    local r2_plain="$tmp_dir/r2.fastq"
    if ! decompress_or_copy "$R1" "$r1_plain" || ! decompress_or_copy "$R2" "$r2_plain"; then
        echo "error: failed to prepare fastq input for fargene" >&2
        return 1
    fi

    # Classes are independent hmmer scans -- don't let one failing class
    # abort the other 21; collect failures and report/exit at the end.
    local failed=()
    for class in "${!models[@]}"; do
        model="${models[$class]}"
        echo ">>> fargene: $class ($model)"
        if ! "$FARGENE_BIN" \
            -i "$r1_plain" "$r2_plain" \
            --hmm-model "$model" \
            --meta \
            -o "$OUTDIR/fargene/$class" \
            -p "$THREADS" \
            --force; then
            echo "warning: fargene failed for class $class ($model)" >&2
            failed+=("$class")
        fi
    done

    if [ "${#failed[@]}" -gt 0 ]; then
        echo "error: fargene failed for classes: ${failed[*]}" >&2
        return 1
    fi
}

run_rgi() {
    check_bin "$RGI_BIN" rgi || return 1
    if [ ! -d "$RGI_LOCALDB_DIR/localDB" ]; then
        echo "error: no localDB/ under RGI_LOCALDB_DIR=$RGI_LOCALDB_DIR -- run" >&2
        echo "  pixi run prepare-rgi-db $RGI_LOCALDB_DIR" >&2
        return 1
    fi
    mkdir -p "$OUTDIR/rgi"
    echo ">>> rgi bwt"
    if ! ( cd "$RGI_LOCALDB_DIR" && \
      "$RGI_BIN" bwt \
          --read_one "$R1" \
          --read_two "$R2" \
          --output_file "$OUTDIR/rgi/sample.bwt" \
          --local \
          --threads "$THREADS" ); then
        echo "error: rgi bwt failed" >&2
        return 1
    fi
}

run_deeparg() {
    check_bin "$DEEPARG_BIN" deeparg || return 1
    mkdir -p "$OUTDIR/deeparg"

    # deeparg derives intermediate filenames (.paired/.unpaired/.merged/
    # .unmerged) by appending suffixes directly onto whatever path it's
    # given -- symlink inputs into outdir/deeparg/ first so those land there
    # instead of next to the original reads, then delete them when done.
    # Prefix with r1_/r2_ so same-named R1/R2 (e.g. both dirs using
    # "reads.fastq.gz") don't collide on a single symlink.
    local r1_link="$OUTDIR/deeparg/r1_$(basename "$R1")"
    local r2_link="$OUTDIR/deeparg/r2_$(basename "$R2")"
    ln -sf "$R1" "$r1_link"
    ln -sf "$R2" "$r2_link"
    CLEANUP_PATHS+=("$r1_link" "$r2_link")

    local hf_flag=()
    if [ -d "$DEEPARG_HF_DIR" ]; then
        hf_flag=(--hf-model-path "$DEEPARG_HF_DIR")
    else
        echo "warning: DEEPARG_HF_DIR=$DEEPARG_HF_DIR doesn't exist yet -- run" >&2
        echo "  pixi run download-deeparg-db $DEEPARG_HF_DIR" >&2
        echo "falling back to a live Hugging Face download instead." >&2
    fi
    echo ">>> deeparg short_reads_pipeline"
    local rc=0
    "$DEEPARG_BIN" short_reads_pipeline \
        --forward_pe_file "$r1_link" \
        --reverse_pe_file "$r2_link" \
        --output_file "$OUTDIR/deeparg/sample" \
        "${hf_flag[@]}" || rc=$?

    rm -f "$OUTDIR"/deeparg/*.paired "$OUTDIR"/deeparg/*.unpaired \
          "$OUTDIR"/deeparg/*.merged "$OUTDIR"/deeparg/*.unmerged \
          "$r1_link" "$r2_link"

    if [ "$rc" -ne 0 ]; then
        echo "error: deeparg short_reads_pipeline failed (exit $rc)" >&2
        return 1
    fi
}

failed_tools=()
run_tool() {
    if ! "run_$1"; then
        failed_tools+=("$1")
    fi
}

case "$ONLY_TOOL" in
    all)     run_tool fargene; run_tool rgi; run_tool deeparg ;;
    fargene) run_tool fargene ;;
    rgi)     run_tool rgi ;;
    deeparg) run_tool deeparg ;;
esac

if [ "${#failed_tools[@]}" -gt 0 ]; then
    echo "error: failed: ${failed_tools[*]} -- see output above" >&2
    exit 1
fi

echo ">>> done -- outputs in $OUTDIR/{fargene,rgi,deeparg}"
