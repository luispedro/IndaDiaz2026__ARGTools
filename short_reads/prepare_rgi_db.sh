#!/usr/bin/env bash
# Usage:
#   ./prepare_rgi_db.sh <work_dir> [card.json]
#
# Downloads CARD v4.0.0 automatically if [card.json] is omitted.
# Set CARD_VERSION to fetch a different version.
# <work_dir> ends up with a localDB/ -- pass it as RGI_LOCALDB_DIR to run_all_tools.sh.
set -euo pipefail

if [ "$#" -lt 1 ] || [ "$#" -gt 2 ]; then
    echo "Usage: $0 <work_dir> [card.json]" >&2
    exit 1
fi

WORK_DIR="$1"
CARD_JSON_ARG="${2:-}"
CARD_VERSION="${CARD_VERSION:-4.0.0}"
CARD_URL="https://card.mcmaster.ca/download/0/broadstreet-v${CARD_VERSION}.tar.bz2"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RGI_BIN="$SCRIPT_DIR/.pixi/envs/rgi/bin/rgi"

if [ ! -x "$RGI_BIN" ]; then
    echo "error: $RGI_BIN not found -- run 'pixi install' in $SCRIPT_DIR first" >&2
    exit 1
fi

if [ -n "$CARD_JSON_ARG" ]; then
    if [ ! -f "$CARD_JSON_ARG" ]; then
        echo "error: card.json not found: $CARD_JSON_ARG" >&2
        exit 1
    fi
    # resolve to an absolute path before we cd into WORK_DIR below
    CARD_JSON_ARG="$(cd "$(dirname "$CARD_JSON_ARG")" && pwd)/$(basename "$CARD_JSON_ARG")"
fi

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

if [ -n "$CARD_JSON_ARG" ]; then
    echo ">>> using provided card.json: $CARD_JSON_ARG"
    cp "$CARD_JSON_ARG" card.json
else
    echo ">>> downloading CARD v$CARD_VERSION from $CARD_URL"
    curl -fsSL -o card_data.tar.bz2 "$CARD_URL"
    tar xjf card_data.tar.bz2 ./card.json
    rm -f card_data.tar.bz2
fi

echo ">>> rgi card_annotation"
"$RGI_BIN" card_annotation -i card.json > card_annotation.log 2>&1
mapfile -t vers < <(ls card_database_v*.fasta 2>/dev/null | grep -v _all | sed -E 's/card_database_v(.*)\.fasta/\1/')
if [ "${#vers[@]}" -ne 1 ]; then
    echo "error: expected exactly one card_database_v*.fasta (excluding _all), found ${#vers[@]} -- see card_annotation.log" >&2
    exit 1
fi
ver="${vers[0]}"

echo ">>> rgi load (CARD v$ver)"
"$RGI_BIN" clean --local
"$RGI_BIN" load --card_json card.json \
    --card_annotation "card_database_v${ver}.fasta" \
    --card_annotation_all_models "card_database_v${ver}_all.fasta" \
    --local

echo ">>> materializing localDB/ with a throwaway seed run"
printf ">seed\nATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT\n" > _seed.fna
"$RGI_BIN" main -a DIAMOND -i _seed.fna -o _seed_out --local --clean -t contig -n 1
rm -f _seed*

echo ">>> done -- localDB/ ready at $WORK_DIR/localDB (CARD v$ver)"
