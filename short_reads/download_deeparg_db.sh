#!/usr/bin/env bash
# Usage:
#   ./download_deeparg_db.sh <dir>
#
# <dir> is what you then pass as DEEPARG_HF_DIR to run_all_tools.sh.
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <dir>" >&2
    exit 1
fi

DIR="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$DIR"
DIR="$(cd "$DIR" && pwd)"
echo ">>> downloading gaarangoa/deeparg (DeepARG v2 model + database bundle) into $DIR"
# Must run in the 'deeparg' pixi environment (has huggingface_hub via
# deeparg-modern) -- the default environment has no python dependency and
# silently falls back to whatever python is on PATH, which won't have it.
# --manifest-path makes this independent of the caller's cwd.
pixi run --manifest-path "$SCRIPT_DIR/pixi.toml" --environment deeparg python -c "
from huggingface_hub import snapshot_download
path = snapshot_download(repo_id='gaarangoa/deeparg', repo_type='model', local_dir='$DIR')
print('done:', path)
"
echo ">>> done -- pass DEEPARG_HF_DIR=$DIR to run_all_tools.sh"
