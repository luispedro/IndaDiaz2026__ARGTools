#!/usr/bin/env python
"""
Usage:
    pixi run download-deeparg-db <dir>

Downloads the gaarangoa/deeparg model+database bundle (DeepARG v2) into
<dir>. <dir> is what you then set as DEEPARG_HF_DIR in config.py.

This has to run under the 'deeparg' pixi environment (huggingface_hub comes
in as a dependency of deeparg-modern there) -- the pixi task below is
defined under [feature.deeparg.tasks] so `pixi run download-deeparg-db`
always picks that environment, never jugfile.py's.
"""
import argparse
import os

from huggingface_hub import snapshot_download


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dir", help="destination directory")
    args = parser.parse_args()

    os.makedirs(args.dir, exist_ok=True)
    dest = os.path.abspath(args.dir)
    print(f">>> downloading gaarangoa/deeparg (DeepARG v2 model + database bundle) into {dest}")
    path = snapshot_download(repo_id="gaarangoa/deeparg", repo_type="model", local_dir=dest)
    print(f"done: {path}")
    print(f">>> done -- set DEEPARG_HF_DIR = {dest!r} in config.py")


if __name__ == "__main__":
    main()
