"""
Runs fargene/rgi/deeparg over the samples listed in config.py.

    pixi run jug-status     # what's done / pending
    pixi run jug-execute    # run everything config.py asks for

Safe to re-run: jug skips tasks it already finished, even if you add new
samples or edit shared settings in between runs (only the affected tasks'
hashes change).

Each tool lives in its own pixi environment (see pixi.toml) with its own
binary -- this process (the 'jug' environment) only needs jug itself and
shells out via `pixi run -e <tool> <tool> ...` for the actual work, the same
way the old run_all_tools.sh / prepare_rgi_db.sh did.
"""
import glob
import gzip
import os
import re
import shutil
import subprocess
import sys

from jug import TaskGenerator

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import config

FARGENE_MODELS = {
    "class_a": "class_a", "class_b1_b2": "class_b_1_2", "class_b3": "class_b_3",
    "class_c": "class_c", "class_d1": "class_d_1", "class_d2": "class_d_2",
    "mph": "mph", "erm_1": "erm_type_a", "erm_2": "erm_type_f",
    "tet_enzyme": "tet_enzyme", "tet_rpg": "tet_rpg",
    "aph2b": "aminoglycoside_model_g", "aph3p": "aminoglycoside_model_h", "aph6": "aminoglycoside_model_i",
    "aac2p": "aminoglycoside_model_a", "aac3_1": "aminoglycoside_model_b", "aac3_2": "aminoglycoside_model_c",
    "aac6p_1": "aminoglycoside_model_d", "aac6p_2": "aminoglycoside_model_e", "aac6p_3": "aminoglycoside_model_f",
    "tet_efflux": "tet_efflux", "qnr": "qnr",
}

VALID_TOOLS = {"fargene", "rgi", "deeparg"}


def pixi_run(tool, args, **kwargs):
    """Run <tool> inside its own pixi environment, installing it if needed.

    This process runs under the 'jug' environment, so it can't see the tools;
    `pixi run -e <tool>` activates theirs, which also lets fargene/rgi/deeparg
    find the sibling binaries they shell out to (deeparg -> trimmomatic, rgi
    -> bamtools). <tool> is passed as a plain command rather than as a pixi
    task because tasks always run from the workspace root, while a plain
    command keeps the cwd -- which `rgi ... --local` needs, as it looks for
    localDB/ there.

    pixi locates pixi.toml by searching up from the cwd, so every directory a
    tool is run from has to sit inside the workspace -- that is why
    config.RGI_LOCALDB_DIR defaults to rgi_db/ right here (gitignored).
    """
    return subprocess.run(["pixi", "run", "-e", tool, tool, *args],
                          check=True, **kwargs)


def is_gzip(path):
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def decompress_or_copy(src, dst):
    if is_gzip(src):
        with gzip.open(src, "rb") as fin, open(dst, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copyfile(src, dst)


@TaskGenerator
def prepare_fastq_pair(r1, r2, outdir):
    for src in (r1, r2):
        if not os.path.isfile(src):
            raise RuntimeError(f"reads file not found: {src}")
    # Shared by all fargene classes below -- decompressed once per sample
    # (fargene can't read .gz directly). Left in place afterwards (not
    # cleaned up): a later jug run adding a class would otherwise get this
    # task's cached path back with nothing on disk behind it.
    tmp_dir = os.path.join(outdir, "fargene", ".decompressed")
    os.makedirs(tmp_dir, exist_ok=True)
    r1_plain = os.path.join(tmp_dir, "r1.fastq")
    r2_plain = os.path.join(tmp_dir, "r2.fastq")
    decompress_or_copy(r1, r1_plain)
    decompress_or_copy(r2, r2_plain)
    return r1_plain, r2_plain


@TaskGenerator
def run_fargene_class(fastq_pair, outdir, class_name, model, threads):
    r1_plain, r2_plain = fastq_pair
    class_outdir = os.path.join(outdir, "fargene", class_name)
    pixi_run("fargene", [
        "-i", r1_plain, r2_plain,
        "--hmm-model", model,
        "--meta",
        "-o", class_outdir,
        "-p", str(threads),
        "--force",
    ])
    return class_outdir


@TaskGenerator
def prepare_rgi_db(work_dir, card_json, card_version):
    # Returns work_dir itself (not work_dir/localDB) -- `rgi ... --local`
    # expects to be run with a localDB/ subdirectory *under* its cwd, so
    # run_rgi below needs to cd into work_dir, one level above localDB/.
    os.makedirs(work_dir, exist_ok=True)
    localdb = os.path.join(work_dir, "localDB")
    if os.path.isdir(localdb):
        # Already built, e.g. by a pre-jug run -- nothing to do.
        return work_dir

    if card_json:
        if not os.path.isfile(card_json):
            raise RuntimeError(f"card.json not found: {card_json}")
        shutil.copyfile(card_json, os.path.join(work_dir, "card.json"))
    else:
        url = f"https://card.mcmaster.ca/download/0/broadstreet-v{card_version}.tar.bz2"
        tarball = os.path.join(work_dir, "card_data.tar.bz2")
        subprocess.run(["curl", "-fsSL", "-o", tarball, url], check=True)
        subprocess.run(["tar", "xjf", os.path.basename(tarball), "./card.json"], cwd=work_dir, check=True)
        os.remove(tarball)

    log_path = os.path.join(work_dir, "card_annotation.log")
    with open(log_path, "wb") as log:
        pixi_run("rgi", ["card_annotation", "-i", "card.json"],
                 cwd=work_dir, stdout=log, stderr=subprocess.STDOUT)

    pattern = re.compile(r"^card_database_v(.+)\.fasta$")
    versions = sorted({
        m.group(1) for fn in os.listdir(work_dir)
        if (m := pattern.match(fn)) and "_all" not in fn
    })
    if len(versions) != 1:
        raise RuntimeError(
            f"expected exactly one card_database_v*.fasta (excluding _all), "
            f"found {len(versions)} -- see {log_path}")
    ver = versions[0]

    pixi_run("rgi", ["clean", "--local"], cwd=work_dir)
    pixi_run("rgi", [
        "load",
        "--card_json", "card.json",
        "--card_annotation", f"card_database_v{ver}.fasta",
        "--card_annotation_all_models", f"card_database_v{ver}_all.fasta",
        "--local",
    ], cwd=work_dir)

    seed_fna = os.path.join(work_dir, "_seed.fna")
    with open(seed_fna, "w") as fh:
        fh.write(">seed\nATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT\n")
    pixi_run("rgi", [
        "main", "-a", "DIAMOND", "-i", "_seed.fna",
        "-o", "_seed_out", "--local", "--clean", "-t", "contig", "-n", "1",
    ], cwd=work_dir)
    for fn in os.listdir(work_dir):
        if fn.startswith("_seed"):
            p = os.path.join(work_dir, fn)
            shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)

    return work_dir


@TaskGenerator
def run_rgi(r1, r2, outdir, rgi_work_dir, threads):
    rgi_outdir = os.path.join(outdir, "rgi")
    os.makedirs(rgi_outdir, exist_ok=True)
    pixi_run("rgi", [
        "bwt",
        "--read_one", r1,
        "--read_two", r2,
        "--output_file", os.path.join(rgi_outdir, "sample.bwt"),
        "--local",
        "--threads", str(threads),
    ], cwd=rgi_work_dir)
    return rgi_outdir


@TaskGenerator
def run_deeparg(r1, r2, outdir, hf_dir):
    deeparg_outdir = os.path.join(outdir, "deeparg")
    os.makedirs(deeparg_outdir, exist_ok=True)

    # deeparg derives intermediate filenames (.paired/.unpaired/.merged/
    # .unmerged) by appending suffixes directly onto whatever path it's
    # given -- symlink inputs into outdir/deeparg/ first so those land there
    # instead of next to the original reads. r1_/r2_ prefixes avoid a
    # collision when both reads files happen to share a basename.
    r1_link = os.path.join(deeparg_outdir, "r1_" + os.path.basename(r1))
    r2_link = os.path.join(deeparg_outdir, "r2_" + os.path.basename(r2))
    try:
        for src, link in ((r1, r1_link), (r2, r2_link)):
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(src, link)

        args = [
            "short_reads_pipeline",
            "--forward_pe_file", r1_link,
            "--reverse_pe_file", r2_link,
            "--output_file", os.path.join(deeparg_outdir, "sample"),
        ]
        if os.path.isdir(hf_dir):
            args += ["--hf-model-path", hf_dir]
        else:
            print(f"warning: DEEPARG_HF_DIR={hf_dir!r} doesn't exist yet -- run "
                  f"'pixi run download-deeparg-db {hf_dir}'; falling back to a "
                  "live Hugging Face download instead.", file=sys.stderr)
        pixi_run("deeparg", args)
    finally:
        for pattern in ("*.paired", "*.unpaired", "*.merged", "*.unmerged"):
            for p in glob.glob(os.path.join(deeparg_outdir, pattern)):
                os.remove(p)
        for link in (r1_link, r2_link):
            if os.path.lexists(link):
                os.remove(link)

    return deeparg_outdir


if len(config.SAMPLES) == 0:
    raise RuntimeError("config.SAMPLES is empty -- add at least one sample")

names = [s["name"] for s in config.SAMPLES]
if len(names) != len(set(names)):
    raise RuntimeError(f"duplicate sample names in config.SAMPLES: {names}")

bad_tools = set(config.TOOLS) - VALID_TOOLS
if bad_tools:
    raise RuntimeError(f"config.TOOLS has unknown entries {bad_tools}, must be a subset of {VALID_TOOLS}")

for sample in config.SAMPLES:
    r1 = os.path.abspath(sample["r1"])
    r2 = os.path.abspath(sample["r2"])
    outdir = os.path.join(os.path.abspath(config.OUTPUT_DIR), sample["name"])

    if "fargene" in config.TOOLS:
        fastq_pair = prepare_fastq_pair(r1, r2, outdir)
        for class_name, model in FARGENE_MODELS.items():
            run_fargene_class(fastq_pair, outdir, class_name, model, config.THREADS)

    if "rgi" in config.TOOLS:
        # Same (work_dir, card_json, card_version) across all samples -- jug
        # de-duplicates identical calls, so this only actually builds once.
        rgi_work_dir = prepare_rgi_db(os.path.abspath(config.RGI_LOCALDB_DIR),
                                      config.CARD_JSON, config.CARD_VERSION)
        run_rgi(r1, r2, outdir, rgi_work_dir, config.THREADS)

    if "deeparg" in config.TOOLS:
        run_deeparg(r1, r2, outdir, os.path.abspath(config.DEEPARG_HF_DIR))
