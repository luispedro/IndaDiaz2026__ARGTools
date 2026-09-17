"""
Runs fargene/rgi/deeparg over the samples listed in the settings below.

    pixi run jug-status
    pixi run jug-execute

Safe to re-run: jug skips tasks it already finished, even if you add new
samples or edit shared settings in between runs (only the affected tasks'
hashes change).

Each tool lives in its own pixi environment (see pixi.toml) with its own
binary -- this process (the 'jug' environment) only needs jug itself and
shells out via `pixi run -e <env> ...` for the actual work.

The one-time database setup is part of the graph too: prepare_rgi_db builds
CARD's localDB and prepare_deeparg_db downloads the DeepARG model bundle,
each once, with the run_* tasks that need them waiting on the result.

Every tool invocation gets a scratch directory of its own (under
TMPDIR), laid out as

    <scratch>/reads/preproc.pair.{1,2}.fq[.gz]   quality-trimmed reads
    <scratch>/out/                               what the tool is told to write
    <scratch>/localDB -> rgi_db/localDB          (rgi only)

The reads are produced there by running preprocess.ngl through ngless, the
tool then runs with the scratch directory as its working directory, and
finally <scratch>/out is copied to OUTPUT_DIR/<sample>/<tool>/ and the whole
scratch directory is deleted. Everything a tool leaves lying around outside
that out/ directory -- deeparg's trimmomatic/vsearch intermediates, rgi's
bowtie2 scratch -- therefore dies with the scratch directory instead of
accumulating in the output tree, and OUTPUT_DIR/<sample>/<tool>/ holds
nothing but that tool's results (plus the ngless QC report for the reads it
was given).
"""
import contextlib
import glob
import os
import shutil
import re
import subprocess
import tempfile

from jug import TaskGenerator

# Old fashioned NGLess-style sample-list+directories organisation
SAMPLES = [line.strip() for line in open("data/samples.txt")]
METAGENOMES_DIR = "data/metagenomes"

OUTPUT_DIR = "output"

# Parent directory for the per-invocation scratch directories, which hold the
# preprocessed reads plus whatever the tool scribbles next to them. None means
# tempfile's default ($TMPDIR, else /tmp). On a cluster, point this at
# node-local scratch: every byte written here is thrown away once the tool's
# results have been copied into OUTPUT_DIR.
TMPDIR = None

THREADS = 4

# `rgi --local` resolves its database as localDB/ relative to the directory it
# is run from; jugfile.py builds it once here and symlinks it into each run's
# scratch directory. Gitignored.
RGI_LOCALDB_DIR = "rgi_db"
CARD_JSON = None
CARD_VERSION = "4.0.0"
DEEPARG_HF_DIR = "deeparg_hf"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PIXI_MANIFEST = os.path.join(SCRIPT_DIR, "pixi.toml")
PREPROCESS_NGL = os.path.join(SCRIPT_DIR, "preprocess.ngl")
DOWNLOAD_DEEPARG_DB = os.path.join(SCRIPT_DIR, "download_deeparg_db.py")


# Only for the stamp-file name; the script's huggingface_hub import lives
# inside its own function, so importing it here (under 'jug') is fine. It is
# *run* below via `pixi run -e deeparg`, never called from this process.
from download_deeparg_db import STAMP_FILE as DEEPARG_DB_STAMP

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

# Basename handed to preprocess.ngl; ngless turns it into <base>.pair.1.<ext>,
# <base>.pair.2.<ext> and (only for samples that carry unpaired reads)
# <base>.singles.<ext>.
PREPROC_BASE = "preproc"


def pixi_run(env, args, **kwargs):
    """Run `args` inside the `env` pixi environment, installing it if needed.

    This process runs under the 'jug' environment, so it can't see the tools;
    `pixi run -e <env>` activates theirs, which also lets fargene/rgi/deeparg
    find the sibling binaries they shell out to (deeparg -> trimmomatic, rgi
    -> bamtools). --manifest-path is what makes `cwd` free to be anything:
    without it pixi locates pixi.toml by searching upwards from the working
    directory, which would rule out the scratch directories below (they live
    outside the workspace, typically on node-local scratch). The commands are
    passed as plain commands rather than as pixi tasks for the same reason:
    tasks always run from the workspace root.
    """
    return subprocess.run(["pixi", "run", "--manifest-path", PIXI_MANIFEST, "-e", env, *args],
                          check=True, **kwargs)


def publish_results(src, dest):
    """Copy the scratch directory's out/ tree to its final home in OUTPUT_DIR.

    Staged next to `dest` and renamed into place so that an interrupted copy
    never leaves a half-written result looking like a finished one, and so
    that re-running an invalidated task replaces the old results wholesale
    rather than merging into them.
    """
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    staging = f"{dest}.incoming.{os.getpid()}"
    superseded = f"{dest}.superseded.{os.getpid()}"
    shutil.rmtree(staging, ignore_errors=True)
    shutil.copytree(src, staging)
    if os.path.exists(dest):
        os.rename(dest, superseded)
    os.rename(staging, dest)
    shutil.rmtree(superseded, ignore_errors=True)
    return dest


@contextlib.contextmanager
def tool_scratch(sample_dir, tool, compressed):
    """Preprocess one sample into a fresh scratch directory, for one tool run.

    Yields (scratch, out_dir, r1, r2): the directory the tool should run in,
    the directory it should be told to write to (copied out by
    publish_results), and the two mate files. `compressed` picks the output
    extension -- fargene cannot read gzipped FastQ, so it asks for plain .fq
    and lets ngless decompress on the way out.

    The whole tree is removed on the way out, including on failure: jug runs
    without --keep-failed, so a failed task is retried from scratch anyway and
    there is nothing here worth keeping. Cleanup errors are deliberately not
    suppressed: failing to remove the tree usually means something else went
    wrong (a tool still holding files open, a read-only file it left behind).
    """
    if not os.path.isdir(sample_dir):
        raise RuntimeError(f"no read directory for this sample: {sample_dir}")
    pixi_run('fargene', [
        'trim_galore', '--version'])
    tmp_parent = os.path.abspath(TMPDIR) if TMPDIR else None
    if tmp_parent:
        os.makedirs(tmp_parent, exist_ok=True)
    with tempfile.TemporaryDirectory(
            prefix=f"{tool}.{os.path.basename(sample_dir)}.",
            dir=tmp_parent) as scratch:
        reads_dir = os.path.join(scratch, "reads")
        out_dir = os.path.join(scratch, "out")
        os.makedirs(reads_dir)
        os.makedirs(out_dir)

        ext = ".fq.gz" if compressed else ".fq"
        pixi_run("ngless", [
            "ngless",
            "--jobs", str(THREADS),
            # Keep ngless' own scratch inside ours, and its QC report inside
            # out/ so it gets copied out next to the tool's results (the
            # default would be a single preprocess.ngl.output_ngless/ next to
            # the script, which concurrent invocations would fight over).
            "--temporary-directory", scratch,
            "-o", os.path.join(out_dir, "ngless-report"),
            PREPROCESS_NGL,
            sample_dir,
            os.path.join(reads_dir, PREPROC_BASE + ext),
        ], cwd=scratch)

        r1 = os.path.join(reads_dir, f"{PREPROC_BASE}.pair.1{ext}")
        r2 = os.path.join(reads_dir, f"{PREPROC_BASE}.pair.2{ext}")
        missing = [p for p in (r1, r2) if not os.path.isfile(p)]
        if missing:
            # ngless omits the paired outputs entirely for a single-end
            # sample, which none of these three tools are wired up for here.
            raise RuntimeError(
                f"preprocessing {sample_dir} did not produce paired reads "
                f"(missing {', '.join(os.path.basename(p) for p in missing)}); "
                f"ngless wrote {sorted(os.listdir(reads_dir))}")

        yield scratch, out_dir, r1, r2


@TaskGenerator
def run_fargene(sample_dir, outdir):
    with tool_scratch(sample_dir, "fargene", compressed=False) as (scratch, out_dir, r1, r2):
        for class_name, model in FARGENE_MODELS.items():
            pixi_run("fargene", [
                "fargene",
                "-i", r1, r2,
                "--hmm-model", model,
                "--meta",
                "-o", os.path.join(out_dir, class_name),
                "-p", str(THREADS),
                "--force",
            ], cwd=scratch)
        # fargene appends to a fargene_analysis.log in its working directory
        # (not in -o); keep it rather than let it go with the scratch tree.
        log = os.path.join(scratch, "fargene_analysis.log")
        if os.path.isfile(log):
            shutil.move(log, os.path.join(out_dir, "fargene_analysis.log"))
        return publish_results(out_dir, os.path.join(outdir, "fargene"))


@TaskGenerator
def prepare_rgi_db(work_dir, card_json, card_version):
    # Returns work_dir itself (not work_dir/localDB) -- run_rgi symlinks
    # work_dir/localDB into its own scratch directory, and the build steps
    # below likewise run one level above localDB/, which is where `rgi
    # ... --local` looks for it.
    os.makedirs(work_dir, exist_ok=True)
    localdb = os.path.join(work_dir, "localDB")
    if os.path.isdir(localdb):
        # Already built -- nothing to do.
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
        pixi_run("rgi", ["rgi", "card_annotation", "-i", "card.json"],
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

    pixi_run("rgi", ["rgi", "clean", "--local"], cwd=work_dir)
    pixi_run("rgi", [
        "rgi", "load",
        "--card_json", "card.json",
        "--card_annotation", f"card_database_v{ver}.fasta",
        "--card_annotation_all_models", f"card_database_v{ver}_all.fasta",
        "--local",
    ], cwd=work_dir)

    seed_fna = os.path.join(work_dir, "_seed.fna")
    with open(seed_fna, "w") as fh:
        fh.write(">seed\nATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGT\n")
    pixi_run("rgi", [
        "rgi", "main", "-a", "DIAMOND", "-i", "_seed.fna",
        "-o", "_seed_out", "--local", "--clean", "-t", "contig", "-n", "1",
    ], cwd=work_dir)
    for fn in os.listdir(work_dir):
        if fn.startswith("_seed"):
            p = os.path.join(work_dir, fn)
            shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)

    return work_dir


@TaskGenerator
def run_rgi(sample_dir, outdir, rgi_work_dir):
    with tool_scratch(sample_dir, "rgi", compressed=True) as (scratch, out_dir, r1, r2):
        # `rgi ... --local` resolves its database as localDB/ under the
        # working directory, so point one there. A symlink (rather than a
        # copy) keeps this cheap and lets concurrent samples share the build.
        os.symlink(os.path.join(rgi_work_dir, "localDB"),
                   os.path.join(scratch, "localDB"))
        pixi_run("rgi", [
            "rgi", "bwt",
            "--read_one", r1,
            "--read_two", r2,
            "--output_file", os.path.join(out_dir, "sample.bwt"),
            "--local",
            "--threads", str(THREADS),
        ], cwd=scratch)
        return publish_results(out_dir, os.path.join(outdir, "rgi"))


@TaskGenerator
def prepare_deeparg_db(hf_dir):
    """Fetch the DeepARG model+database bundle into `hf_dir`; returns it.

    The download is a separate script run under the 'deeparg' environment
    because it needs huggingface_hub, which this process (the 'jug'
    environment) doesn't have -- see download_deeparg_db.py, which `pixi run
    download-deeparg-db <dir>` runs by hand in exactly the same way.

    A directory carrying the stamp file is taken as already downloaded and
    left alone, so a bundle staged by hand (or by an earlier, since
    invalidated, run of this task) is not re-fetched.
    """
    if os.path.isfile(os.path.join(hf_dir, DEEPARG_DB_STAMP)):
        return hf_dir

    pixi_run("deeparg", ["python", DOWNLOAD_DEEPARG_DB, hf_dir])
    if not os.path.isfile(os.path.join(hf_dir, DEEPARG_DB_STAMP)):
        raise RuntimeError(
            f"{DOWNLOAD_DEEPARG_DB} exited successfully but left no "
            f"{DEEPARG_DB_STAMP} in {hf_dir}")
    return hf_dir


@TaskGenerator
def run_deeparg(sample_dir, outdir, hf_dir):
    with tool_scratch(sample_dir, "deeparg", compressed=True) as (scratch, out_dir, r1, r2):
        pixi_run("deeparg", [
            "deeparg", "short_reads_pipeline",
            "--forward_pe_file", r1,
            "--reverse_pe_file", r2,
            "--output_file", os.path.join(out_dir, "sample"),
            # Without this deeparg downloads the bundle itself, once per run;
            # prepare_deeparg_db has already put it here.
            "--hf-model-path", hf_dir,
        ], cwd=scratch)

        # deeparg builds its trimmomatic/vsearch intermediates by appending
        # suffixes to the paths it is handed. The ones derived from the reads
        # land in reads/ and go with the scratch tree; these are the ones that
        # would otherwise be published alongside the actual results.
        for pattern in ("*.paired", "*.unpaired", "*.merged", "*.unmerged"):
            for p in glob.glob(os.path.join(out_dir, pattern)):
                os.remove(p)
        return publish_results(out_dir, os.path.join(outdir, "deeparg"))


# Prepared once, ahead of the samples; every run_rgi/run_deeparg task takes
# the corresponding directory as an argument and so waits for it.
rgi_work_dir = prepare_rgi_db(os.path.abspath(RGI_LOCALDB_DIR),
                              CARD_JSON, CARD_VERSION)
deeparg_hf_dir = prepare_deeparg_db(os.path.abspath(DEEPARG_HF_DIR))

metagenomes_dir = os.path.abspath(METAGENOMES_DIR)
output_dir = os.path.abspath(OUTPUT_DIR)

# Adding a fourth tool is a run_<tool> task above plus one line here: jug
# schedules the new tasks and leaves every result already on disk untouched.
for sample in SAMPLES:
    # Absolute, and passed as task arguments rather than looked up inside the
    # tasks, so that moving METAGENOMES_DIR invalidates the affected tasks.
    sample_dir = os.path.join(metagenomes_dir, sample)
    outdir = os.path.join(output_dir, sample)

    run_fargene(sample_dir, outdir)
    run_rgi(sample_dir, outdir, rgi_work_dir)
    run_deeparg(sample_dir, outdir, deeparg_hf_dir)
