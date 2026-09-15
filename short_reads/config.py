# One entry per sample. Reads are picked up from METAGENOMES_DIR/<name>/ by
# ngless' load_fastq_directory, which accepts .fq/.fastq (optionally .gz/.bz2/
# .xz) and pairs them up on a .1/.2 or _1/_2 suffix, e.g.
#     data/metagenomes/s1/s1.pair.1.fq.gz
#     data/metagenomes/s1/s1.pair.2.fq.gz
SAMPLES = ["s1", "s2"]
METAGENOMES_DIR = "data/metagenomes"

OUTPUT_DIR = "output"

# Parent directory for the per-invocation scratch directories, which hold the
# preprocessed reads plus whatever the tool scribbles next to them. None means
# tempfile's default ($TMPDIR, else /tmp). On a cluster, point this at
# node-local scratch: every byte written here is thrown away once the tool's
# results have been copied into OUTPUT_DIR.
TMP_DIR = None

THREADS = 4
TOOLS = ["fargene", "rgi", "deeparg"]

# `rgi --local` resolves its database as localDB/ relative to the directory it
# is run from; jugfile.py builds it once here and symlinks it into each run's
# scratch directory. Gitignored.
RGI_LOCALDB_DIR = "rgi_db"
CARD_JSON = None
CARD_VERSION = "4.0.0"
DEEPARG_HF_DIR = "deeparg_hf"
