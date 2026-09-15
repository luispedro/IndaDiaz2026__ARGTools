SAMPLES = [
    {"name": "s1", "r1": "/tmp/jugtest/r1.fastq.gz", "r2": "/tmp/jugtest/r1.fastq.gz"},
    {"name": "s2", "r1": "/tmp/jugtest/r1.fastq.gz", "r2": "/tmp/jugtest/r1.fastq.gz"},
]
OUTPUT_DIR = "/tmp/jugtest/output"
THREADS = 4
TOOLS = ["fargene", "rgi", "deeparg"]
# Kept inside the workspace (and gitignored): `rgi --local` resolves its
# database as localDB/ relative to the cwd it is run from, and pixi needs
# to find pixi.toml by searching up from that same directory.
RGI_LOCALDB_DIR = "rgi_db"
CARD_JSON = None
CARD_VERSION = "4.0.0"
DEEPARG_HF_DIR = "deeparg_hf"
