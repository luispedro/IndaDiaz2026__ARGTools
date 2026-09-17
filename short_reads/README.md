# Short-read ARG detection

Runs three ARG detection tools directly on metagenomic short reads (rather
than on the GMGC10 gene catalog, as in the rest of this repository):

| Tool    | Version                                  | Invocation                      |
|---------|------------------------------------------|---------------------------------|
| fARGene | [indajuan/fargene@3be196d][fargene]      | `fargene --meta`, one run per HMM model (22 models) |
| RGI     | 6.0.3, CARD 4.0.0                        | `rgi bwt --local`               |
| DeepARG | [deeparg-modern v2.0.0][deeparg] (DeepARG v2 bundle) | `deeparg short_reads_pipeline` |

[fargene]: https://github.com/indajuan/fargene
[deeparg]: https://github.com/gaarangoa/deeparg

The pipeline is a [jug](https://jug.readthedocs.io/) script and every tool
lives in its own [pixi](https://pixi.sh) environment.

## Files

- `jugfile.py` — the pipeline. All settings (input/output directories,
  threads, scratch location, CARD version, fARGene models) are at the top.
- `preprocess.ngl` — [NGLess](https://ngless.readthedocs.io/) script that
  quality-trims a sample's reads (`substrim`, min. quality 25, discard reads
  shorter than 45 bp; singletons are dropped).
- `download_deeparg_db.py` — downloads the DeepARG model/database bundle from
  Hugging Face. Run automatically by the pipeline; can also be run by hand.
- `pixi.toml` / `pixi.lock` — environments `fargene`, `rgi`, `deeparg`,
  `ngless`, and `jug` (the last one runs the pipeline itself).

## Input

Samples use the classic NGLess layout, relative to this directory:

```
data/samples.txt                       # one sample name per line
data/metagenomes/<sample>/*.1.fq.gz    # paired-end reads, as understood by
data/metagenomes/<sample>/*.2.fq.gz    # ngless' load_fastq_directory
```

Only paired-end samples are supported.

## Running

```bash
pixi install
pixi run jug-status     # what is done / pending
pixi run jug-execute    # run everything (can be started several times in parallel)
```

Pixi tasks always run from this directory, which matters because
`jugfile.py` resolves `data/samples.txt` (and the other settings' paths)
relative to the working directory. If you call `jug` directly, `cd` here first.

The first run also builds the prerequisites, once:

- **CARD localDB** in `rgi_db/` (downloaded from card.mcmaster.ca, or built
  from a local `card.json` if `CARD_JSON` is set);
- **DeepARG bundle** in `deeparg_hf/`. To stage it ahead of time (e.g. on a
  login node with internet access), run
  `pixi run download-deeparg-db deeparg_hf`.

`jug-execute` uses `--keep-going`, so a failing sample does not stop the
others; failed tasks are retried on the next run.

## How each task runs

For every (sample, tool) pair, `jugfile.py`:

1. creates a scratch directory under `TMPDIR` (default: the system temporary
   directory; point it at node-local scratch on a cluster);
2. runs `preprocess.ngl` to write the trimmed reads there;
3. runs the tool with the scratch directory as working directory;
4. copies the tool's `out/` directory atomically to `output/<sample>/<tool>/`
   and deletes the scratch directory.

Thus, reads are preprocessed separately for each tool (fARGene needs
uncompressed FastQ, the others get gzipped files), and tool intermediates
never reach the output tree.

## Output

```
output/<sample>/fargene/<model>/          # one fARGene output directory per model
output/<sample>/fargene/fargene_analysis.log
output/<sample>/rgi/sample.bwt.*
output/<sample>/deeparg/sample.*
output/<sample>/<tool>/ngless-report/     # QC report of the preprocessing
```

`output/`, `rgi_db/`, `deeparg_hf/`, and `jugfile.jugdata/` are gitignored.

## Caveats

- Only settings passed to the tasks as arguments (directories, `CARD_JSON`,
  `CARD_VERSION`) are part of jug's task hashes. After changing `THREADS`,
  `TMPDIR`, or `FARGENE_MODELS`, results are **not** recomputed automatically;
  use `jug invalidate` if needed.
- `prepare_rgi_db` and `prepare_deeparg_db` skip the work if their target
  directory already looks complete (`rgi_db/localDB/` exists, or
  `deeparg_hf/.download-complete` exists). Delete the directory to force a
  rebuild (e.g. after changing `CARD_VERSION`).
