# nf-core/mag hackathon roadmap: replace bin-level QUAST with BBMap `statswrapper.sh`

_Last revised: 2026-10-15 (Barcelona hackathon, 2 x 6 h, one contributor)._
_Base: `upstream/dev` @ `264a7fe2`. Latest release: v5.5.0. Working branch: `quast-statswrapper`._

This file is both the plan and a learning guide. Every step says **what** to do, **why**, and
**which command to run yourself** (nothing here is run for you). Steps marked 🔎 are things I
could not verify from the code alone: check them before relying on them.

---

## 0. Ground rules

| Rule | Source |
| --- | --- |
| PRs target `dev`, never `main` | `.github/CONTRIBUTING.md` |
| Check for / open an issue first, then fork + branch | CONTRIBUTING |
| Keep PRs small and focused; no unrelated refactors; review everything generated | CONTRIBUTING (AI policy) |
| Tests: `nf-test`; lint: `nf-core pipelines lint`; CI must pass before review | CONTRIBUTING |
| New step: params in `nextflow.config` + schema + docs (`usage.md`, `output.md`, `CITATIONS.md`) + `docs/usage/resource_guidance.md` | CONTRIBUTING |
| Soft budget: **< 100 changed files per PR** (each snapshot/`meta.yml`/test counts) | plan-nf-core-releases skill |
| One theme per PR | same |
| Changing columns of `bin_summary.tsv` or the output structure is **breaking** under nf-core semver | nf-core spec, 🔎 confirm by grepping the specs (see 7.3) |

Tag legend: `#n` = GitHub issue/PR on `nf-core/mag` unless stated.

---

## 1. Overview

| PR | Theme | Type | Est. files | Blocked by |
| --- | --- | --- | --- | --- |
| Issue | Open "Replace QUAST_BINS with BBMap statswrapper" issue, comment on #789 / #748 | n/a | 0 | none |
| **PR 0** | `bbmap/statswrapper` module (local first, nf-core/modules later) | additive | ~4 (local) / ~5 (upstream) | experiment in section 3 |
| **PR 1** | Add statswrapper bin stats to `bin_summary.tsv`, alongside QUAST | additive (minor) | ~15-25 | PR 0 |
| **PR 2** | Remove `QUAST_BINS` (keep or drop assembly-level `QUAST`: open question) | **breaking** | ~20-40 | PR 1 + maintainer decision |
| Stretch | Bakta as optional alternative to Prokka (#346) | additive | ~25-35 | nothing |
| Discuss | Bins entry point (#623) | n/a | 0 | maintainers |

Realistic target for two 6-hour days, working alone and learning: **Issue + PR 0 + PR 1**.
PR 2 only needs doing once the remote maintainer has agreed to the breaking change; it can
be drafted but should not be merged at the hackathon.

---

## 2. Background: what exists today (verified on `dev`)

- `modules/local/quast_bins/main.nf` loops over every bin and runs
  `metaquast.py --max-ref-number 0 --rna-finding --gene-finding`, then concatenates each
  `transposed_report.tsv` into `QUAST/<prefix>-quast_summary.tsv`. Gene finding is the slow
  part (#789).
- `workflows/mag.nf:551-567`: runs `QUAST_BINS`, concatenates the per-group tables with
  `CONCAT_QUAST_SUMMARY` (nf-core `qsv/cat`), gives `ch_quast_bins_summary`.
- `workflows/mag.nf:628-641`: `BIN_SUMMARY` receives `ch_quast_bins_summary` as its second input.
- `bin/combine_tables.py:235-252`: reads the table, strips `.gz` from the `Assembly` column,
  **exits with an error if its bin names do not exactly match the depths table**, suffixes every
  column with `_quast`, merges on `bin == Assembly_quast`.
- `workflows/mag.nf:729-733`: QUAST outputs go to MultiQC.
- `--skip_quast` (`nextflow.config:77`) switches off **both** the assembly-level `QUAST`
  process and `QUAST_BINS`. Also referenced in `subworkflows/local/utils_nfcore_mag_pipeline/main.nf`
  (bigMAG check at ~419, tool list at ~504, citation at ~633).
- Also mentions QUAST: `conf/modules.config` (publishDir), `conf/base.config`
  (`QUAST_BINS` memory), `docs/output.md`, `docs/usage/new_to_mag.md`,
  `docs/usage/resource_guidance.md`, metromap SVGs, test snapshots (`tests/*.nf.test.snap`),
  `ro-crate-metadata.json`.
- Maintainer position in #789 (jfy133): the loss of the "# predicted rRNA genes" column is the
  only worry; Prokka's `.txt` could supply it. Aim: all MIMAG-relevant info in `bin_summary.tsv`.
- Local modules live in `modules/local/<name>/` (`main.nf`, `environment.yml`, tests).
  🔎 Look at one existing local module, e.g. `modules/local/tiara_classify/`, as a template.

What `statswrapper.sh` is: a BBMap script that runs `stats.sh` on many FASTA files and
prints **one row per file** (usage text: "Runs stats.sh on multiple assemblies to produce one
output line per file"). `format=3` is a tab-delimited header + data row. It does no gene
finding, so it is fast, but it also cannot replace QUAST's gene/rRNA counts.

There is **no** `bbmap/statswrapper` in nf-core/modules (checked 2026-10-15 against
`modules/nf-core/bbmap/`; I found no open PR). 🔎 Re-check on the day.

---

## 3. Step 1: learn the tool before writing any Nextflow (about 45 min)

Goal: find out the exact output columns so we know how big PR 1 really is.

### 3.1 Install BBMap locally

Pick one. (Conda/mamba is closest to what the pipeline's `environment.yml` will use.)

```bash
# option A: conda/mamba
mamba create -n bbmap -c conda-forge -c bioconda bbmap
mamba activate bbmap

# option B: macOS homebrew
brew install bbtools    # 🔎 check the formula name with: brew search bbmap
```

Check: `statswrapper.sh` (no arguments) should print the usage text.

### 3.2 Get a few bins to test on

Any FASTA files work. Easiest: take bins from a previous mag run (`results/GenomeBinning/...`),
or published MAGs. Put 3-5 files in a scratch folder, e.g. `~/scratch/bins/`. Include one
gzipped file (`.fa.gz`) because mag passes bins to modules gzipped.

### 3.3 Run it and look at the output

```bash
cd ~/scratch/bins
statswrapper.sh *.fa *.fa.gz format=3 > stats.tsv       # 🔎 does it accept .gz? does it accept a list of files this way?
statswrapper.sh in=a.fa,b.fa format=3                  # the documented comma form
cat stats.tsv
```

Write down answers to these (they decide the design):

1. Is the header printed once or once per file? (If once per file, the module must strip them.)
2. What is the file-name column called, and does it contain the full file name incl. `.gz`?
   (`combine_tables.py` matches bin names exactly.)
3. Which columns appear (N50, total length, GC, contig counts, ...)? Compare with the
   `*_quast` columns in a real `bin_summary.tsv`.
4. Does it work on gzipped input? If not, the module must `zcat` or decompress first.
5. How long does it take vs. QUAST on the same bins? (`time statswrapper.sh ...`)
6. Memory: the script defaults to `-Xmx200m`; bigger bins may need `-Xmx`. 🔎 See the `-Xmx` handling in the script.

Paste the results into the issue (section 4); this is exactly the evidence a reviewer wants.

---

## 4. Step 2: open the issue (before any code)

nf-core wants an issue first. Search for an existing one, then open one titled e.g.
**"Replace QUAST_BINS with BBMap statswrapper.sh"**:

- Link #789 (QUAST_BINS slow / reconfigure BIN_SUMMARY) and #748 (replace local QUAST modules).
- Paste the section 3 findings: columns gained and lost, runtime comparison.
- Ask the questions in section 9 (breaking vs. staged, assembly-level QUAST, rRNA column).
- Mention you will work on it at the Barcelona hackathon, and ping the maintainer who suggested it.

Also post in the `#mag` Slack channel with the issue link.

---

## 5. PR 0: the `bbmap/statswrapper` module

### Decision: local or nf-core/modules first?

| Option | Pro | Con |
| --- | --- | --- |
| **Local module first** (`modules/local/bbmap_statswrapper/`) (default) | No dependency on another repo's review; unblocks PR 1 immediately | Later needs moving to nf-core/modules (a normal nf-core practice) |
| nf-core/modules PR first | Cleanest end state | Review in a second repo may not finish in two days; you need a fork of `nf-core/modules` and its `nf-core modules create/test` workflow |

Default: **local now, upstream later** (listed as open question 3). Use the
`write-nextflow-local-modules` skill when we get here: it covers the modern conventions.

### 5.1 Files (local version, about 4 files)

```
modules/local/bbmap_statswrapper/
  main.nf                 # the process
  environment.yml         # bbmap pinned to an exact version
  tests/main.nf.test      # nf-test for the process
  tests/main.nf.test.snap # generated by nf-test, do not hand-write
```

### 5.2 Process design (sketch, to be checked against section 3 findings)

- **Input:** `tuple val(meta), path(fasta)`: all bins of one group as a list, mirroring
  what `QUAST_BINS` receives now (`workflows/mag.nf:553-555` sorts them first).
- **Command:** one call, `statswrapper.sh ${fasta} format=3 ... > ${prefix}.tsv`,
  not a bash loop. 🔎 `-Xmx` from `task.memory`.
- **Output:** `tuple val(meta), path("*.tsv"), emit: stats`, plus version via the modern
  topic-channel approach used by other local modules in this repo (copy from an existing one).
- **Labels / resources:** `process_single` or `process_low`; mag also documents resources in
  `docs/usage/resource_guidance.md` (required by CONTRIBUTING).
- `ext.prefix` and `ext.args` from `conf/modules.config`, as in other modules.

### 5.3 How you test it

```bash
nf-test generate process modules/local/bbmap_statswrapper/main.nf   # 🔎 check this subcommand exists in your nf-test version
nf-test test modules/local/bbmap_statswrapper/tests/main.nf.test --profile +docker --update-snapshots
```

`--profile +docker` assumes Docker is installed; use `+conda` or `+singularity` otherwise.
Test input: a small FASTA from `tests/` or a file in the nf-core test-datasets. 🔎 Look at what
existing `modules/local/*/tests/` use.

Done when: the test passes, snapshot is committed, `nextflow lint` is clean on the module.

---

## 6. PR 1: add statswrapper stats to `bin_summary.tsv` (additive)

Why additive first: removing QUAST changes outputs for every user (breaking); adding
columns is minor. It also lets the maintainer compare both side by side, and answers #789's
concern about missing columns with data.

### 6.1 Changes

| # | File | Change |
| --- | --- | --- |
| 1 | `workflows/mag.nf` | `include` the module; call it next to the QUAST block (~551); concatenate per-group tables (reuse `CONCAT_*` pattern with `qsv/cat`, or have the module produce one table per run, a simpler design if possible); pass result into `BIN_SUMMARY` |
| 2 | `modules/local/bin_summary/main.nf` | new input `path statswrapper_sum` and `--statswrapper_summary` flag, same pattern as the others (note `ch_dereplication_cluster_tsv` was added as the 9th input recently: mirror how it was done) |
| 3 | `bin/combine_tables.py` | new `--statswrapper_summary` argument; read, normalise bin name (strip `.gz`), validate against `bins`, `add_suffix("_statswrapper")`, merge on `bin` |
| 4 | `nextflow.config` | new parameter. Naming convention for toggles is `skip_*` (existing: `skip_quast`, `skip_gtdbtk`, `skip_binqc`...). Proposal: `skip_bin_stats = false`. 🔎 Check the list of booleans (command in 7.2) |
| 5 | `nextflow_schema.json` | regenerate: `nf-core pipelines schema build` |
| 6 | `conf/modules.config` | publishDir + `ext.args`; `conf/base.config`: resources if needed |
| 7 | `subworkflows/local/utils_nfcore_mag_pipeline/main.nf` | citation text (BBMap) in the tool list/references, as done for QUAST (~504, ~633) |
| 8 | `CITATIONS.md` | BBMap entry |
| 9 | `docs/output.md`, `docs/usage/resource_guidance.md` | document new columns/outputs and resource use |
| 10 | `CHANGELOG.md` | entry under the `dev` heading, with PR number |
| 11 | `tests/*.nf.test.snap` | update snapshots (`--update-snapshots`) and review the diff: only the expected new files/columns should change |
| 12 | `modules.json` | only if an nf-core module is installed (not for local) |

Estimated about 15-25 files, mostly snapshots.

### 6.2 Things that will bite you

- **Bin name matching.** `combine_tables.py` errors if bin names differ (see `:241-244`).
  statswrapper's file-name column may include `.gz` or a path. Normalise it like the QUAST code does.
- **`skip_*` interactions.** The `if` at `workflows/mag.nf:628` decides when `BIN_SUMMARY` runs,
  and `utils_nfcore_mag_pipeline/main.nf` (~419) contains an error check that includes `skip_quast`.
  Decide how the new toggle should interact; do not silently change existing behaviour.
- **Dereplicated vs. all bins.** QUAST currently runs on `ch_input_for_postbinning`. Keep the same input.
- **Gzip:** see section 3.3 question 4.

### 6.3 How you test it

```bash
nf-core pipelines lint                                   # from the repo root; fix everything it flags
nf-test test --tag test --profile +docker --verbose      # CONTRIBUTING's command
nf-test test --tag test --profile +docker --update-snapshots   # only after reading the snapshot diff
nextflow run . -profile test,docker --outdir /tmp/out_mag      # a quick manual run
```

Open `/tmp/out_mag/GenomeBinning/bin_summary.tsv` and confirm the new `*_statswrapper` columns exist.

Done when: tests + lint pass locally, CI is green on the PR, docs and changelog updated,
PR description links the issue and shows before/after columns.

---

## 7. How-to reference

### 7.1 Branching and PR flow
```bash
git fetch upstream
git switch -c bbmap-statswrapper-module upstream/dev    # one branch per PR
# ... edit, commit small ...
git push -u origin bbmap-statswrapper-module
gh pr create -R nf-core/mag --base dev --fill           # or via the GitHub web UI
```
The PR template has a checklist; tick it honestly. Never develop on `dev` or `main`.

### 7.2 Check existing boolean parameter names (for the `skip_*` convention)
```bash
python3 -c "import json;s=json.load(open('nextflow_schema.json'))
[print(k,v.get('default')) for g in s['\$defs'].values() for k,v in g['properties'].items() if v.get('type')=='boolean']"
```

### 7.3 Check the nf-core semver rules at the source (don't trust memory)
```bash
git clone -q --depth 1 --filter=blob:none --sparse https://github.com/nf-core/website.git
cd website && git sparse-checkout set sites/docs/src/content
grep -rn -i "breaking" sites/docs/src/content/docs/specifications | head
```

### 7.4 When something fails
- Lint failure: the message links the relevant docs page; fix, rerun lint.
- Snapshot mismatch: look at the diff first. Update only when the change is intended.
- Process failed: `work/` hash directories in the error message contain `.command.sh`,
  `.command.err`; `cd` there and read them.
- Stuck: ask in Slack `#mag` with the issue link and the exact error.

---

## 8. Day plan

**Day 1**
1. (30 min) Sync done; read this file; skim `quast_bins/main.nf`, `bin_summary/main.nf`, `combine_tables.py`.
2. (60 min) Section 3: install BBMap, run on a few bins, write down answers.
3. (30 min) Section 4: open the issue and Slack post.
4. (150 min) PR 0: local module + nf-test.
5. (30 min) Open PR 0 (or fold into PR 1 if the maintainer prefers one PR).

**Day 2**
1. (180 min) PR 1 steps 1-3 and 4-6, run the test profile, inspect `bin_summary.tsv`.
2. (90 min) Docs, CITATIONS, changelog, schema, snapshots, lint.
3. (60 min) Open PR 1, answer CI, post in Slack. If time is left: draft PR 2 or the Bakta issue triage.

If the module (PR 0) takes longer than planned, **cut PR 1's polish, not its tests.**

---

## 9. Open questions (for the remote maintainer)

1. **Staged or direct?** Add statswrapper next to QUAST first (PR 1, additive), then remove
   `QUAST_BINS` in a later breaking release (PR 2)? Default assumed: staged.
2. **Assembly-level `QUAST`** (`modules/local/quast_run`): replace too (statswrapper works on
   assemblies as well, but loses misassembly-style metrics only QUAST gives, though with
   `--max-ref-number 0` these are limited anyway) or keep it?
   Default assumed: keep, out of scope.
3. **Module home:** local now and nf-core/modules later, or nf-core/modules first?
4. **rRNA column:** is it acceptable to drop `# predicted rRNA genes_quast`, with Prokka's output
   as the future source (#789)?
5. **Toggle:** is `--skip_quast` reused for the new tool, or add a new toggle (`--skip_bin_stats`)?
6. **Column set:** keep all statswrapper columns or a curated subset in `bin_summary.tsv`?

---

## 10. Backlog and other issues

| Item | Notes |
| --- | --- |
| Bakta as optional alternative to Prokka (#346) | Decision: optional, Prokka stays default (Bakta covers bacteria only; `meta.domain` could route bacterial bins). Needs Bakta DB handling. Stretch only. |
| fairy (#796) | Module merged in nf-core/modules#11204; no committed mag contributor. Good first follow-up. |
| myloasm (#1023) | No design discussion yet. |
| binette (#714) | **Claimed by JeanMainguy.** Skipped. |
| Bins entry point (#623) | Discussion, not a PR. `prototaxites` suggests a separate pipeline; `erikrikarddaniel` has a prototype (`--bins_input`); `dialvarezs` open to iterating. |
| Replace local QUAST modules with nf-core ones (#748) | Superseded in part by this work (no nf-core QUAST module is used if QUAST is removed). |
| Hostile / CoverM / dRep etc. | In progress by others (#857, #997); don't touch. |

## 11. Traceability

| Source | Disposition |
| --- | --- |
| Hackathon task: replace QUAST | Issue, PR 0, PR 1, PR 2 |
| Maintainer idea: statswrapper | PR 0, PR 1 |
| #789 QUAST_BINS slow | PR 1 and PR 2 |
| #748 official QUAST module | Superseded (see backlog) |
| #346 Bakta | Stretch |
| #796 fairy, #1023 myloasm | Backlog |
| #714 binette | Skipped (claimed) |
| #623 entry point | Discussion |

## 12. Decisions log

| Date | Decision | By |
| --- | --- | --- |
| 2026-10-15 | Scope is mostly QUAST/statswrapper; solo | user |
| 2026-10-15 | Bakta, if done, is an optional alternative with Prokka as default | user |
| 2026-10-15 | binette skipped for now | user |
| 2026-10-15 | Assistant does not run commands unless asked | user |
| 2026-10-15 | Staged approach (additive PR 1, breaking PR 2) and local-module-first are **assumed defaults**, pending maintainer answers | assistant |
