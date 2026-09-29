# HERD

[![PyPI](https://img.shields.io/pypi/v/immunoHERD.svg)](https://pypi.org/project/immunoHERD/)
[![Python](https://img.shields.io/pypi/pyversions/immunoHERD.svg)](https://pypi.org/project/immunoHERD/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/isi1000/<REPO-NAME>/blob/main/HERD_colab.ipynb)

**H**ost-aware **E**stimation of antigen immunogenicity and **R**anking through
**D**eep-learning.

HERD scores proteins for immunogenicity **in a given host species**. You give it
protein sequences and a host; it returns one score per sequence, so a proteome
can be ranked for that host. Fifteen hosts are supported, six of them validated
as core.

**[Website](https://isi1000.github.io/HERDweb/) · [Run in the browser](https://isi1000-immunoherd.hf.space) · [Run in Colab](https://colab.research.google.com/github/isi1000/<REPO-NAME>/blob/main/HERD_colab.ipynb) · Paper in preparation**

Three ways to run it, in increasing order of size: the **web app** for a handful
of sequences with no setup, **Colab** for a candidate list on a free GPU, and
this **package** for a whole proteome or a reproducible pipeline.

## Installation

```bash
pip install immunoHERD
```

Python 3.11 or newer. The package installs as `immunoHERD`, imports as `herd`,
and the console command is `herd`.

The model weights ship inside the package. The ESM-2 weights (~2.5 GB) are
downloaded automatically the first time an embedding is computed, and cached
afterwards.

## Quick start

```bash
herd predict proteins.fasta --host Bos_taurus --out results.csv
```

or, from Python:

```python
from herd import predict_fasta

df = predict_fasta("proteins.fasta", host="Bos_taurus")
print(df.head())
```

## Usage

### From Python

Three entry points, all returning the same table:

```python
from herd import predict, predict_fasta, predict_text

predict_fasta("proteins.fasta", host="Bos_taurus")           # a FASTA file
predict_text(">p1\nMKVLAAGIVGLNLGG...", host="Canis_sp")     # pasted text
predict([("p1", "MKVLAAGIVGLNLGG...")], host="Sus_scrofa")   # (id, sequence) pairs
```

`predict_text` also accepts a bare sequence with no header.

Optional arguments on all three: `model` (`"integrated"` or `"individual"`),
`threshold`, `device`, `batch_seqs`.

### From the command line

```bash
herd predict proteins.fasta --host Bos_taurus --out results.csv
herd predict proteins.fasta --host Canis_sp --model individual --out dog.csv
herd hosts        # list the available hosts
herd --version
```

`--host` and `--out` are both required: a long run cannot be lost by forgetting
to save it. Other options are `--model`, `--threshold`, `--batch-seqs` and
`--device`.

Progress and warnings go to standard error, so the CSV stays clean when
redirected.

## Output

A `DataFrame` (or CSV) sorted by descending score:

<!-- PENDIENTE (1/3): salida real de una ejecución neutra.
     El CSV que tengo es el panel de antígenos vacunales, y decidiste que ese
     panel no entra en el README. Necesito otra salida real: corre
     `herd predict <cualquier fasta> --host Bos_taurus --out x.csv` y pásame
     las primeras filas. No la invento. -->

| id | length | host | model | score | prediction | warnings |
|---|---|---|---|---|---|---|
| … | … | … | … | … | … | … |

| Column | Meaning |
|---|---|
| `id` | The FASTA identifier, up to the first space. |
| `length` | Length of the sequence as submitted, before any truncation. |
| `host` | The host the score refers to. |
| `model` | `integrated` or `individual`. |
| `score` | Immunogenicity score in `[0, 1]`, rounded to four decimals. |
| `prediction` | `probable immunogen`, `probable non immunogen`, or `ERROR`. |
| `warnings` | Empty when the sequence was used as given. |

Sequences that fail validation are returned last, with `score` set to `NaN` and
`prediction` set to `ERROR`; the reason is in `warnings`. The same column
records when a sequence was *modified* rather than rejected: non-canonical
residues replaced by `X`, or truncation to 1022 residues.

## Interpreting the score

The score is a sigmoid output in `[0, 1]`, read as *how much this protein looks
like the immunogenic proteins seen for this host during training*. The default
decision threshold is 0.5 and can be changed with `--threshold`.

It is a **ranking** score. Its most reliable use is ordering a set of candidate
antigens to decide what to test first, not reading an individual value as a
probability of immunogenicity — nothing here is calibrated for that.

### Performance

Two evaluation schemes were used, and the difference between them matters:

- **Homology-aware nested cross-validation.** Sequences are clustered by
  similarity and a cluster is kept inside a single fold, so no fold is tested on
  a close homologue of something it trained on. Hyperparameters are tuned within
  each outer fold and never see the test fold.
- **Leave-one-pathogen-out (LOPO).** Every source genus is withheld in turn.
  This is the harder and more realistic number: it estimates performance on a
  pathogen the model has never seen.

**The integrated model**, the default, under homology-aware nested
cross-validation:

| Host | Positives | MCC |
|---|---|---|
| Homo sapiens | 3562 | 0.725 |
| Equus caballus | 335 | 0.633 |
| Sus scrofa | 528 | 0.577 |
| Gallus gallus | 460 | 0.544 |
| Bos taurus | 441 | 0.520 |
| Canis sp. | 338 | 0.488 |

Across all fifteen hosts it reaches a pooled out-of-fold AUROC of 0.916
(AUPRC 0.922): 0.902 over the six core hosts, 0.929 over the nine auxiliary
ones.

**The individual models**, under both schemes:

| Host | Positives | AUROC | AUROC (LOPO) | Genera held out | MCC |
|---|---|---|---|---|---|
| Homo sapiens | 3562 | 0.916 | 0.897 | 66 | 0.693 |
| Equus caballus | 335 | 0.896 | 0.855 | 10 | 0.670 |
| Sus scrofa | 528 | 0.847 | 0.806 | 15 | 0.547 |
| Gallus gallus | 460 | 0.830 | 0.744 | 8 | 0.540 |
| Canis sp. | 338 | 0.794 | 0.740 | 9 | 0.441 |
| Bos taurus | 441 | 0.776 | 0.784 | 16 | 0.420 |

Three things are worth reading off those tables.

**Pooling hosts helps.** Training the six core datasets together, alongside nine
further IEDB-only species, raised AUROC in every core host and MCC in five of
the six. The largest gain was for *Bos taurus*, previously the weakest model
(AUROC 0.776 → 0.826, MCC 0.420 → 0.520); even *Homo sapiens*, with by far the
largest dataset, improved (0.916 → 0.927, MCC 0.693 → 0.725). *Equus caballus*
is the exception on MCC, where its own model stays slightly ahead
(0.670 vs 0.633).

**Performance holds on unseen pathogens.** LOPO AUROC stays between 0.74 and
0.90, so most of the signal is not a memorised pathogen signature. The drop is
small for most hosts and absent for *Bos taurus*. *Gallus gallus* falls the
furthest (0.830 → 0.744), which is consistent with its having only eight source
genera to withhold.

**Dataset size is not everything.** Some of the smallest heads are among the
strongest, and the relationship between how much data a host has and how well it
scores is weaker than it looks.

### What HERD does not model

- **Sequence only.** No structure, no post-translational modifications, no
  epitope accessibility, no adjuvant, formulation, dose or route. A protein can
  score high and still fail experimentally for any of these reasons.
- **Antibody-mediated immunogenicity.** The positives come from B-cell assays
  and from immunoproteomic, ELISA and microarray studies in the literature.
  T-cell responses are not modelled.
- **Natural antigen sequence space.** Engineered constructs — fusions,
  chimeras, carrier-linked designs — sit outside what the model has seen, and
  can score low for that reason alone rather than for any biological one.
- **Protein level.** The training set was filtered to sequences of at least
  40 residues. No minimum is enforced at prediction time, but short peptides
  fall outside the training distribution.
- **The first 1022 residues.** Longer sequences are truncated to the ESM-2
  limit and only that prefix is seen. The `warnings` column records it.
- **Nine of the fifteen hosts are auxiliary heads.** They score well in
  cross-validation (AUROC 0.86–0.96), but they were trained on IEDB records
  alone, without the curated literature set behind the six core hosts, and were
  included to broaden the shared antigen repertoire rather than as validated
  predictors in their own right. Treat their scores as exploratory.

## Models and hosts

Two models, chosen with the `model` argument:

- **`integrated`** (default, recommended) — a pan-species model: one shared
  ESM-2 backbone with fifteen species-specific output heads.
- **`individual`** — one model trained per species. Available only for the six
  core hosts.

| | Hosts |
|---|---|
| **Core** (both models) | `Homo_sapiens`, `Bos_taurus`, `Sus_scrofa`, `Gallus_gallus`, `Canis_sp`, `Equus_caballus` |
| **Auxiliary** (integrated only) | `Camelidae`, `Capra_hircus`, `Cavia_porcellus`, `Macaca_sp_`, `Mus_musculus`, `Non_human_primate`, `Oryctolagus_cuniculus`, `Ovis_aries`, `Rattus_sp_` |

Host names accept common aliases (`Canis_familiaris`, `Homo`, `Rattus`…) and a
loose genus match. `herd hosts` prints the list.

## How it works

1. **Cleaning.** Any character outside `ACDEFGHIKLMNPQRSTVWY` becomes `X`, and
   the sequence is truncated to 1022 residues. Both events are reported in
   `warnings`.
2. **Embedding.** ESM-2 (`esm2_t33_650M_UR50D`, layer 33) through `fair-esm`,
   mean-pooled over the real residues only: a 1280-dimensional vector.
3. **Host-aware head** on top of that embedding. Normalization lives inside the
   `.keras` file, so the raw embedding is fed in unscaled.
4. **Threshold** at 0.5 by default.

The same cleaning and pooling recipe is used at training and at prediction
time; that is what makes the shipped weights valid.

## Training data

The models were trained on balanced datasets of immunogens and matched
non-immunogens, built from two sources:

- a corpus of **experimentally supported immunogens curated by hand from the
  primary literature**, favouring proteome-wide screens where available;
- **IEDB** records, restricted to B-cell assays, which complement the curated
  corpus — notably in viral antigens, which gel- and array-based proteomics
  tend to under-sample.

Each positive is paired with a length-matched negative drawn from the proteome
of its own source organism, and required to show no detectable homology to any
positive in the dataset.

The curated dataset will be deposited on publication. Meanwhile, this
repository ships the **recipe** rather than the data: the numbered pipeline in
`development/dataset-construction/` rebuilds the set from IEDB and UniProt, and
documents the curation decisions behind it. Anyone using data obtained that way
must cite IEDB; see [their citation policy](https://www.iedb.org/citation_v3.php).

## Advanced

### Where the weights are loaded from

By default, the `models/` directory inside the installed package. Set
`HERD_MODELS_DIR` to override it:

```bash
export HERD_MODELS_DIR=/path/to/my/models
```

The expected layout is `<HERD_MODELS_DIR>/integrated.keras` and
`<HERD_MODELS_DIR>/individual/<host>.keras`.

### Speed

Sequences are batched **sorted by length**, which removes almost all padding:
on a real proteome of 4,570 sequences this halves the GPU work, with identical
results — each mean is taken over the real residues of its own sequence.

A GPU is not required, but embedding on CPU is considerably slower. `--device`
forces `cpu` or `cuda`; by default it is chosen automatically.

## Repository contents

Beyond the installable package:

- `HERD_colab.ipynb` — a ready-to-run Colab notebook, for when there is no
  local GPU.
- `development/dataset-construction/` — the numbered pipeline that rebuilds the
  training set from IEDB and UniProt.
- `development/model-training/` — the notebooks that trained the models whose
  weights ship in the package.
- `tests/` — the test suite.

See [`development/README.md`](development/README.md) for how the dataset was
built and how to reconstruct it.

## Tests

```bash
pip install -e .
pytest
```

The `src/` layout means the tests always import the *installed* package, never
the repository folder.

## License

MIT, see [`LICENSE`](LICENSE). The model weights are released under the same
terms.

ESM-2 is used through `fair-esm` and carries its own licence. Data obtained by
running the pipeline in `development/` comes from IEDB and UniProt and must be
cited accordingly.

## Citation

A paper describing HERD is in preparation; the reference will be added here on
publication. In the meantime, please cite this repository, along with **ESM-2**
for the protein embeddings and **IEDB** for the epitope data underlying the
training set.
