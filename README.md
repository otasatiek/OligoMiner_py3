# OligoMiner (Python 3 ports)

This repository provides **AI-assisted, Python 3–compatible ports** of selected **OligoMiner** scripts (originally developed for Python 2.7).

## Install

### 1) Create a conda environment (Miniforge recommended)

```bash
# (Optional) channel setup (Miniforge users typically already have conda-forge)
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --set channel_priority strict

# Create env
mamba create -n oligominer \
  python=3.11 \
  biopython numpy scipy scikit-learn \
  bowtie2 kmer-jellyfish \
  -y

conda activate oligominer
```

### 2) Quick sanity check

```bash
python -c "import Bio, numpy, scipy, sklearn; print('python deps: OK')"
bowtie2 --version
jellyfish --version
```

---

## Requirements

### Python packages

* biopython
* numpy
* scipy
* scikit-learn (required if you use the LDA-based specificity model)

### External tools

* Bowtie2 (required for genome alignment / off-target screening)
* Jellyfish (required only if you use k-mer abundance screening in your workflow)

### Reference data you must prepare

* Genome FASTA
* Bowtie2 index (recommended: build from **unmasked** genome FASTA)
* (Optional) Jellyfish `.jf` dictionary (recommended: build from **unmasked** genome FASTA)

---

## Usage

> The canonical flow is: **blockParse → Bowtie2 alignment → outputClean → (optional) k-mer filter → probeRC**.

### 0) Prepare Bowtie2 index (once per genome)

```bash
# Example
bowtie2-build genome.fa genome_index_prefix
```

(Recommended: use an unmasked genome FASTA to build the index.)

### 1) Generate candidate probes from a FASTA region (blockParse)

```bash
python blockParse_py3.py -f target.fa
# -> produces target.fastq (candidate probes)
```

(The original guide recommends using **repeat-masked** FASTA as input to probe mining.)

### 2) Align candidates to the genome (Bowtie2)

```bash
# Example command patterns from the original guide
bowtie2 -x /path/to/genome_index_prefix -U target.fastq --no-hd -t -k 100 --very-sensitive-local -S target_u.sam
# or
bowtie2 -x /path/to/genome_index_prefix -U target.fastq --no-hd -t -k 2 --local -D 20 -R 3 -N 1 -L 20 -i C,4 --score-min G,1,4 -S target.sam
```

(Replace `/path/to/genome_index_prefix` with your Bowtie2 index prefix.)

### 3) Specificity filtering (outputClean)

**A. Unique-only mode (exactly 1 alignment, no LDA):**

```bash
python outputClean_py3_ver2.py -u -f target_u.sam
# -> target_u_probes.bed (name depends on script options)
```

(Unique mode exists in the base script.)

**B. LDA mode (thermodynamics-aware model; requires scikit-learn):**

```bash
python outputClean_py3_ver2.py -T 42 -f target.sam
```

(The guide describes the `-T 42` LDA usage and notes sklearn is required.)

### 4) (Optional) Build a Jellyfish dictionary and screen high-abundance k-mers

Create a dictionary (example from the original guide; adjust `-s`, `-m`, and filenames for your genome/compute):

```bash
jellyfish count -s 3300M -m 18 -o genome_18.jf --out-counter-len 1 -L 2 genome.fa
```

If you have a k-mer filtering script in this repo, run it here using your `.jf`.

```bash
python kmerFilter.py -f target_probes.bed -m 18 -j 18 -j genome_18.jf -k 4
```

### 5) (Optional) Reverse-complement probes (probeRC)

```bash
python probeRC_py3.py -f target_probes.bed
# -> target_probes_RC.bed
```

(The original pipeline includes this conversion step.)

---

## Notes: `outputClean_py3_ver2.py` differences vs the base `outputClean_py3.py`

### 1) New “Unique+Zero” mode (`-uz / --unique_zero`)

`outputClean_py3_ver2.py` adds a new CLI option to **return probes with either 0 alignments or exactly 1 alignment**, without using the LDA model.
This is useful when you want to keep:

* **Uniquely-mapped** candidates, and also
* **Unmapped** candidates (e.g., when screening against a genome index that lacks a transgene sequence, or for certain exogenous targets)

The new flag is handled explicitly in the SAM filtering branch.

### 2) Logging / reporting text updated for Unique+Zero

`outputClean_py3_ver2.py` prints and reports messages that explicitly mention “Unique+Zero mode active.”
The base `outputClean_py3.py` prints only “unique”, “zero”, or “LDA model” messages and does not include this combined mode.

---

## Citation

If you use OligoMiner concepts or reproduce the workflow, **please cite the original paper**:

* Beliveau BJ et al. *OligoMiner: A rapid, flexible environment for the design of genome-scale oligonucleotide in situ hybridization probes.* PNAS. doi:10.1073/pnas.1714530115

(Also see the OligoMiner supplemental installation / usage guide bundled with the original distribution.)

---

## License

This repository is released under the **MIT License** (see `LICENSE`).

Portions of this repository are **modified versions** (Python 3–compatible ports and related changes) of scripts originally distributed as part of **OligoMiner**, developed by the **Molecular Systems Lab** (Wyss Institute for Biologically-Inspired Engineering, Harvard University). The upstream MIT license text and copyright notices are preserved in the headers of the modified source files where applicable.

Modified files include a statement such as:
> “This file is a modified version of ‘<original>.py’ originally part of OligoMiner. Modified for Python 3 and modern Biopython compatibility by Keita Sato.”

**Third-party dependencies:** Some workflows may require external tools and Python packages (e.g., Bowtie2, Jellyfish, NumPy, SciPy, scikit-learn, Biopython). These third-party components are **not distributed** with this repository and are governed by their respective licenses.

For additional attribution and notes, see `NOTICE`.

---

## Disclaimer

This is an **unofficial** port. It is provided **as-is**, without warranty. Validate results for your specific genome, alignment settings, and hybridization conditions. (The upstream project also explicitly provides the software without warranty under MIT.)
