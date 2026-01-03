# LinearCapR: Linear-time computation of local structural-context probabilities for genome-scale RNA without span constraints

LinearCapR is a **linear-time implementation of CapR** that computes,
for each nucleotide in an RNA sequence without a span limit, the posterior probabilities of
belonging to six local structural contexts:

- Stem (S)
- Hairpin loop (H)
- Bulge loop (B)
- Internal loop (I)
- Multiloop (M)
- Exterior loop (E)

The output format is compatible with CapR: you can drop LinearCapR into
existing pipelines that expect CapR-style structural profiles.

Internally, LinearCapR combines CapR's stochastic context-free grammar (SCFG) with beam search, similar in spirit to LinearFold and LinearPartition. This yields **practical linear-time scaling in sequence length** without imposing a span limit on base pairs, so long-range interactions (e.g. end-to-end pairing in mRNAs or compact viral genomes) are retained.

If you use this software in your research, please cite the LinearCapR
paper (see [Citation](#citation)).

---

## Features

- Linear-time scaling in sequence length (for fixed beam width)
- No explicit span/window limit on base pairs
- Same 6 structural contexts and output format as CapR
- Supports Turner 2004 (default) and Turner 1999 energy parameters
- Optional output of ensemble free energy (`g_ensemble`)
- Suitable for genome- and transcriptome-scale analyses

---

## Build

### Requirements

- C++17 compliant compiler (e.g. `g++` ≥ 9, `clang++` ≥ 10)
- `make`
- Standard C++ STL and libc

### Compilation

```bash
$ make
$ ls
LinCapR  ...
```

## Basic usage
```sh
$ ./LinCapR <input_fasta> <output_file> <beam_size> [options]
```

 - input_fasta FASTA file containing the RNA sequence. The first sequence in the file is used. Sequence characters should be standard RNA bases (A, C, G, U); any non-canonical characters are treated as unpaired.
 - output_file
Path to the output structural profile file. The format is identical to CapR (see Output format).
 - beam_size
Beam width $b$ used for beam pruning in the dynamic programming:
 - $b > 0$: keep at most $b$ highest-scoring states per non-terminal
and position
 - $b = 0$: disable beam pruning (effectively “infinite beam”); this
is only recommended for short sequences because time and memory
scale cubically without pruning.
 - options
	 - `-e`
	   - also print the ensemble free energy g_ensemble to stdout.
	 - `--energy`
	   - choose Energy model: turner2004 (default) or turner1999
	 - `--help`
	   - show this description

### Output format
The output structural profile is line-based and compatible with CapR:
1.	1st line:
Sequence name (the header line in the input FASTA, without >).

2.	Next 6 lines:
For each of the six structural contexts, one line containing N floating-point numbers separated by spaces, where N is the sequence length.

The contexts appear in the following order:
   1. Bulge (B)
   2. Exterior (E)
   3. Hairpin (H)
   4. Internal (I)
   5. Multiloop (M)
   6. Stem (S)

Each number is the posterior probability that the corresponding
nucleotide belongs to that context. For every position i,
```math
p_S(i) + p_H(i) + p_B(i) + p_I(i) + p_M(i) + p_E(i) = 1
```

### Example
For a sequence of length 5, the output looks like:
```txt
example_seq
Bulge 0.80 0.75 0.10 0.05 0.02
Exterior 0.05 0.10 0.60 0.10 0.05
Hairpin 0.05 0.05 0.10 0.60 0.10
Internal 0.05 0.05 0.10 0.10 0.60
Multiloop 0.05 0.05 0.10 0.15 0.23
Stem 0.00 0.00 0.00 0.00 0.00
```
(Values here are illustrative only.)

### Recommended settings

The optimal beam size depends on sequence length and available memory.

Based on the experiments in the LinearCapR paper:
- For short RNAs (≤ ~1 kb)
   - beam_size = 50 is usually sufficient.
   - For maximum accuracy, you can try beam_size = 100.
- For medium RNAs (~1–10 kb)
   - beam_size = 100 is a good default.
   - Increasing the beam beyond 200 yields only marginal AUC improvements on bpRNA datasets.
- For long RNAs / viral genomes (~10–30 kb)
   - beam_size = 100–200 is recommended.
   - On the 29.9 kb SARS-CoV-2 genome, the ensemble free energy and context distributions stabilize for beam_size >= 100, while runtime and memory scale linearly with the beam width.
- For very long RNAs (hundreds of kb)
   - Use beam_size in the range 50–150, depending on your memory
budget.
   - Larger beams improve accuracy but also increase memory usage.

If you are unsure, start with:
```sh
beam_size = 100
```
and adjust upward (for accuracy) or downward (for speed/memory) based on your application.

### Example commands

Compute structural profiles with beam size 100:
```sh
$ ./LinCapR example.fa example.profile 100
```

Compute structural profiles and print the ensemble free energy:
```sh
$ ./LinCapR example.fa example.profile 100 -e
# g_ensemble is printed to stdout
```

Use the output as input to downstream tools (e.g. custom scripts, R/Python analysis, or CapR-compatible pipelines) by parsing the 6 context-probability lines.

## Notes on the algorithm
 - LinearCapR uses a beam-pruned inside–outside dynamic program over
a CapR-like SCFG for pseudoknot-free RNA secondary structures.
 - The grammar is truncated so that consecutive unpaired nucleotides in
loops are capped by a parameter C (default 30), which covers more
than 99% of multiloop unpaired runs in the bpRNA-1m(90) dataset.
 - All dynamic programming values are accumulated in log space with a
numerically stable log-sum-exp approximation, following CONTRAfold
and LinearPartition.

Details and pseudocode are provided in the main text and Supplementary
Information of the LinearCapR paper.

## Citation

If you use LinearCapR in published work, please cite:

```bibtex
@ARTICLE{Otagaki2025-cu,
  title       = "{LinearCapR}: Linear-time computation of per-nucleotide
                 structural-context probabilities of {RNA} without base-pair
                 span limits",
  author      = "Otagaki, Takumi and Hosokawa, Hiroaki and Fukunaga, Tsukasa and
                 Iwakiri, Junichi and Terai, Goro and Asai, Kiyoshi",
  journal     = "bioRxiv",
  institution = "bioRxiv",
  pages       = "2025.12.26.696559",
  abstract    = "Motivation: RNA molecules adopt dynamic ensembles of secondary
                 structures, where the local structural context of each
                 nucleotide--such as whether it resides in a stem or a specific
                 type of loop--strongly shapes molecular interactions and
                 regulatory function. Structural-context probabilities therefore
                 provide a more functionally informative view of RNA folding
                 than the minimum free energy structures or base-pairing
                 probabilities. However, existing tools either require O(N 3 )
                 time or employ span-restricted approximations that omit
                 long-range base-pairs, limiting their applicability to large
                 and biologically important RNAs. Results: We introduce
                 LinearCapR, enabling linear-time, span-unrestricted computation
                 of structural-context marginalized probabilities, using
                 beam-pruned Stochastic Context Free Grammar-based computation.
                 LinearCapR retains global ensemble features lost by
                 span-limited methods and yields superior predictive power on
                 bpRNA-1m(90) dataset, especially for multiloops and exterior
                 regions, as well as long-distance stems. LinearCapR supports
                 analysis of long RNAs, demonstrated on the full genome of
                 SARS-CoV-2. Conclusions: LinearCapR provides the first
                 base-pair-span-unrestricted, linear-time framework for RNA
                 structural context analysis, retaining key thermodynamic
                 ensemble features essential for functional interpretation. It
                 enables large-scale studies of viral genomes, long non-coding
                 RNAs, and downstream analyses such as RNA-binding protein site
                 prediction. Availability: The source code of LinearCapR is
                 available at https://github.com/hoget157/LinearCapR.",
  month       =  dec,
  year        =  2025,
  language    = "en"
}

```

If you want to reproduce all the graphs and tables in the paper, please visit the following repository:
 - https://github.com/TakumiOtagaki/LinearCapR_ComputationalExperiments

## Contact
For questions, bug reports, or feature requests, please open an issue in this repository or contact the corresponding author:
 - Takumi Otagaki: takumiotagaki@gmail.com
