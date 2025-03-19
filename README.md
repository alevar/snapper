# SNAPPER: BAM/SAM/CRAM Intron adjustment guided by reference intron positions.

[![PyPI version](https://badge.fury.io/py/snapper-av.svg)](https://pypi.org/project/snapper-av/)
[![GitHub Downloads](https://img.shields.io/github/downloads/alevar/snapper/total.svg)](https://github.com/alevar/SNAPPER/releases/latest)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://opensource.org/licenses/GPL-3.0)

## Introduction

SNAPPER is a method designed to adjust introns in spliced alignments to best match their reference transcriptomic positions.
Aligners such as minimap2 are built to only consider guide annotations on the genome to which sequences are being aligned to.
However, in case of transcript alignments where sequences have often been extracted from the reference genome using an annotation,
there is an additional unused prior in the form of reference exon boundaries in the sequences themselves.
it is only logical to use this prior information to guide the alignment, alas such functionality is not present in most aligners.
To circumvent this issue, SNAPPER was developed to adjust the intron positions in the alignment to best match the reference transcriptomic positions.
The method reconstructs the alignment DP matrix, augments it with weights for reference intron positions and performs a traceback penalizing mismatches
with respect to the original alignment and prioritizing intron positions from the reference annotation. The result is a more consistent alignment
especially when it comes to alignment between poorly conserved genomes.

## Publications

Coming soon...

## Documentation

## Installation

### Via PyPI

The easiest way to install SNAPPER is through PyPI:

```bash
$ pip install snapper-av
$ snapper --help
```

To uninstall SNAPPER:

```bash
$ pip uninstall snapper-av
```

### Building from source

To build from source, clone the git repository:

```bash
$ git clone https://github.com/alevar/snapper.git --recursive
$ cd snapper
$ pip install -r requirements.txt
$ pip install .
```

### Requirements

| Requirement | Details |
| ----------- | ------- |
| Language support | Python ≥ 3.6 |
| Dependencies | - |

## Getting started

### Usage

```bash
snapper [-h] -s SAM -r REFERENCE [-o OUTPUT] [--qry_intron_match_score QRY_INTRON_MATCH_SCORE] 
        [--trg_pos_match_score TRG_POS_MATCH_SCORE] [--trg_pos_mismatch_score TRG_POS_MISMATCH_SCORE]
```

### Options

| Option | Description |
| ------ | ----------- |
| `-s, --sam` | Path to the SAM/BAM alignment file. Read names in the alignment are expected to match corresponding transcript_id in the reference annotationPath to the query GTF/GFF annotation file. |
| `-r, --reference` | Path to the reference annotation. transcript_id field is expected to match read names in the sam alignment. |
| `-o, --output` | Path to the output SAM/BAM file. |
| `--qry_intron_match_score` | Score for matching query introns. |
| `--trg_pos_match_score` | Score for matching target positions. |
| `--trg_pos_mismatch_score` | Score for mismatching target positions. |

### Help Options

| Option | Description |
| ------ | ----------- |
| `-h, --help` | Prints help message. |

## Example Data

Sample datasets are provided in the "example" directory to test and get familiar with SNAPPER.

The included example can be run with the following command from the root directory of the repository:

```bash
snapper --sam example/example.gtf --reference example/example.gtf --output example/output.sam
```