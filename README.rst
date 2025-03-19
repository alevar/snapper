SNAPPER: BAM/SAM/CRAM Intron adjustment guided by reference intron positions.
=======================================================================================================================================================================

.. image:: https://badge.fury.io/py/snapper.svg
    :target: https://pypi.org/project/snapper/
    :alt: SNAPPER Install
.. image:: https://img.shields.io/github/downloads/alevar/snapper/total.svg
    :target: https://github.com/alevar/SNAPPER/releases/latest
    :alt: Github All Releases
.. image:: https://img.shields.io/badge/License-GPLv3-blue.svg
    :target: https://opensource.org/licenses/GPL-3.0
    :alt: GPLv3 License

.. contents::
    :local:
    :depth: 2

Introduction
^^^^^^^^^^^^

SNAPPER is a method designed to adjust introns in spliced alignments to best match their reference transcriptomic positions. 
Aligners such as minimap2 are built to only consider guide annotations on the genome to which sequences are being aligned to.
However, in case of transcript alignments where sequences have often been extracted from the reference genome using an annotation,
there is an additional unitilized prior in the form of reference exon boundaries in the sequences themselves.
it is only logical to use this prior information to guide the alignment, alas such functionality is not present in most aligners.
To circumvent this issue, SNAPPER was developed to adjust the intron positions in the alignment to best match the reference transcriptomic positions.
The method reconstructs the alignment DP matrix, augments it with weights for reference intron positions and performs a traceback penalizing mismatches 
with respect to the original alignment and prioritizing intron positions from the reference annotation. The result is a more consistent alignment 
especially when it comes to alignment between poorly conserved genomes.


Publications
^^^^^^^^^^^^
Coming...

Documentation
^^^^^^^^^^^^

Installation
^^^^^^^^^^^^

PyPi
""""""""""""

By far the easiest way to install SNAPPER is by using PyPi.

::

        $ pip install snapper

Building from source
""""""""""""""""""""

If you want to build it from source, we recommend cloning the git repository as shown below.

::

    $ git clone https://github.com/alevar/snapper.git --recursive
    $ cd snapper
    $ pip install -r requirements.txt
    $ pip install .

.. list-table:: **Requirements**
   :widths: 15 35

   * - Language support
     - Python ≥ 3.6

Getting started
^^^^^^^^^^^^^^^

Usage: vira [-h] -s SAM -r REFERENCE [-o OUTPUT] [--qry_intron_match_score QRY_INTRON_MATCH_SCORE] [--trg_pos_match_score TRG_POS_MATCH_SCORE]
                   [--trg_pos_mismatch_score TRG_POS_MISMATCH_SCORE]

Options:

  -s SAM, --sam SAM     Path to the SAM/BAM alignment file. Read names in the alignment are expected to match corresponding transcript_id in the reference annotation
  -r REFERENCE, --reference REFERENCE
                        Path to the reference annotation. transcript_id field is expected to match read names in the sam alignment.
  -o OUTPUT, --output OUTPUT
                        Path to the output SAM/BAM file
  --qry_intron_match_score QRY_INTRON_MATCH_SCORE
                        Score for matching query introns
  --trg_pos_match_score TRG_POS_MATCH_SCORE
                        Score for matching target positions
  --trg_pos_mismatch_score TRG_POS_MISMATCH_SCORE
                        Score for mismatching target positions

Help options:

  -h, --help            Prints help message.

Data
^^^^

Sample datasets are provided in the "example" directory to test and get familiar with SNAPPER.
The included examples can be run with the following base commands from the root directory of the repository:

1. snapper --sam example/example.gtf --reference example/example.gtf --output example/output.sam
