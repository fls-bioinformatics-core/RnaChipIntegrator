RnaChipIntegrator: analysis of genes with peak data
===================================================

.. image:: https://readthedocs.org/projects/rnachipintegrator/badge/?version=latest
   :target: https://rnachipintegrator.readthedocs.io

.. image:: https://badge.fury.io/py/RnaChipIntegrator.svg
   :target: https://pypi.python.org/pypi/rnachipintegrator/

.. image:: https://travis-ci.org/fls-bioinformatics-core/RnaChipIntegrator.png?branch=master
   :target: https://travis-ci.org/fls-bioinformatics-core/RnaChipIntegrator

``RnaChipIntegrator`` is a utility that performs integrated analyses
of 'gene' data (a set of genes or other genomic features) with 'peak'
data (a set of regions, for example ChIP peaks) to identify the genes
nearest to each peak, and vice versa.

The program was originally written specifically for ChIP-Seq and RNA-Seq
data but works equally well for ChIP-chip and microarray expression data,
and can also be used to integrate any set of genomic features (e.g.
canonical genes, CpG islands) with peak data.

Quick Start
***********

Install the latest version of the program from the Python Package Index
(PyPI)::

    pip install RnaChipIntegrator

The simplest use of the program is::

    RnaChipIntegrator GENES PEAKS

where ``GENES`` and ``PEAKS`` are tab-delimited files containing the
'gene' and 'peak' data respectively.

This will output two files with the nearest genes for each peak
("peak-centric" analysis), and the nearest peaks for each gene
("gene-centric" analysis).

Full documentation can be found at:

 * http://rnachipintegrator.readthedocs.io/en/latest/

Licensing
*********

This software is licensed under the Artistic License 2.0.
