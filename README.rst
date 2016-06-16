RnaChipIntegrator: analysis of genes with peak data
===================================================

.. image:: https://readthedocs.org/projects/rnachipintegrator/badge/?version=latest
   :target: https://rnachipintegrator.readthedocs.io

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

Full documentation can be found at ReadTheDocs:

 * http://rnachipintegrator.readthedocs.io/en/latest/

See the ``INSTALL`` file for complete installation instructions.

Developers
**********

The source code for the development version of the program is hosted
on GitHub in the ``devel`` branch:

 * https://github.com/fls-bioinformatics-core/RnaChipIntegrator/tree/devel

and can be installed directly from GitHub using ``pip``::

    pip install git+https://github.com/fls-bioinformatics-core/RnaChipIntegrator.git@devel

The program depends on the Python ``xlwt``, ``xlrd`` and ``xlutils``
libraries, which should be installed automatically if using ``pip``.

Documentation based on ``sphinx`` is available under the ``docs`` directory.

To build do either::

    python setup.py sphinx_build

or::

    cd docs
    make html

both of which create the documentation in the ``docs/_build`` subdirectory.

Running Tests
*************

The Python unit tests can be run using::

    python setup.py test

Note that this requires the ``nose`` package.

Examples
********

Example data files can be found in the ``examples`` subdirectory, which
can be used as input to the program for test or demonstration purposes; see
the ``README`` file in the same directory for more information.

Licensing
*********

This software is licensed under the Artistic License 2.0; see the ``LICENSE``
document.
