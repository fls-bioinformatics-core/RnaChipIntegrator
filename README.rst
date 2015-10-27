RnaChipIntegrator: analysis of genomic features with peak data
==============================================================

``RnaChipIntegrator`` is a utility that performs integrated analyses
of feature data (any set of genomic features, for example gene expression
data or canonical gene lists) with 'peak' data (a set of regions, for
example ChIP peaks) to identify the features to each peak, and vice versa.

The program was originally written specifically for ChIP-Seq and RNA-Seq
data but works equally well for ChIP-chip and microarray expression data,
and can also be used to integrate any set of genomic features (e.g.
canonical genes, CpG islands) with peak data.

Basic usage
***********

The simplest use of the program is::

    RnaChipIntegrator FEATURES PEAKS

where ``FEATURES`` and ``PEAKS`` are tab-delimited files containing
the genomic feature and peak data respecitively.

This will output two files with the nearest features for each peak,
and the nearest peaks for each feature.

Use::

    RnaChipIntegrator -h

to see the full set of options for controlling the analyses and the
outputs.

Installation
************

It is recommended to use::

    pip install .

from within the top-level source directory to install the package.

To use the package without installing it first you will need to add the
directory to your ``PYTHONPATH`` environment.

To install directly from github using ``pip``::

    pip install git+https://github.com/fls-bioinformatics-core/RnaChipIntegrator.git

The program depends on the Python ``xlwt``, ``xlrd`` and ``xlutils``
libraries, which should be installed automatically if using ``pip``.

Documentation
*************

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
