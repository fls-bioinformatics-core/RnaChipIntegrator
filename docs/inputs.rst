.. _inputs:

Input files
===========

``RnaChipIntegrator`` expects two input files: a list of genomic
features and a list of peaks.

.. _features_data_file:

Features data file
------------------

The features data file must be a tab-delimited file with at least
5 columns of data for each genomic feature (one per line)::

    ID  chr  start  end  strand

where:

* ``chr`` is chromosome the feature appears on
* ``start`` and ``end`` define the limits of the feature
* ``strand`` is the strand direction (either ``+`` or ``-``)
* ``ID`` is a name which is used to identify the genomic feature
  in the output.

Optionally there can be a sixth column::

    ID  chr  start  end  strand  DE_flag

If ``DE_flag`` is present then it can be used to indicate whether
the feature should be considered to be differentially expressed
(``DE_flag`` = 1) or not (``DE_flag`` = 0);
see :ref:`using_differential_expression_data`.

Note that any additional columns will be are discarded.

Peaks data file
---------------

The peaks data file must be a tab-delimited file with at least 3
columns of data for each peak (one per line). By default the
first 3 columns should be::

    chrom  start  end

where:

* ``chrom`` is the chromosome that the peak appears on
* ``start`` and ``end`` define the limits of the peak region

The ``--peak_cols`` option can be used to specify an arbitrary
set of three columns to use for the chromosome and start and end
positions. For example::

    --peak_cols=2,4,5

will use the values from the 2nd, 4th and 5th columns for ``chrom``,
``start`` and ``end`` respectively.

In either case, any additional columns are ignored.

.. note::

   In previous versions of ``RnaChipIntegrator`` a distinction was
   made between peak 'regions' and peak 'summits', depending on
   whether the ``start`` and ``end`` positions defined a region of
   width 1 (i.e. a summit) or greater than 1 (i.e. a region).

   For this version of the program no distinction is made and the
   same analyses are performed regardless of whether the data
   define summits or regions.
