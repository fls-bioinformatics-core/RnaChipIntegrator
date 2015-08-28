Input files
===========

``RnaChipIntegrator`` expects two input files: a list of genomic
features and a list of peaks.

Features data file
******************

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

If ``DE_flag`` is present then it is used to indicate whether the
feature should be considered to be differentially expressed.

Note that any additional columns will be are discarded.

Peaks data file
***************

The peaks data file must be a tab-delimited file with at least 3
columns of data for each peak (one per line)::

    chr  start  end

where:

* ``chr`` is the chromosome that the peak appears on
* ``start`` and ``end`` define the limits of the peak region

Any additional columns will be are discarded.

.. note::

   In previous versions of ``RnaChipIntegrator`` a distinction was
   made between peak 'regions' and peak 'summits', depending on
   whether the ``start`` and ``end`` positions defined a region of
   width 1 (i.e. a summit) or greater than 1 (i.e. a region).

   For this version of the program no distinction is made and the
   same analyses are performed regardless of whether the data
   define summits or regions.
