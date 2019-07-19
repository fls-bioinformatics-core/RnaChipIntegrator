Version History and Changes
===========================

--------------------------
Version 1.2.0 (2019-07-19)
--------------------------

* Implement batch mode operation which allows a single run of
  ``RnaChipIntegrator`` to iterate over any combination of
  multiple cutoffs (via new ``--cutoffs`` option) and multiple
  peaks and genes files (via new ``--peaks`` and ``--genes``
  options)
  (`PR #58 <https://github.com/fls-bioinformatics-core/RnaChipIntegrator/pull/58>`_).
* Drop support for Python 2.6
  (`PR #58 <https://github.com/fls-bioinformatics-core/RnaChipIntegrator/pull/66>`_).

--------------------------
Version 1.1.0 (2018-08-06)
--------------------------

* Fix bug in gene-centric analysis when using ``--edge=both``,
  which meant valid peaks could sometimes be dropped from the
  results
  (`PR #60 <https://github.com/fls-bioinformatics-core/RnaChipIntegrator/pull/60>`_).
* New ``--peak_id`` option enables peak names to be read from
  an ``ID`` column in the peaks file, which will be carried
  through to the output files
  (`PR #56 <https://github.com/fls-bioinformatics-core/RnaChipIntegrator/pull/56>`_).
* New ``--analyses`` option enables specific analyses to be
  selected (i.e. either peak-centric only, gene-centric only,
  or both peak- and gene-centric)
  (`PR #46 <https://github.com/fls-bioinformatics-core/RnaChipIntegrator/pull/46>`_).

--------------------------
Version 1.0.3 (2018-01-23)
--------------------------

* Update description of 'direction' in the documentation and
  the XLSX file, to clarify the meaning of the terms 'upstream'
  and downstream' in the program outputs.
* Update documentation to explain how comments, header lines,
  and lines with 'bad' data, are handled by the program when
  reading in the input files.
* Add a section on 'known problems' to the documentation.

--------------------------
Version 1.0.2 (2016-08-17)
--------------------------

* Bugfix release: fixes bug which crashed the program if
  ``--cutoff=0`` was specified together with ``--xlsx``.

--------------------------
Version 1.0.1 (2016-06-16)
--------------------------

* Bugfix release: fixes bug which crashed the program if
  ``--cutoff=0`` was specified.

--------------------------
Version 1.0.0 (2016-02-24)
--------------------------

* Complete reimplementation of RnaChipIntegrator to unify internal
  algorithms, simplify usage and substantially update the
  documentation.
* Installs using ``pip``; can be obtained directly from PyPI using
  ``pip install RnaChipIntegrator``.
* Executable now installs as ``RnaChipIntegrator`` (i.e. no ``.py``
  extension)
* Documentation now available via ReadTheDocs:
  http://rnachipintegrator.readthedocs.org/en/latest/
* No distinction is now made between 'summits' and 'peaks'; the
  same algorithm is applied in each case.
* The program always finds the nearest genes to each peak, and
  vice versa. The same distance cutoff and maximum number of hits
  are applied to both and can be specified using the ``--cutoff``
  and ``--number`` options.
* By default all pairs within the cut-off distance are reported
  unless the user explicitly restricts this to a subset by
  specifying the ``--number`` option (i.e. ``--number`` now turned
  off by default).
* By default nearest distances between peaks and gene are
  calculated from the TSS of the feature to whichever of the peak
  edges are closer; alternatively distances can be calculated
  between the nearest pair of peak/gene edges by specifying the
  ``--edge=both`` option.
* Any differential expression flags in the input genes file
  are ignored unless the ``--only-DE`` option is specified, in which
  case only the differentially expressed genes are considered
  in the analyses.
* By default each peak/gene pair is reported on a separate
  line; the ``--compact`` option reports all nearest gene/peaks
  on a single line of output.
* New ``direction`` field in output indicates whether hits are
  up- or downstream from reference.
* Specify arbitrary columns from input peaks file using new
  ``--peak_cols`` options to set chromosome, start and end.
* Output file names now end with ``gene-centric`` and
  ``peak-centric``.
* Excel output is only produced if the ``--xlsx`` option is
  specified; spreadsheets are now output in XLSX format (instead
  of XLS).
* Summary output is only produced if ``--summary`` is specified.
* The ``rearrange_columns.py`` utility has been dropped.

----------------------------------
Version 0.5.0-alpha.8 (2016-02-18)
----------------------------------

* Fix typo in XLSX 'notes' sheet.

----------------------------------
Version 0.5.0-alpha.7 (2016-02-03)
----------------------------------

* Update ``--xls`` option to ``--xlsx`` and generate XLSX
  files (instead of XLS); as XLSX has much greater limits on
  the number of rows and columns allowed in a worksheet
  this should address previous problems with having data
  split over multiple sheets.
* Correct headers and placeholders now output when using
  ``--compact`` option when ``--number`` is not specified.
* Pre-existing output files are explicitly removed before
  analysis is run (rather than relying on overwrite).

----------------------------------
Version 0.5.0-alpha.6 (2016-01-29)
----------------------------------

* Fix broken ``--xls`` option (crashed program if specified).

----------------------------------
Version 0.5.0-alpha.5 (2016-01-29)
----------------------------------

* Rename 'feature' to 'gene' in program output, documentation etc
  NB this doesn't affect the program function.
* By default all pairs within the cut-off distance are reported
  unless the user explicitly restricts this to a subset by
  specifying the ``--number`` option (i.e. ``--number`` now turned off
  by default).
* Output file names changed to ``feature-centric`` and ``peak-centric``.
* Options are grouped into subsets when displayed by ``-h/--help``.
* Parameter defaults are also given in the documentation.
* Peaks in the input have 'start' and 'end' positions which
  aren't at least 1bp apart cause the program to raise an error.

----------------------------------
Version 0.5.0-alpha.4 (2015-12-01)
----------------------------------

* Fix the broken ``--promoter_region`` option which was being
  ignored.

----------------------------------
Version 0.5.0-alpha.3 (2015-12-20)
----------------------------------

* ``--compact`` now only changes the output format from "multi-line"
  (i.e. one hit pair per line) to "single-line" (i.e. all hits on
  the same line). The same fields are reported in both modes.
* The explanatory text for the dist_closest field has been updated
  to make it clearer what this means.

----------------------------------
Version 0.5.0-alpha.2 (2015-10-28)
----------------------------------

* Executable now installs as ``RnaChipIntegrator`` (i.e. no ``.py``
  extension)
* Specify feature type (e.g. ``gene``, ``transcript`` etc) to be used
  in output using ``--feature`` option.
* New ``direction`` field in output indicates whether hits are
  up- or downstream from reference.
* Specify arbitrary columns from input peaks file using new
  ``--peak_cols`` options to set chromosome, start and end.
* ``--pad`` option is automatically implied by the ``--compact``
  option (i.e. single line output is always padded).

----------------------------------
Version 0.5.0-alpha.1 (2015-09-01)
----------------------------------

* Complete reimplementation of RnaChipIntegrator to unify internal
  algorithms, simplify usage and substantially update the
  documentation.
* No distinction is now made between 'summits' and 'peaks'; the
  same algorithm is applied in each case.
* The program always finds the nearest features to each peak, and
  vice versa. The same distance cutoff and maximum number of hits
  are applied to both and can be specified using the ``--cutoff``
  and ``--number`` options.
* By default nearest distances between peaks and features are
  calculated from the TSS of the feature to whichever of the peak
  edges are closer; alternatively distances can be calculated
  between the nearest pair of peak/feature edges by specifying the
  ``--edge=both`` option.
* Any differential expression flags in the input features file
  are ignored unless the ``--only-DE`` option is specified, in which
  case only the differentially expressed features are considered
  in the analyses.
* By default each peak/feature pair is reported on a separate
  line; the ``--compact`` option reports all nearest features/peaks
  on a single line of output.
* Excel output is only produced if the ``--xls`` option is specified;
  summary output is only produced if ``--summary`` is specified.
* The ``rearrange_columns.py`` utility has been dropped.

--------------------------
Version 0.4.4 (2015-06-10)
--------------------------

* Use ``/usr/bin/env`` rather than ``/bin/env`` to invoke Python
  interpreter in RnaChipIntegrator.py (was broken for e.g. Ubuntu
  linux).

--------------------------
Version 0.4.3 (2014-05-08)
--------------------------

* Update ``--pad`` output so that requested number of lines appears
  for peaks even when there are no hits, and "empty" lines contain
  the chromosome, start and end positions for the peak in question.

--------------------------
Version 0.4.2 (2014-05-02)
--------------------------

* Truncate worksheet titles if they exceed maximum length as defined by
  the spreadsheet writing libraries.

--------------------------
Version 0.4.1 (2014-01-20)
--------------------------

* Add ``--pad`` option: for 'NearestTranscriptsToPeakEdge' and
  'NearestTSSToPeakEdge' analyses, where necessary adds blank lines to
  output files and spreadsheet so that each reported peak has the same
  number of lines associated regardless of the number of hits.

--------------------------
Version 0.4.0 (2014-01-20)
--------------------------

* Fixed bug in overlap determination, which manifested when a gene was on
  the negative strand *and* was also wider than the peak. In those cases
  the start and end of the gene were being assigned incorrectly way around.

  (The bug didn't affect results for other genes on the negative strand
  which were narrower than the peak.)

  Note that this bug would have a similar effect on determining whether a
  peak was within the promoter region of a gene on the negative strand.
  However the lists of nearest genes/peaks were not affected and the results
  should otherwise have been correct.

--------------------------
Version 0.3.3 (2012-02-16)
--------------------------

* Added explanatory text to the "notes" page of the output XLS spreadsheet
  and standardised naming of output files to match XLS page titles.
* Minor updates to READMEs/documentation.

--------------------------
Version 0.3.2 (2012-01-27)
--------------------------

* Output files now use ``<Rna-Seq-file>_vs_<ChIP-Seq-file>``
  as the default basename (unless overridden by the
  ``--project`` option).
* Added example data files in new ``examples`` directory.

--------------------------
Version 0.3.1 (2012-01-20)
--------------------------

* Added ``setup.py`` into an installable Python package.
* Updated documentation.

--------------------------
Version 0.3.0 (2012-01-05)
--------------------------

* Rename ``ID`` column to ``geneID`` (using ``ID`` has the
  potential to clash with other programs where this is a
  reserved word).
* Various improvements to some of the column descriptions
  on the "notes" page of the output XLS file.
* In all analyses, now only use those genes flagged as
  differentially expressed (use all if no flag was specified
  on the input gene data).

--------------------------
Version 0.2.0 (2011-12-19)
--------------------------

* Only performs analyses which are appropriate for the supplied ChIP peak
  data i.e. ignore "region"-based analyses if ChIP data are summits, or
  summit-based analyses if data are regions.

--------------------------
Version 0.1.4 (2011-12-08)
--------------------------

* Program will stop if it encounters any 'bad' lines in the RNA-seq/transcipt
  input data, with the exception of the first line (which is treated as a
  header and skipped if it contains bad data).
* New option ``--no-xls``: suppresses output of XLS spreadsheet.

--------------------------
Version 0.1.3 (2011-12-07)
--------------------------

* Skip input transcripts where 'start' position is higher than 'end'.
* In output spreadsheet, splits the lists of ``transcripts inbetween``
  across multiple columns in the ``TSSToSummits`` sheet if they exceed 250
  characters, and creates multiple sheets for result sets that exceed 65536
  rows.

--------------------------
Version 0.1.2 (2011-12-05)
--------------------------

* Fixed failure when using with Python 2.4 (``optparse.OptionParser``
  "epilog" argument is unsupported)

--------------------------
Version 0.1.1 (2011-11-24)
--------------------------

* Updated to use ``optparse`` library to process command line arguments,
  and substantially expanded help text (available using ``-h`` or
  ``--help`` option).

--------------------------
Version 0.1.0 (2011-11-21)
--------------------------

* Baseline version of ``RnaChIPIntegrator.py``.


