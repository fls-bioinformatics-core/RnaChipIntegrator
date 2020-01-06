.. _additional_info:

Additional information
======================

Technical details
-----------------

``RnaChipIntegrator`` has been tested against the following versions of
Python:

 * 2.7
 * 3.5
 * 3.6
 * 3.7
 * 3.8

It requires the external ``xlsxwriter`` library in order to generate the
``.xlsx`` Excel spreadsheets:

 * http://xlsxwriter.readthedocs.io/index.html

This library will be installed automatically if using ``pip``.

The source code for the program is hosted on GitHub at

 * https://github.com/fls-bioinformatics-core/RnaChipIntegrator

The ``devel`` branch holds the developmental code.

The documentation under the ``docs`` subdirectory is generated using the
``sphinx`` package, and can be built by doing either::

    python setup.py sphinx_build

or::

    cd docs
    make html

both of which create the documentation in the ``docs/_build``
subdirectory.

A copy of the documentation is also hosted on ReadTheDocs at
http://rnachipintegrator.readthedocs.io/en/latest/

The unit tests for the code can be run using::

    python setup.py test

(Note that this requires the ``nose`` package.)

The tests are also run on the Travis-CI continuous integration
server each time a change is made to the code; the test results
can be found at
https://travis-ci.org/fls-bioinformatics-core/RnaChipIntegrator/


Credits
-------

``RnaChipIntegrator`` was written by Peter Briggs, Ian Donaldson
and Leo Zeef in the Bioinformatics Core Facility (BCF) in the
Faculty of Life Sciences, University of Manchester, with
additional contributions from Casey Bergman.


Licensing
---------

This software is licensed under the Artistic License 2.0; see
the ``LICENSE`` document.


Version history and changes
---------------------------

See the :doc:`CHANGELOG <changes>`.
