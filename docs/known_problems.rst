.. _known_problems:

Known problems
==============

* :ref:`command_installed_as_lower_case`

.. _command_installed_as_lower_case:

Command installs as 'rnachipintegrator' not 'RnaChipIntegrator'
---------------------------------------------------------------

When installed, the package should provide a command called
``RnaChipIntegrator``; however in some circumstances it appears that
the command name is converted to all lower-case, and is installed as
``rnachipintegrator`` instead.

It's not clear why this happens but may be related to the version
of ``pip`` that is used to install the software: the behaviour has been
observed when using ``pip`` version 7.1.2 but not with version 9.0.1.

In these cases the workaround is to use ``rnachipintegrator`` rather
``RnaChipIntegrator`` when running the program.

(See `issue #48 <https://github.com/fls-bioinformatics-core/RnaChipIntegrator/issues/48>`_)
