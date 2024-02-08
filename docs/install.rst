.. _install:

Getting started
===============

``RnaChipIntegrator`` is Python software which runs with Python
version 3.6 or higher.

---------------------------------------------------------------
Installing RnaChipIntegrator from PyPI using virtualenv and pip
---------------------------------------------------------------

The recommended way to get the latest version of ``RnaChipIntegrator``
is to create a Python virtual environment, and then install the
software using the ``pip`` utility.

For example: to create and activate a virtual environment called
``venv.rci`` using the ``virtualenv`` utility:

::

   virtualenv venv.rci
   source venv.rci/bin/activate

``RnaChipIntegrator`` can then be installed using ``pip``, by
running:

::

   pip install RnaChipIntegrator

which will make the ``RnaChipIntegrator`` program available.

.. note::

   If using ``RnaChipIntegrator`` from a virtual environment in
   this way, make sure to activate the environment each time
   before using it, for example:

   ::

      source venv.rci/bin/activate

To deactive the virtual environment afterwards, do ``deactivate``.

.. note::

   For an introduction to ``pip`` and ``virtualenv``, see for example:

   * http://www.dabapps.com/blog/introduction-to-pip-and-virtualenv-python/
   * https://www.biostars.org/p/109179/

To update an existing version of the program to a newer one, use:

::

    pip install -U RnaChipIntegrator

from within the activated virtual environment.

For other ways of installing please refer to the ``INSTALL`` document
included with the distribution.

----------------------------------------
Installing RnaChipIntegrator using Conda
----------------------------------------

Another approach for installing ``RnaChipIntegrator`` to use
`Conda <http://conda.pydata.org/docs/>`__
(most easily obtained via the
`Miniconda Python distribution <http://conda.pydata.org/miniconda.html>`__).

Once you have Conda installed you can create a new Conda environment
with ``RnaChipIntegrator`` installed using the following command:

::

   conda create -c bioconda -n rci rnachipintegrator

Alternatively you can install ``RnaChipIntegrator`` into an existing
Conda environment using:

::

   conda install -c bioconda rnachipintegrator

.. warning::

   It's recommended that ``RnaChipIntegrator`` be installed into a
   new Conda environment to avoid issues with incompatible packages
   (which is possible for example when trying to install directly
   into the Conda's ``base`` distribution).
