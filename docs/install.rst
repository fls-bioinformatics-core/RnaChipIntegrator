.. _install:

Getting started
===============

The easiest way to get the latest version of ``RnaChipIntegrator`` is to
use Python's ``pip`` utility to install the latest version of the program
directly from the `Python Package Index (PyPI)
<https://pypi.python.org/pypi/>`_, by doing::

    pip install RnaChipIntegrator

.. note::

   You may need to have root privileges to install to the system
   directories, in which case preface this command with ``sudo``
   i.e.::

       sudo pip install RnaChipIntegrator

   or you can do::

       pip install --user RnaChipIntegrator

   to install it under your home area.

Alternatively you can use Python's ``virtualenv`` mechanism to install
a non-root version (this example creates one under ``.venv``)::

    virtualenv .venv; . .venv/bin/activate
    pip install RnaChipIntegrator

.. note::

   For an introduction to ``pip`` and ``virtualenv``, see for example:

   * http://www.dabapps.com/blog/introduction-to-pip-and-virtualenv-python/
   * https://www.biostars.org/p/109179/

To update an existing version of the program to a newer one, use::

    pip install -U RnaChipIntegrator

For other ways of installing please refer to the ``INSTALL`` document
included with the distribtion.

