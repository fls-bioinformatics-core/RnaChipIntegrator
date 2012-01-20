"""Description

Setup script to install RnaChipIntegrator: analyses of ChIP-Seq and RNA-Seq data

Copyright (C) University of Manchester 2011-12 Peter Briggs, Leo Zeef & Ian Donaldson

This code is free software; you can redistribute it and/or modify it
under the terms of the Artistic License 2.0 (see the file LICENSE included
with the distribution).

This will install RnaChipIntegrator and its dependencies and provide the
'RnaChipIntegrator' command.
"""

from setuptools import setup
import RnaChipIntegrator
setup(
    name = "RnaChipIntegrator",
    version = RnaChipIntegrator.__version__,
    maintainer = 'Peter Briggs',
    maintainer_email = 'peter.briggs@manchester.ac.uk',
    license = 'Artistic License 2.0',
    py_modules = ['RnaChipIntegrator','Spreadsheet'],
    install_requires = ['xlwt >= 0.7.2',
                        'xlrd >= 0.7.1',
                        'xlutils >= 1.4.1'],
    entry_points = { 'console_scripts':
                         ['RnaChipIntegrator = RnaChipIntegrator:main']
                     },
    data_files = [ ('RnaChipIntegrator', ['README.markdown',
                                          'LICENCE',
                                          'INSTALL',
                                          'ChangeLog']),
                   ('RnaChipIntegrator/doc', ['doc/README.markdown',
                                              'doc/rnachipintegrator_nearestEdgetoPeak.png',
                                              'doc/rnachipintegrator_nearestTSStoSummit.png',
                                              'doc/rnachipintegrator_nearestTSStoPeak.png']) ]
)
