"""Description

Setup script to install RnaChipIntegrator: analyses of ChIP-Seq and RNA-Seq data

Copyright (C) University of Manchester 2011-15 Peter Briggs, Leo Zeef &
Ian Donaldson

"""

# Hack to acquire all scripts that we want to
# install into 'bin'
from glob import glob
scripts = []
for pattern in ('bin/*.py',):
    scripts.extend(glob(pattern))

# Setup for installation etc
from setuptools import setup
import rnachipintegrator
setup(
    name = "RnaChipIntegrator",
    version = rnachipintegrator.get_version(),
    maintainer = 'Peter Briggs',
    maintainer_email = 'peter.briggs@manchester.ac.uk',
    packages = ['rnachipintegrator'],
    license = 'Artistic License',
    py_modules = ['RnaChipIntegrator','Spreadsheet'],
    install_requires = ['xlwt >= 0.7.2',
                        'xlrd >= 0.7.1',
                        'xlutils >= 1.4.1'],
    test_suite = 'nose.collector',
    tests_require = ['nose'],
    scripts = scripts,
    data_files = [ ('RnaChipIntegrator',
                    ['README.markdown',
                     'LICENSE',
                     'INSTALL',
                     'ChangeLog']),
                   ('RnaChipIntegrator/doc',
                    ['doc/MANUAL.markdown',
                     'doc/rnachipintegrator_nearestEdgetoPeak.png',
                     'doc/rnachipintegrator_nearestTSStoSummit.png',
                     'doc/rnachipintegrator_nearestTSStoPeak.png']) ]
)
