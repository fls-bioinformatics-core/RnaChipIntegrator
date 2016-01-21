"""Description

Setup script to install RnaChipIntegrator: analyses of peak data with
genomic feature data

Copyright (C) University of Manchester 2011-16 Peter Briggs, Leo Zeef &
Ian Donaldson

"""

from setuptools import setup
import rnachipintegrator

# Package version
version = rnachipintegrator.get_version()

download_url = 'https://github.com/fls-bioinformatics-core/RnaChipIntegrator/archive/v' + version + '.tar.gz'

# Hack to acquire all scripts that we want to
# install into 'bin'
from glob import glob
scripts = []
for pattern in ('bin/*.py',):
    scripts.extend(glob(pattern))

# Setup for installation etc
setup(
    name = "RnaChipIntegrator",
    version = version,
    description = "Integrate genomic features with expression data",
    long_description = """
Utility for integrating sets of genomic features (e.g. canonical genes,
CpG islands) with expression data.
""",
    url = 'https://github.com/fls-bioinformatics-core/RnaChipIntegrator',
    download_url = download_url,
    author = 'Peter Briggs, Ian Donaldson, Leo Zeef',
    author_email = 'peter.briggs@manchester.ac.uk',
    maintainer = 'Peter Briggs',
    maintainer_email = 'peter.briggs@manchester.ac.uk',
    packages = ['rnachipintegrator'],
    entry_points = { 'console_scripts': [
        'RnaChipIntegrator = rnachipintegrator.cli:main',]
    },
    license = 'Artistic License',
    keywords = ['RnaChipIntegrator',],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Environment :: Console',
        'License :: OSI Approved :: Artistic License',
        'Operating System :: POSIX',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
    ],
    install_requires = ['xlwt >= 0.7.2',
                        'xlrd >= 0.7.1',
                        'xlutils >= 1.4.1'],
    test_suite = 'nose.collector',
    tests_require = ['nose'],
    scripts = scripts,
    data_files = [ ('RnaChipIntegrator',
                    ['README.rst',
                     'LICENSE',
                     'INSTALL',
                     'ChangeLog']),],
    include_package_data=True,
    zip_safe = False
)
