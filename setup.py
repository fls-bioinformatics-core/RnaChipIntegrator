"""Description

Setup script to install RnaChipIntegrator: analyses of peak data with
genomic feature data

Copyright (C) University of Manchester 2011-20,2024 Peter Briggs,
Ian Donaldson, Leo Zeef

"""

from setuptools import setup
import rnachipintegrator

# Package version
version = rnachipintegrator.get_version()

download_url = 'https://github.com/fls-bioinformatics-core/RnaChipIntegrator/archive/v' + version + '.tar.gz'

with open('README.rst') as fh:
    readme = fh.read()

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
    description = "Analyse genes against peak data, and vice versa",
    long_description = readme,
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
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    install_requires = ['xlsxwriter >= 0.8.4',],
    test_suite = 'nose.collector',
    tests_require = ['nose'],
    scripts = scripts,
    data_files = [ ('RnaChipIntegrator',
                    ['README.rst',
                     'LICENSE',
                     'INSTALL',
                     'CHANGELOG.rst']),],
    include_package_data=True,
    zip_safe = False
)
