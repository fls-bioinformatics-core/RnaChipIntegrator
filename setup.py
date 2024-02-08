"""
Setup script to install RnaChipIntegrator (analyses of peak data with
genomic feature data)
"""

from setuptools import setup
import codecs
import os.path

# Installation requirements
install_requires = ['xlsxwriter >= 0.8.4']

# Acquire package version for installation
# (see https://packaging.python.org/guides/single-sourcing-package-version/)
def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

RNACHIPINTEGRATOR_VERSION = get_version("rnachipintegrator/__init__.py")

DOWNLOAD_URL = "https://github.com/fls-bioinformatics-core/"\
               "RnaChipIntegrator/archive/v%s.tar.gz" %\
               RNACHIPINTEGRATOR_VERSION

# Get long description from the project README
readme = read('README.rst')

# Setup for installation etc
setup(
    name = "RnaChipIntegrator",
    version = RNACHIPINTEGRATOR_VERSION,
    description = "Analyse genes against peak data, and vice versa",
    long_description = readme,
    url = 'https://github.com/fls-bioinformatics-core/RnaChipIntegrator',
    download_url = DOWNLOAD_URL,
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
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Environment :: Console',
        'License :: OSI Approved :: Artistic License',
        'Operating System :: POSIX',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Natural Language :: English',
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    install_requires = install_requires,
    test_suite = 'nose.collector',
    tests_require = ['nose'],
    data_files = [ ('RnaChipIntegrator',
                    ['README.rst',
                     'LICENSE',
                     'INSTALL',
                     'CHANGELOG.rst']),],
    include_package_data=True,
    zip_safe = False
)
