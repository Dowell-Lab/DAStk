#!/usr/bin/env python

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='DAStk',
    version='0.1.1',
    description='Differential ATAC-seq toolkit',
    long_description=long_description,
    license='BSD',
    url='https://biof-git.colorado.edu/dowelllab/DAStk',
    author='Ignacio Tripodi',
    author_email='ignacio.tripodi@colorado.edu',
    python_requires='>=2.7',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
    ],

    keywords='bioinformatics genomics chromatin ATAC-seq motif transcription_factor',

    package_dir={'DAStk' : 'DAStk'},
    packages=['DAStk'],

    install_requires=[
        'argparse',
        'datetime',
        'numpy',
        'matplotlib',
        'adjustText',
        'scipy',
        'multiprocessing',
        'pandas',
    ],

    entry_points={
        'console_scripts': [
            'process_atac=DAStk:process_atac',
            'differential_md_score=DAStk:differential_md_score',
        ],
    },

    project_urls={
        'Bug Reports': 'https://biof-git.colorado.edu/dowelllab/DAStk/issues',
        'Source': 'https://biof-git.colorado.edu/dowelllab/DAStk',
        'Human Motif Sites': 'http://dowell.colorado.edu/pubs/DAStk/human_motifs.tar.gz',
        'Mouse Motif Sites': 'http://dowell.colorado.edu/pubs/DAStk/mouse_motifs.tar.gz',
    },
)
