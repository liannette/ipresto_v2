#!/usr/bin/env python
import os

from setuptools import setup
from setuptools import find_packages

here = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(here, "ipresto", "__version__.py")) as f:
    exec(f.read(), version)

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="iPRESTO",
    version=version["__version__"],
    author="Joris Louwen",
    author_email="jorislouwen@hotmail.com",
    description="Detection of biosynthetic sub-clusters",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://git.wageningenur.nl/bioinformatics/iPRESTO",
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
          'biopython',
          'matplotlib',
          'networkx',
          'numpy',
          'gensim',
          'pyLDAvis',
          'pandas',
          'scipy',
          'seaborn',
          'statsmodels',
          'sympy'
      ],
    scripts=[os.path.join(here, 'ipresto.py')]
)
