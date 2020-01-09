import setuptools
import glob

with open("README.md", "r") as fh:
    long_description = fh.read()

#add all scripts to path, to be able to run from commandline
all_scripts = glob.glob('**/*[A-Za-z0-9].py', recursive=True)

setuptools.setup(
    name="iPRESTO",
    version="1.0.3",
    author="Joris Louwen",
    author_email="jorislouwen@hotmail.com",
    description="Detection of biosynthetic sub-clusters",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://git.wageningenur.nl/bioinformatics/iPRESTO",
    packages=setuptools.find_packages(),
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
    scripts=all_scripts
)
