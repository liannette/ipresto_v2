import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="iPRESTO",
    version="1.0",
    author="Joris Louwen",
    author_email="jorislouwen@hotmail.com",
    description="Detection of biosynthetic sub-clusters",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://git.wageningenur.nl/bioinformatics/iPRESTO",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
          'Bio',
          'matplotlib',
          'networkx',
          'numpy',
          'gensim',
          'pyLDAvis',
          'pandas',
          'scipy',
          'seaborn'
          'statsmodels',
          'sympy'
          
      ],
)
