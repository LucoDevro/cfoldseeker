#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md") as readme:
    long_description = readme.read()

setup(name = "cfoldseeker",
      version = "0.0.0",
      author="Lucas De Vrieze",
      author_email="lucas.devrieze@kuleuven.be",
      license = "MIT",
      description = "Find gene clusters via structural similarity",
      packages = find_packages(),
      long_description = long_description,
      long_description_content_type = "text/markdown",
      python_requires = ">=3.12.0",
      classifiers = [
          "Programming Language :: Python :: 3.12",
          "License :: OSI Approved :: MIT License",
          "Operating System :: OS Independent",
          "Intended Audience :: Science/Research",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
      ],
      entry_points = {"console_scripts": ['cfoldseeker = cfoldseeker.main:main',
                                          'cfoldseeker-cds = cfoldseeker.build_cds_db:main']},
      install_requires=[
          "biopython",
          "cblaster >=1.3.20",
          "polars >=1.0.0",
          "networkx",
          "requests",
          "fastexcel"
      ],
      )
