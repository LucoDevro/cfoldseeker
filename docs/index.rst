.. cfoldseeker documentation master file, created by
   sphinx-quickstart on Tue Apr 14 13:00:12 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

cfoldseeker
=======================================

**Welcome to cfoldseeker's documentation!**

|docs|
|downloads|
|bioconda|
|docker|
|pypi|
|doi|

.. |docs| image:: https://img.shields.io/readthedocs/cfoldseeker/latest?style=flat-square&maxAge=600&logo=readthedocs
   :target: https://cfoldseeker.readthedocs.io/en/latest/

.. |bioconda| image:: https://img.shields.io/conda/vn/bioconda/cfoldseeker?style=flat-square&maxAge=3600&logo=anaconda
   :target: https://anaconda.org/bioconda/cfoldseeker

.. |docker| image:: https://img.shields.io/docker/v/lucodevro/cfoldseeker?sort=semver&label=quay.io&logo=docker
   :target: https://quay.io/repository/biocontainers/cfoldseeker

.. |pypi| image:: https://img.shields.io/pypi/v/cfoldseeker?sort=semver&logo=pypi
   :target: https://pypi.org/project/cfoldseeker/

.. |downloads| image:: https://anaconda.org/bioconda/cfoldseeker/badges/downloads.svg
   :target: https://bioconda.github.io/recipes/cfoldseeker/README.html#download-stats

.. |doi| image:: https://zenodo.org/badge/904110273.svg
   :target: https://doi.org/10.5281/zenodo.20746182

*cfoldseeker: Identify gene clusters using protein structure similarity.*

``cfoldseeker`` brings protein structure similarity to gene cluster mining! Starting from a query gene cluster of which you have all protein structure models, **cfoldseeker** searches all genomically colocalised structural homologs in a protein structure database.

**The remote mode** leverages the FoldSeek webserver to search against the AlphaFoldDB, and uses several cross-referencing APIs (UniProt ID mapping, KEGG, ENA GenPept) to retrieve the genomic location of every hit.

**The local mode** launches FoldSeek searches against a local structure database, and uses a genomic context database prepared from your local sequences using its helper tool ``cfoldseeker-cds``.

**The local-clustered mode** extends local mode searches for huge protein structure databases by searching against representative proteins obtained through preclustering the original sequences with ``mmseqs2``.
   
If you find ``cfoldseeker`` useful, please cite:

::

	In preparation


.. toctree::
   :caption: User Guide
   :maxdepth: 2

   guide/index
   guide/install
   guide/usage
   guide/tutorial


Comprehensive documentation for all API exposed by ``cfoldseeker``:

.. toctree::
   :caption: API Reference
   :maxdepth: 1

   api_ref

