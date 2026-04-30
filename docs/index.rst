.. CAGEcleaner documentation master file, created by
   sphinx-quickstart on Tue Apr 14 13:00:12 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

cfoldseeker
=======================================

**Welcome to cfoldseeker's documentation!**

*cfoldseeker: Identify gene clusters using protein structure similarity.*

``cfoldseeker`` brings protein structure similarity to gene cluster mining! Starting from a query gene cluster of which you have all protein structure models, **cfoldseeker** searches all genomically colocalised structural homologs in a protein structure database.

**The remote mode** leverages the FoldSeek webserver to search against the AlphaFoldDB, and uses several cross-referencing APIs (UniProt ID mapping, KEGG, ENA GenPept) to retrieve the genomic location of every hit.

**The local mode** launches FoldSeek searches against a local structure database, and uses a genomic context database prepared from your sequences using its helper tool ``cfoldseeker-cds``.

**The local-clustered mode** facilitates searches against huge structure databases by searching against the representative proteins of a target database preclustered with ``mmseqs2``.
   
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

