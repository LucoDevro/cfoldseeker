User guide
============

``cfoldseeker`` has several search modes and helper tools, each one requiring different prior work to be done with MMseqs and/or FoldSeek.

.. tip::

   Each ``cfoldseeker`` command prints logs at the stdout. You can save them in a log file by ``tee``-ing.

   .. code-block:: bash

   		cfoldseeker ... | tee log_file.log

Remote mode
-----------
The remote mode requires from you a set of protein structural models, and our copy of the UniProt ID mapping table (only the KEGG and GenPept rows of `the official mapping table <https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz>`_).

Prior work
~~~~~~~~~~
Get structural models as CIF files of your query proteins, either experimentally or computationally (AlphaFold, ESMFold, OpenFold...). Collect them all in one folder ``query_models``.

The search
~~~~~~~~~~
Run ``cfoldseeker`` in remote mode using our UniProt mapping table ``uniprot_kegg_genpept.gz`` and save results in a new folder ``results``. Search the AFDB50, the AFDB-Proteome, the AFDB-SwissProt databases, using 4 workers for the cross-referencing APIs.

.. code-block:: bash

	cfoldseeker -m remote -q query_models -o results -rdb afdb50 afdb-proteome afdb-swissprot -w 4 -uma uniprot_kegg_genpept.gz

Local mode
----------
The local mode requires from you a set of query protein structures, a genomic context table made using ``cfoldseeker-cds``, and a FoldSeek target database.

Prior work
~~~~~~~~~~
**1.** Get structural models as CIF files of your query proteins, either experimentally or computationally (AlphaFold, ESMFold, OpenFold). Collect them all in one folder ``query_models``.

**2.** Build a genomic context table using our provided helper tool ``cfoldseeker-cds``, which builds a context table directly from a set of NCBI or Bakta GFF files, or from a folder holding an extracted NCBI package of GFF files.

.. tip::

   ``cfoldseeker-cds`` usually produces context DBs populated with accession codes. Although using accession codes standardises analysis outputs greatly, they are not very human-friendly. Enable the ``-tn`` flag to make ``cfoldseeker-cds`` populate the context DB with **readily readable taxon names**. These names will recur in the outputs of a ``cfoldseeker`` analysis.

   In both NCBI modes, ``cfoldseeker-cds`` will fetch the taxon names from NCBI via the Entrez API. In Bakta mode, it will generate a generic taxon name. In TSV mode, it will trust the user's inputs.

To produce a compressed genomic context DB ``ncbi_package_db.gz`` from a default NCBI package (a folder ``ncbi_dataset``, which has an identically named subfolder):

.. code-block:: bash

	cfoldseeker-cds -i ncbi-dataset -o ncbi_package_db.gz -m ncbi-package -gz

To produce a compressed context DB ``ncbi_files_db.gz`` populated with taxon names from a folder of NCBI GFF files ``gffs``:

.. code-block:: bash

	cfoldseeker-cds -i gffs -o ncbi_files_db.gz -m ncbi-gff -gz -tn

.. note::

   Keep in mind that ``cfoldseeker-cds`` gets the taxa names from the GFF filenames, or from the subfolders in the NCBI package. **Give your GFF files a unique name** (e.g. NCBI accession ID)!

.. tip::

   cfoldseeker-cds DBs can be concatenated using ``cat``. No need to rerun the builder!

**3.** Generate the target DB ``target_DB`` using ``foldseek createdb``. This prepares a FoldSeek DB from your folder containing the target set of protein structures (``input``). You probably don't have thousands of protein structures lingering around, so you will need FoldSeek's builtin support of **ProstT5**, a LLM that directly translates amino acid sequences into FoldSeek's internal 3Di alphabet, skipping protein model prediction. **This is the key step that allows searching sequence databases using structural similarity.**

First make sure you have downloaded ProstT5's weights.

.. code-block:: bash

	foldseek databases ProstT5 <path-to-prostt5-weights> tmp

Then you can build the FoldSeek DB directly from your protein sequences in ``input``.

.. code-block:: bash

	foldseek createdb input target_DB --prostt5-model <path-to-prostt5-weights>

.. warning::

   This is a **time- and computationally demanding** task! Consider using a GPU (by adding the ``--gpu 1`` flag if you have the hardware configured).

.. tip::

   FoldSeek offers a command to **merge existing DBs** ``foldseek concatdbs``. Use it to concat existing target DBs, and save time and computational work.

The search
~~~~~~~~~~
Run ``cfoldseeker`` in local mode using FoldSeek DB ``target_DB``, and context DB ``cds_db.gz``. Save results in a new folder ``results``.

.. code-block:: bash

	cfoldseeker -m local -q query_models -o results -ldb target_DB/target_DB -cdb cds_db.gz

Local-clustered mode
--------------------
The local-clustered mode requires a set of query protein structures, a genomic context table made using ``cfoldseeker-cds``, a MMseqs2 clustering TSV file, and a FoldSeek target database of the representative proteins.

Prior work
~~~~~~~~~~
**1.** Get structural models as CIF files of your query proteins, either experimentally or computationally (AlphaFold, ESMFold, OpenFold...). Collect them all in one folder ``query_models``.

**2.** Build a genomic context table using our provided helper tool ``cfoldseeker-cds``, which builds a context table directly from a set of NCBI or Bakta GFF files, or from a folder holding an extracted NCBI package of GFF files. *(See the prior work section of local search above)*

**3.** Cluster your target sequences in the folder ``input_all`` using ``mmseqs``. You can do this using ``easy-cluster``, or ``easy-linclust`` for huge sequence databases. We recommend to use an identity and a coverage threshold of 90 % to ensure all proteins in a protein cluster have identical functions. Use other thresholds at your own risk!

.. code-block:: bash

	mmseqs easy-linclust input_all clustered tmp --min-seq-id 0.9 -c 0.9

This will, among others, produce a fasta file ``clustered_rep_seq.fasta`` containing the amino acid sequences of the representative protein of each cluster.

**4.** Generate the target DB from the representative sequences using FoldSeek and ProstT5 (see also the prior work section of local search above). Make sure you have downloaded ProstT5's weights.

.. code-block:: bash

	foldseek createdb clustered_rep_seq.fasta target_DB --prostt5-model <path-to-prostt5-weights>

The search
~~~~~~~~~~
Run ``cfoldseeker`` in local mode using FoldSeek DB ``target_DB``, context DB ``cds_db.gz``, and preclustering table ``clustered_table.tsv``. Save results in a new folder ``results``.

.. code-block:: bash

	cfoldseeker -m local_clustered -q query_models -o results -ldb target_DB/target_DB -cdb cds_db.gz -scl clustered_table.tsv

Specifying search options
-------------------------
``cfoldseeker`` offers several filtering thresholds to refine your hit set and find relevant gene clusters.

General options
~~~~~~~~~~~~~~~
The general search options are a mix of what ``cblaster`` and ``foldseek`` offer. These options are listed below.

+-------------------+------------------------------------------------------------+
| **filter**        | **Description**                                            |
+-------------------+------------------------------------------------------------+
| ``--max-eval``    | Maximum E-value of a Foldseek hit                          |
+-------------------+------------------------------------------------------------+
| ``--min-score``   | Minimum bitscore of a Foldseek hit                         |
+-------------------+------------------------------------------------------------+
| ``--min-seqid``   | Minimum sequence identity between a hit and a query (in %) |
+-------------------+------------------------------------------------------------+
| ``--min-qcov``    | Minimum query coverage of the hit (in %)                   |
+-------------------+------------------------------------------------------------+
| ``--min-tcov``    | Minimum target coverage of the hit (in %)                  |
+-------------------+------------------------------------------------------------+
| ``--max-gap``     | Maximum gap between two hits on the same scaffold (in bp)  |
+-------------------+------------------------------------------------------------+
| ``--max-length``  | Maximum length of a cluster (in bp)                        |
+-------------------+------------------------------------------------------------+
| ``--min-hits``    | Minimum number of hits in a cluster                        |
+-------------------+------------------------------------------------------------+
| ``--min-cov-qrs`` | Minimum number of queries represented in a cluster         |
+-------------------+------------------------------------------------------------+
| ``--require``     | Queries that must have a hit in a cluster                  |
+-------------------+------------------------------------------------------------+

Getting all cluster layouts
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Sometimes, a protein may match with multiple query proteins, e.g. when you have two paralogs among your query proteins. This makes it tricky to determine what the query layout of an identified cluster is. For example, if two proteins in a cluster both match with query proteins 1 and 2, cluster layout *12* is equally correct as layout *21*. By default, if ``cfoldseeker`` encounters the identical cluster passing the filtering thresholds with different layouts, it will only keep the one with the highest cluster score.

If you are interested in all possible cluster layouts passing your filtering thresholds rather than only the highest-scoring one, you can turn on the ``all-layouts`` flag and keep all passing configurations.

Using taxon filters
~~~~~~~~~~~~~~~~~~~
The Foldseek webserver offers an interface to filter hits by a taxonomic filter. ``cfoldseeker`` exposes this interface in its remote mode via the ``-tf`` or ``--taxon-filter`` flag. Foldseek expects NCBI taxon IDs for its taxonomic filter. If you are not sure what the exact taxon ID of your taxa group is, you can check it at `the NCBI Taxonomy website <https://www.ncbi.nlm.nih.gov/datasets/taxonomy/browser/>`_.

Another option is to add as a taxonomic filter under *settings* on the Foldseek webserver. Click then the *API* button on the top right, which will give you a pop-up with code to submit your query via the terminal, including the taxonomy ID of your taxa group among the settings.

Specifying outputs
------------------
``cfoldseeker`` can produce several outputs. By default, it produces only hit and cluster overview tables in TSV format, but several cblaster-style outputs and the raw Foldseek hit tables can be returned on request.

Cluster table
~~~~~~~~~~~~~
The ``clusters.tsv`` file is a tab-separated file summarising the properties of the identified clusters. It comprises the 10 columns described below.

+-------------+------------------------------------------------------+
| **Column**  | **Description**                                      |
+-------------+------------------------------------------------------+
| number      | Arbitrary unique number                              |
+-------------+------------------------------------------------------+
| hits        | IDs of the hits part of this cluster                 |
+-------------+------------------------------------------------------+
| start       | Starting coordinate of the entire cluster            |
+-------------+------------------------------------------------------+
| end         | End coordinate of the entire cluster                 |
+-------------+------------------------------------------------------+
| length      | Sum of the lengths of all exons part of this cluster |
+-------------+------------------------------------------------------+
| score       | Sum of the Foldseek bitscores of all cluster members |
+-------------+------------------------------------------------------+
| scaff       | ID of the scaffold/contig harbouring this cluster    |
+-------------+------------------------------------------------------+
| strand      | Strand location of the cluster                       |
+-------------+------------------------------------------------------+
| taxon_name  | Name of the taxon (e.g. NCBI Assembly ID)            |
+-------------+------------------------------------------------------+
| taxon_id    | Unique taxon ID (e.g. NCBI taxon ID)                 |
+-------------+------------------------------------------------------+

Hit table
~~~~~~~~~
The ``hits.tsv`` file gathers metadata about all hits part of the identified clusters. It contains the 16 columns described below.

+-----------------+----------------------------------------------------------------------------+
| **Column**      | **Description**                                                            |
+-----------------+----------------------------------------------------------------------------+
| db_id           | Unique hit ID                                                              |
+-----------------+----------------------------------------------------------------------------+
| query           | ID of the query with which the hit matches                                 |
+-----------------+----------------------------------------------------------------------------+
| scaff           | ID of the scaffold/contig harbouring this hit                              |
+-----------------+----------------------------------------------------------------------------+
| strand          | Strand location of the hit                                                 |
+-----------------+----------------------------------------------------------------------------+
| coords          | Comma-separated list of coordinate intervals for all exons of this protein |
+-----------------+----------------------------------------------------------------------------+
| db              | DB in which this hit was found (for local DBs: *local*)                    |
+-----------------+----------------------------------------------------------------------------+
| crossref_id     | ID of the cross-referenced record (same as db_id for local DBs)            |
+-----------------+----------------------------------------------------------------------------+
| crossref_method | Cross-referencing method (*KEGG*, *GenPept*, *WGS-GenPept*, or *local*)    |
+-----------------+----------------------------------------------------------------------------+
| name            | Protein annotation                                                         |
+-----------------+----------------------------------------------------------------------------+
| taxon_name      | Name of the taxon (e.g. NCBI Assembly ID)                                  |
+-----------------+----------------------------------------------------------------------------+
| taxon_id        | Unique taxon ID (e.g. NCBI taxon ID)                                       |
+-----------------+----------------------------------------------------------------------------+
| evalue          | Hit e-value                                                                |
+-----------------+----------------------------------------------------------------------------+
| score           | Hit bitscore                                                               |
+-----------------+----------------------------------------------------------------------------+
| seqid           | Sequence identity between hit and query protein                            |
+-----------------+----------------------------------------------------------------------------+
| qcov            | Query coverage of the hit (which fraction of the query fits)               |
+-----------------+----------------------------------------------------------------------------+
| tcov            | Target coverage of the hit (which fraction of the target fits)             |
+-----------------+----------------------------------------------------------------------------+

Foldseek output
~~~~~~~~~~~~~~~
``cfoldseeker`` can return the raw Foldseek output from which it started. This can be useful if you want to track down why a certain hit was not found being part of a cluster, or how many hits have been found for a certain query protein. In remote mode, ``cfoldseeker`` returns the json files it received from the Foldseek webserver. In the local modes, it will return the tab-separated text file returned by the local ``foldseek`` call.

*cblaster* outputs
~~~~~~~~~~~~~~~~
``cfoldseeker`` has tightly integrated ``cblaster``. All results are cast into a cblaster session, from which familiar outputs can be obtained, such as the summary table, the binary table, the hit plot, and the clinker alignment. See `the cblaster documentation <https://cblaster.readthedocs.io/en/latest/guide/search_module.html#specifying-output>`_ for specifics on these outputs.


