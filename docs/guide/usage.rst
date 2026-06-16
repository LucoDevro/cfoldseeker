User guide
============

``cfoldseeker`` has several search modes and helper tools, each one requiring different prior work to be done with MMseqs and/or FoldSeek.

.. tip::

   Each ``cfoldseeker`` command prints logs at the stdout. You can capture them in a log file by ``tee``-ing.

   .. code-block:: bash

   		cfoldseeker ... | tee log_file.log

Executing searches
------------------

Remote mode
~~~~~~~~~~~
The remote mode requires from you a set of protein structural models, and `our reduced copy <https://github.com/LucoDevro/cfoldseeker/releases>`_ of the UniProt ID mapping table (`uniprot_kegg_genpent.gz`, which has only the KEGG and GenPept rows of `the official mapping table <https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz>`_).

Prior work
^^^^^^^^^^
Get structural models as **CIF files** of your query proteins, either experimentally or computationally (AlphaFold, ESMFold, OpenFold...). Collect them all in one folder ``query_models``.

The search
^^^^^^^^^^
Run ``cfoldseeker`` in remote mode using our UniProt mapping table ``uniprot_kegg_genpept.gz`` and save results in a new folder ``results``. Search the AFDB50, the AFDB-Proteome, the AFDB-SwissProt databases, using 4 workers for the cross-referencing APIs.

.. code-block:: bash

	cfoldseeker -m remote -q query_models -o results -rdb afdb50 afdb-proteome afdb-swissprot -w 4 -uma uniprot_kegg_genpept.gz

Local mode
~~~~~~~~~~
The local mode requires from you a set of query protein structures, a genomic context table made using ``cfoldseeker-cds``, and a FoldSeek target database.

Prior work
^^^^^^^^^^
**1.** Get structural models as CIF files of your query proteins, either experimentally or computationally (AlphaFold, ESMFold, OpenFold). Collect them all in one folder ``query_models``.

**2.** Build a genomic context table using our provided helper tool ``cfoldseeker-cds``, which builds a context table directly from a folder containing NCBI or Bakta Genbank files, or from an NCBI package of Genbank files.

.. tip::

   ``cfoldseeker-cds`` usually produces context DBs populated with filelabels sourced from the Genbanks' filenames. Typically, these are some database's accession codes. Although using accession codes standardises analysis outputs greatly, they are not very human-friendly. You can use one of the ``-tn*`` flags to make ``cfoldseeker-cds`` populate the context DB with **readily readable taxon names**. These names will then be used in the outputs of a ``cfoldseeker`` analysis.

   In both NCBI modes, ``cfoldseeker-cds`` will fetch the taxon names, either from NCBI via the Entrez API (``-tna``) or from a local mapping file (``-tnf``). In Bakta mode, it will generate a generic taxon name or read them in from a local mapping file. In TSV mode, it will trust the user's inputs as documented in the TSV files.

.. warning::

   **Replacing taxon names breaks the direct link between the sequence files and the hits!** Don't use this when you intend to do downstream analyses (e.g. extract gene cluster Genbank files later on, hit dereplication). You can keep two versions of the context DB: one with intact filelabels for downstream processing, and one with human-readable taxon names for reporting.

   For more information about ``cfoldseeker-cds``, head over to `its user guide <https://cfoldseeker.readthedocs.io/en/latest/guide/usage.html#context-DB-construction>`_.

To produce a compressed genomic context DB ``ncbi_package_db.gz`` from a default NCBI package (a folder ``ncbi_dataset``, which has an identically named subfolder):

.. code-block:: bash

	cfoldseeker-cds -i ncbi-dataset -o ncbi_package_db.gz -m ncbi-package -gz

To produce a compressed context DB ``ncbi_files_db.gz`` populated with taxon names fetched from NCBI from a folder of NCBI Genbank files ``gbffs``:

.. code-block:: bash

	cfoldseeker-cds -i gbffs -o ncbi_files_db.gz -m ncbi-gbff -gz -tna

Or with taxon names fetched from a local file:

.. code-block:: bash

	cfoldseeker-cds -i gbffs -o ncbi_files_db.gz -m ncbi-gbff -gz -tnf <local-taxon-name-file>

.. note::

   Keep in mind that ``cfoldseeker-cds`` gets the taxon names from the Genbank filenames, or from the subfolders in the NCBI package. **Give your Genbank files a unique name** (e.g. NCBI accession ID)!

.. tip::

   Context DBs can be concatenated using ``cat``. No need to rerun the builder!

**3.** Generate the target DB ``target_DB`` using ``foldseek createdb``. This prepares a FoldSeek DB from your folder containing the target set of protein structures (``input``). You probably don't have thousands of protein structures lingering around, so for large-scale analyses, you will need FoldSeek's builtin support of **ProstT5**, a LLM that directly translates amino acid sequences into FoldSeek's internal 3Di alphabet, skipping protein model prediction. **This is the key step that makes searching sequence databases using structural similarity computationally tractable.**

First make sure you have downloaded ProstT5's weights.

.. code-block:: bash

	foldseek databases ProstT5 <path-to-prostt5-weights> tmp

Then you can build the FoldSeek DB directly from your protein sequences in ``input``.

.. code-block:: bash

	foldseek createdb input target_DB --prostt5-model <path-to-prostt5-weights>

.. warning::

   This is still a **time- and computationally demanding** task! Consider using a GPU (by adding the ``--gpu 1`` flag if you have the hardware configured).

.. tip::

   FoldSeek offers a command to **merge existing DBs**: ``foldseek concatdbs``. Use it to concat existing target DBs, and save time and computational work.
   
   **You need to run it thrice, once for every database type (recognisable by the subscripts; none, ``_h``, and ``_ss``.)**

The search
^^^^^^^^^^
Run ``cfoldseeker`` in local mode using FoldSeek DB ``target_DB``, and context DB ``cds_db.gz``. Save results in a new folder ``results``.

.. code-block:: bash

	cfoldseeker -m local -q query_models -o results -ldb target_DB/target_DB -cdb cds_db.gz

Local-clustered mode
~~~~~~~~~~~~~~~~~~~~
The local-clustered mode requires a set of query protein structures, a genomic context table made using ``cfoldseeker-cds``, a MMseqs2 clustering TSV file, and a FoldSeek target database of the representative proteins.

Prior work
^^^^^^^^^^
**1.** Get structural models as CIF files of your query proteins, either experimentally or computationally (AlphaFold, ESMFold, OpenFold...). Collect them all in one folder ``query_models``.

**2.** Build a genomic context table using our provided helper tool ``cfoldseeker-cds``, which builds a context table directly from a set of NCBI or Bakta Genbank files, or from a folder holding an NCBI package of Genbank files. *(See the prior work section of local search above)*

**3.** Cluster your target sequences in the folder ``input_all`` using ``mmseqs``. You can do this using ``easy-cluster``, or ``easy-linclust`` for huge sequence databases. We recommend to use an identity and a coverage threshold of 90 % to ensure all proteins in a protein cluster have identical functions. **Use other thresholds at your own risk!**

.. code-block:: bash

	mmseqs easy-linclust input_all clustered tmp --min-seq-id 0.9 -c 0.9

This will, among others, produce a fasta file ``clustered_rep_seq.fasta`` containing the amino acid sequences of the representative protein of each cluster, and a clustering table ``clustered_table.tsv`` outlining the members and representatives of each sequence cluster.

**4.** Generate the target structure DB from the representative sequences using FoldSeek and ProstT5 (see also the prior work section of local search above). Make sure you have downloaded ProstT5's weights.

.. code-block:: bash

	foldseek createdb clustered_rep_seq.fasta target_DB --prostt5-model <path-to-prostt5-weights>

The search
^^^^^^^^^^
Run ``cfoldseeker`` in local_clustered mode using FoldSeek DB ``target_DB``, context DB ``cds_db.gz``, and preclustering table ``clustered_table.tsv``. Save results in a new folder ``results``.

.. code-block:: bash

	cfoldseeker -m local_clustered -q query_models -o results -ldb target_DB/target_DB -cdb cds_db.gz -scl clustered_table.tsv

Specifying search options
~~~~~~~~~~~~~~~~~~~~~~~~~
``cfoldseeker`` offers several filtering thresholds to refine your hit set and find relevant gene clusters.

General options
^^^^^^^^^^^^^^^
The general search options are a mix of what ``cblaster`` and ``foldseek`` offer. The available options are listed below.

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Sometimes, a protein may match with multiple query proteins, for example when you have two paralogs among your query proteins. This makes it tricky to determine what the query layout of an identified cluster is. For example, if two proteins in a cluster both match with query proteins 1 and 2, cluster layout *12* is equally correct as layout *21*. By default, if ``cfoldseeker`` encounters an identical cluster passing the filtering thresholds with different layouts, it will keep the one with the highest cluster score.

If you are interested in all possible cluster layouts passing your filtering thresholds rather than only the highest-scoring one, turn on the ``all-layouts`` flag to keep all passing configurations.

Using taxon filters
^^^^^^^^^^^^^^^^^^^^
The Foldseek webserver offers an interface to filter hits by a taxonomic filter. ``cfoldseeker`` exposes this interface in its remote mode via the ``-tf`` or ``--taxon-filter`` flag. Foldseek expects NCBI taxon IDs for its taxonomic filter.

If you are not sure what the exact taxon ID of your taxa group is, you can check it at `the NCBI Taxonomy website <https://www.ncbi.nlm.nih.gov/datasets/taxonomy/browser/>`_. Another option is to make the Foldseek webserver look it up for you. Browse to the webserver, and add your taxon as a taxonomic filter under *settings* at the website. Click then the *API* button on the top right, which will give you a pop-up with code to submit your query via the terminal. This code will include the taxonomy ID of your taxa group among the outlined settings.

Specifying outputs
~~~~~~~~~~~~~~~~~~~
``cfoldseeker`` can produce several outputs. By default, it produces only hit and cluster overview tables in TSV format, but several cblaster-style outputs and the raw Foldseek hit tables can be returned on request as well.

Cluster table
^^^^^^^^^^^^^^
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
^^^^^^^^^
The ``hits.tsv`` file gathers metadata about all hits part of the identified clusters. It contains the 16 columns described below.

+-----------------+----------------------------------------------------------------------------+
| **Column**      | **Description**                                                            |
+-----------------+----------------------------------------------------------------------------+
| db_id           | Unique hit ID (e.g. UniProt accession code)                                |
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
| crossref_method | Cross-referencing method (*KEGG*, *GenPept*, *WGS-GenPept*; or *local*)    |
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
| seqid           | Sequence identity between hit and query protein (in %)                     |
+-----------------+----------------------------------------------------------------------------+
| qcov            | Query coverage of the hit (which fraction of the query fits) (in %)        |
+-----------------+----------------------------------------------------------------------------+
| tcov            | Target coverage of the hit (which fraction of the target fits) (in %)      |
+-----------------+----------------------------------------------------------------------------+

Foldseek output
^^^^^^^^^^^^^^^
``cfoldseeker`` can return the raw Foldseek output from which it started. This can be useful if you want to track down why a certain hit was not found being part of a cluster, or how many unfiltered hits had been found for a certain query protein. In remote mode, ``cfoldseeker`` returns the json files it received from the Foldseek webserver. In the local modes, it will return the tab-separated text file returned by the local ``foldseek`` call.

*cblaster* outputs
^^^^^^^^^^^^^^^^^^
``cfoldseeker`` has tightly integrated ``cblaster``. All results are cast into a cblaster session, from which familiar outputs can be obtained, such as the summary table, the binary table, the hit plot, and the clinker alignment. See `the cblaster documentation <https://cblaster.readthedocs.io/en/latest/guide/search_module.html#specifying-output>`_ for specifics on these outputs.

Constructing context DBs
------------------------
Context DBs are crucial as they hold the genomic location of each protein's CDS, and hence facilitate identifying gene clusters. ``cfoldseeker`` context DBs are (compressed) headerless tab-separated text files constructed using ``cfoldseeker-cds``. They contain the columns described below.

+--------------+-----------------------------------------------------------------------+
| **Column**   | **Description**                                                       |
+--------------+-----------------------------------------------------------------------+
| Protein ID   | Unique protein identifier (e.g. NCBI Protein ID)                      |
+--------------+-----------------------------------------------------------------------+
| Description  | Description or annotation of the protein                              |
+--------------+-----------------------------------------------------------------------+
| Contig ID    | Unique identifier of the contig/scaffold harbouring the protein's CDS |
+--------------+-----------------------------------------------------------------------+
| Strand       | The strand coding for the protein                                     |
+--------------+-----------------------------------------------------------------------+
| Location     | Genomic coordinates of the CDS (format: "start..end")                 |
+--------------+-----------------------------------------------------------------------+
| Taxon ID     | Identifier of the taxon (e.g. NCBI taxon ID)                          |
+--------------+-----------------------------------------------------------------------+
| Taxon name   | Name of the taxon (identical to filelabel at default settings)        |
+--------------+-----------------------------------------------------------------------+
| Filelabel    | Filename of the associated sequence file                              |
+--------------+-----------------------------------------------------------------------+

``cfoldseeker-cds`` provides four parsing modes to construct this genomic context DB, depending on the source of your Genbank files and the annotation source or flexibility you prefer.

Typical ``cfoldseeker-cds`` commands look like the one below.

.. code-block:: bash

	cfoldseeker-cds -i <path-to-input-folder> -m <parsing-mode> -o <context-DB> (-gz)


NCBI Genbank package
~~~~~~~~~~~~~~~~~~~~~~~
This mode parses a package of NCBI Genbank files. After downloading and extracting a set of NCBI Genbank files from the Datasets portal, you typically have a folder ``ncbi_dataset`` with the following file structure.

.. code-block:: text

	ncbi_dataset/
    ├── ncbi_dataset/
    │   └── data/
    │       ├── <accession1>/
    │       |    └── genomic.gbff
    │       ├── <accession2>/
    │       |    └── genomic.gbff
    ...     ...
    │       ├── <accessionN>/
    │       |    └── genomic.gbff
    │       ├── assembly_data_report.jsonl
    │       ├── dataset_catalog.json
    │       └── data_summary.tsv
    ├── md5sum.txt
    └── README.md

When using this parsing mode, ``cfoldseeker-cds`` expects the input path to point to the parent folder ``ncbi_dataset``. It will then parse all Genbank files inside ``ncbi_dataset/ncbi_dataset/data``, and by default will use the name of the accession subfolders (typically an accession ID) as taxon names for the CDSes in that Genbank file. It will also keep track of every NCBI taxon ID it finds in the Genbank files. You can override the taxon names by making it fetch the taxon names associated with the taxon IDs using NCBI Entrez (``-tna``). However, fetching from Entrez may be cumbersome when the number of files is large. Therefore, we also provide the option to use a local rename file (``-tnf``). This rename file is a simple headerless tab-separated text file, with as first column the current taxon name (typically an accession ID), and as second column the new taxon name.

I usually prepare my rename files from the metadata table downloaded from that same NCBI Datasets portal I downloaded my Genbank files from (typically ``ncbi_dataset.tsv``). Download the metadata table for the accessions you've selected earlier for your Genbank files (download dropdown menu). This TSV table contains all information you need for a rename file. You only need the `Assembly accession` and the `Organism Name` columns. For bacteria, the `Organism Infraspecific Names Strain` column may also be relevant.

So, for bacteria, I usually cut-and-paste my rename file (`taxon_name_mapping`) together using the following bash command. The first column is for `Assembly accession`, and the second is for the space-joined columns `Organism Name` and `Organism Infraspecific Names Strain`.

.. code-block:: bash

	paste \
	<(cut -f 2 ncbi_dataset.tsv) \
	<(paste -d ' ' <(cut -f 4 ncbi_dataset.tsv) <(cut -f 6 ncbi_dataset.tsv)) | \
	tail -n +2 > taxon_name_mapping


NCBI Genbank files
~~~~~~~~~~~~~~~~~~~
This mode parses a folder of Genbank files assuming that they are standardised NCBI files. By default, this mode sources the taxon name from the filename. However, you can again override this by making ``cfoldseeker-cds`` fetch taxon names for the NCBI taxon IDs it found in the Genbank files using Entrez, or by supplying a rename file.

Bakta Genbank files
~~~~~~~~~~~~~~~~~~~~
This mode parses a folder of Genbank files you have generated yourself using a genome annotation tool like Bakta. However, since these files typically don't contain a taxon ID, ``cfoldseeker-cds`` will use a generic taxon ID. Again, the taxon name is sourced from the filename by default, but this is overridable. If you select automatic override, ``cfoldseeker-cds`` will derive a generic taxon name from the earlier generated taxon ID. Using a rename file is more preferable here.

Manually annotated TSVs
~~~~~~~~~~~~~~~~~~~~~~~
This mode offers maximum flexibility to annotate your CDSes. You can use multiple TSVes. Each file at least requires the columns specified below. **Make sure to include the header.**

+--------------+------------------------------------------------------+
| **Column**   | **Description**                                      |
+--------------+------------------------------------------------------+
| `gene_tag`   | Unique tag for the gene/CDS.                         |
+--------------+------------------------------------------------------+
| `name`       | Description or annotation of the gene.               |
+--------------+------------------------------------------------------+
| `contig`     | Unique identifier of the contig harbouring the gene. |
+--------------+------------------------------------------------------+
| `start`      | Start nucleotide coordinate of the CDS.              |
+--------------+------------------------------------------------------+
| `end`        | End nucleotide coordinate of the CDS.                |
+--------------+------------------------------------------------------+
| `strand`     | The strand coding for the protein.                   |
+--------------+------------------------------------------------------+
| `taxon_id`   | Identifier of the taxon.                             |
+--------------+------------------------------------------------------+
| `filelabel`  | Filename of the associated sequence file.            |
+--------------+------------------------------------------------------+

Extracting gene clusters
-------------------------
To facilitate downstream analyses, ``cfoldseeker-seqs`` facilitates getting sequence files for all identified gene clusters. Therefore, it fetches the nucleotide and amino acid sequences of every hit involved from your Genbank files. ``cfoldseeker-seqs`` supports various filters to pinpoint the the gene clusters you prefer (e.g. cluster numbers, organism filters, scaffold filters, score threshold, the top x scoring hits).

A gonna-catch-them-all ``cfoldseeker-seqs`` commands looks like the one below.

.. code-block:: bash

    cfoldseeker-seqs -s <path-to-session> -o <output-folder> -gb <genbanks-folder>

To extract clusters 1-10 and cluster 25 use

.. code-block:: bash

    cfoldseeker-seqs -s <path-to-session> -o <output-folder> -gb <genbanks-folder> --cluster-numbers 1-10 25

To extract all clusters with a score above 1000:

.. code-block:: bash

    cfoldseeker-seqs -s <path-to-session> -o <output-folder> -gb <genbanks-folder> --score-threshold 1000

To extract clusters only from specific organisms (regular expressions):

.. code-block:: bash

    cfoldseeker-seqs -s <path-to-session> -o <output-folder> -gb <genbanks-folder> --organisms "Aspergillus.*"

To extract clusters only from a specific range on scaffold_123 and all clusters on scaffold_234 (**NOTE: expects unique scaffold names**):

.. code-block:: bash

    cfoldseeker-seqs -s <path-to-session> -o <output-folder> -gb <genbanks-folder> --scaffolds scaffold_123:1-80000 scaffold_234


