Tutorial
========

This tutorial will guide you through the practicalities of running a ``cfoldseeker`` search, including all prior work with other tools. This is an example running the `local_clustered` mode, the mode requiring the most prior work. This is probably the mode that can do the most comprehensive analyses if taken good care of the target database.

We will download GFF and protein Fasta files from NCBI, build a genomic context DB using ``cfoldseeker-cds``, generate ProstT5 encodings using ``foldseek`` for the representative proteins of the downloaded sequence database clustered with ``MMseqs2``, and finally run ``cfoldseeker``.

.. note::

   All tools mentioned in the tutorial below (except the NCBI Datasets CLI) are present in the ``cfoldseeker`` conda environment. You can just run everything inside this environment.

.. note::

   You will need **GPU acceleration** for building the target database in this tutorial!

Searching thioalbimide BGC against *Bacillaceae*
-------------------------------------------------

In this example, you will search the thioalbamide biosynthetic gene cluster (BGC) from `Amycolatopsis alba DSM 44262` (assembly accession GCF_000384215.1) against all *Bacillaceae* proteomes in NCBI.

Preparing query proteins
~~~~~~~~~~~~~~~~~~~~~~~~

First, get yourself the amino acid sequences of your query proteins. I have already done so for this BGC: you can find them in the `thioalbamide.faa` fasta file in the example output of the GitHub repo. I predicted the protein structures for each sequence using `the AlphaFold3 webserver <https://alphafoldserver.com/welcome>`_. For each protein, extract the top-ranked structure (model 0) from the downloaded results zip and gather them all in a `query_models` folder (mine is also in the repo).

Preparing context DB
~~~~~~~~~~~~~~~~~~~~

To build the context DB, we need all `Bacillaceae` GFF files from NCBI. I usually first get a list of accession IDs from the NCBI website. In this case, I searched for `Bacillaceae` (taxonomy ID 186817) and started browsing their genomes. I applied some gentle filtering (RefSeq-annotated genomes, excluding atypical ones). As of 17th April 2026, there were 16.394 genomes.

Now, to download the GFF files, select all genomes, and start downloading a package (download dropdown menu). Select only RefSeq as file source and make sure only a GFF file is selected as file type.

Alternatively, you can also download these files using `the NCBI Datasets CLI <https://github.com/ncbi/datasets>`_. This is quicker and more robust than downloading in your browser. First, from the website, download a table (same download dropdown menu), copy-paste the column with the RefSeq accession IDs (GCF_*) in a notepad program and save it as a new text file `accessions.txt`. Then fire up a terminal and run the following commands in the folder containing `accessions.txt`.

.. code-block:: bash

	datasets download genome accession --inputfile accessions.txt --include gff3 --dehydrated
	unzip ncbi_dataset.zip && rm ncbi_dataset.zip
	datasets rehydrate --directory ncbi_dataset
	mv ncbi_dataset gff_package

This will get you a folder `ncbi_dataset` holding a package of NCBI GFF files.

Finally, run ``cfoldseeker-cds`` to construct the genomic context DB in compressed form (`bacillaceae_cds.gz`).

.. code-block:: bash

	cfoldseeker-cds -i gff_package -m ncbi-package -o bacillaceae_cds.gz -gz

Preparing target DB
~~~~~~~~~~~~~~~~~~~

To build the target DB, we need all `Bacillaceae` protein Fasta files from NCBI. Downloading these can be done similarly as for the GFFs for the context DB, yet do not forget to check Protein Fasta.

Via the NCBI Datasets CLI, you can reuse your earlier made `accessions.txt` using the following commands.

.. code-block:: bash

	datasets download genome accession --inputfile accessions.txt --include protein --dehydrated
	unzip ncbi_dataset.zip && rm ncbi_dataset.zip
	datasets rehydrate --directory ncbi_dataset
	mv ncbi_dataset faa_package

To make my life easier, I collect all the protein fasta files in this NCBI package into one new folder `faas` using this bash oneliner.

.. code-block:: bash

	mkdir faas
	dir -1 faa_package/ncbi_dataset/data | xargs -I % mv faa_package/ncbi_dataset/data/%/protein.faa faas/%.faa

Together, these files may easily contain more than 40M protein sequences. So, to reduce later computational work spent generating protein models, we cluster them first using ``mmseqs easy-linclust`` at an identity and coverage threshold of 90 %. Using 32 cores on a HPC, this took about 15 minutes, resulting in 5.157.432 clusters, or **an eightfold reduction** in proteins for which we need to generate protein models! 

.. code-block:: bash

	mmseqs easy-linclust faas/* clustered tmp --min-seq-id 0.9 -c 0.9

``MMseqs2`` will produce three files: a fasta file with all sequences (`clustered_all_seqs.fasta`), one with only the representative sequences (`clustered_rep_seq.fasta`), and a clustering table in TSV format (`clustered_cluster.tsv`). The latter one is the one ``cfoldseeker`` will need later on, while the second one is the input for the protein model generation.

.. tip::

   Although you can definitely run ``cfoldseeker`` against a set of protein structure models, it is currently computationally intractable to get full protein structures at the same scale as we do for sequences in the NCBI databases.

   `ProstT5 <https://academic.oup.com/nargab/article/6/4/lqae150/7901286>`_ is a LLM that mitigates this by directly translating amino acid sequences to Foldseek's 3Di alphabet, skipping the expensive structure prediction step. ProstT5 is integrated in ``foldseek``.

We will now generate prostT5 3Di encodings for our 5M+ representative proteins using ``foldseek createdb``.

First make sure you have downloaded ProstT5's weights to a folder `weights` by running the code line below

.. code-block:: bash

	foldseek databases ProstT5 weights tmp

Then start generating 3Di encodings using GPU acceleration.

.. code-block:: bash

	mkdir DB	
	foldseek createdb clustered_rep_seq.fasta DB/Bacillaceae --gpu 1 --prostt5-model weights/

Using two NVIDIA H200 GPUs (Hopper generation) on an HPC, this took xx hours.

.. warning::

   **This is a computationally demanding task!** It is highly recommended to use GPU acceleration with an NVIDIA GPU of at least the Ampere generation!

   You can get GPU-compatible binaries `here <https://dev.mmseqs.com/foldseek/>`_ if there are no binaries compiled on your (HPC) system.

``foldseek`` will have generated 11 files in the folder `DB`, all starting with the prefix `Bacillaceae`. This is your local target structure DB.

Running ``cfoldseeker``
~~~~~~~~~~~~~~~~~~~~~~~

We now have all prerequisites to run ``cfoldseeker`` in local-clustered mode.

The following command runs ``cfoldseeker`` at default search settings using 14 cores, and makes it produce every supported output file in the new folder `cfoldseeker_search`. By appending a ``tee`` pipe, you can capture the logs in a log file.

.. code-block:: bash

	cfoldseeker \
	-m local_clustered \
	-c 14 \
	-q query_models \
	-o cfoldseeker_search \
	-ldb DB/Bacillaceae \
	-cdb bacillaceae_cds.gz \
	-scl clustered_cluster.tsv \
	--session --summary --binary --plot --clinker --foldseek | \
	tee cfoldseeker.log 

All intermediary and output files of this tutorial can also be found in the example output of the ``cfoldseeker`` GitHib repo. The large files (context DB, target DB) are only available in the Zenodo copy.
