What does cfoldseeker do?
===============================

``cfoldseeker`` finds clusters of proteins with colocalised coding sequences.

Starting from a set of protein structures of interest, it searches protein structure databases for structurally homologous proteins of which the coding sequences are clustered together in the genome. This tool is developed as the protein structure similarity-based sister tool of ``cblaster``, and aims to bring the more sensitive protein structure similarity to gene cluster mining, aiding in finding coordinated groups of proteins and/or protein domains.

``cfoldseeker`` has several useful features:

- **The remote mode** facilitates querying the AlphaFoldDB and fetching genomic contexts via various APIs in parallel (KEGG, UniProt, ENA GenPept).
- **The local mode** allows querying a local protein structure database using a premade genomic context table constructed from GFF files or TSV files with ``cfoldseeker-cds``.
- **The local_clustered** enables querying a database of representative protein structures generated from a preclustered sequence database, as a means to reduce the computational workload of the protein model generation. If the representative protein of a sequence cluster is picked up, all members of that sequence cluster are added to the hit set.
- Tight integration with ``cblaster`` by handling its session files, allowing the use of existing processing and visualisation workflows such as ``clinker`` and ``CAGEcleaner``.

For detailed information on how to use the several modes and helper tools, and how to make it work in conjunction with other relevant tools, head over to the user guide and the tutorial.
