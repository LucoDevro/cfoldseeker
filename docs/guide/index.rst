What does cfoldseeker do?
===============================

``cfoldseeker`` finds clusters of proteins with colocalised coding sequences.

Starting from a set of protein structures of interest, it searches protein structure databases for structurally homologous proteins of which the coding sequences are clustered together in the genome. This tool is developed as the protein structure similarity-based sister tool of `cblaster <https://github.com/gamcil/cblaster>`_, and aims to bring more sensitive protein structure similarity to gene cluster mining, aiding in finding colocalised groups of proteins and/or protein domains.

``cfoldseeker`` has several useful features:

- **The remote mode** facilitates searching against the AlphaFoldDB and fetching genomic contexts via various APIs in parallel (KEGG, UniProt, ENA GenPept).
- **The local mode** allows searching against a local protein structure database using a premade genomic context table constructed with ``cfoldseeker-cds``.
- **The local_clustered mode** enables searching against a database of representative protein structures from a preclustered sequence database, as a means to reduce the computational workload of the protein model generation. If the representative protein of a sequence cluster is picked up, all members of its sequence cluster are added to the hit set.
- **Tight integration with ``cblaster``** by casting results into its session files, allowing the use of existing processing and visualisation workflows such as ``clinker`` and ``CAGEcleaner``.

For detailed information on how to use the several modes and helper tools, and how to make it work in conjunction with other relevant tools, head over to the user guide and the tutorial.
