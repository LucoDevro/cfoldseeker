# cfoldseeker

[![Docs](https://img.shields.io/readthedocs/cfoldseeker/latest?style=flat-square&maxAge=600&logo=readthedocs)](https://cfoldseeker.readthedocs.io/en/latest/)
[![Downloads](https://anaconda.org/bioconda/cfoldseeker/badges/downloads.svg)](https://bioconda.github.io/recipes/cfoldseeker/README.html#download-stats)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/cfoldseeker?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/cfoldseeker)
[![Docker Image Version](https://img.shields.io/docker/v/lucodevro/cfoldseeker?sort=semver&label=docker&logo=docker)](https://hub.docker.com/r/lucodevro/cfoldseeker)
[![PyPI version](https://img.shields.io/pypi/v/cfoldseeker?sort=semver&logo=pypi)](https://pypi.org/project/cfoldseeker/)

## Description
`cfoldseeker` finds homologous gene clusters via protein structural similarity. It searches structural homologs for your query protein structures using `foldseek` (both local and remote target databases supported) and identifies the genomically colocalised hits among these by fetching the genomic location of each protein's coding sequence (fetched from various remote cross-referencing APIs, or from a locally prepared database).

`cfoldseeker` has been designed as the structural similarity-driven sister tool of `cblaster`, which it tighly integrates for generating outputs. As such, `cfoldseeker` can naturally produce `cblaster`-style output and `clinker` visualisations.

> [!TIP]
> Although `cfoldseeker` can be used as a stand-alone tool, it is the structural similarity-based discovery engine of the ✨ [`csuite`](https://github.com/LucoDevro/csuite) ✨, our new integrated toolbox featuring streamlined workflows for both sequence and protein structure-based gene cluster mining. Try it out!

![workflow](workflow.png)

## Features

- **A remote search mode** for searches against the AlphaFoldDB, leveraging the [Foldseek webserver](https://search.foldseek.com) and various cross-referencing APIs for fetching genomic locations ([`kegg_pull`](https://github.com/MoseleyBioinformaticsLab/kegg_pull), [UniProt ID mapping](https://www.uniprot.org/id-mapping), [ENA Browser API](https://www.ebi.ac.uk/ena/browser/api/)).
- **A local search mode** for searches against a local protein structure DB prepared with [`foldseek`](https://github.com/steineggerlab/foldseek).
- **A local-clustered search mode** for searches against a local `foldseek` DB of representative proteins derived from a sequence set preclustered with [`MMseqs2`](https://github.com/soedinglab/MMseqs2)). If the representative protein of a sequence cluster is identified as a homolog, all cluster members are added to the hit set.
- **A helper tool** to construct local genomic context databases: `cfoldseeker-cds`
- **Tight integration with `cblaster`**, facilitating similar output and interactive `clinker` visualisations

## Installation, documentation and more
For installation instructions, usage, a tutorial and more, head over to the [`cfoldseeker` docs](https://cfoldseeker.readthedocs.io/en/latest)!

## Citations
If you found `cfoldseeker` useful, please cite our manuscript:

```
De Vrieze, L., Masschelein, J. (2026) In preparation
```

## License
`cfoldseeker` is freely available under an MIT license.

Use of the third-party software, libraries or code referred to in the References section above may be governed by separate terms and conditions or license provisions. Your use of the third-party software, libraries or code is subject to any such terms and you should check that you can comply with any applicable restrictions or terms and conditions before use.
