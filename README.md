[![DOI](https://img.shields.io/badge/STAR%20Protocols-10.1016%2Fj.xpro.2022.101490-%23007DBC)](https://star-protocols.cell.com/protocols/1750)

[![DOI](https://img.shields.io/badge/Molecular%20Cell-10.1016%2Fj.molcel.2020.11.046-%23007DBC)](https://www.cell.com/molecular-cell/fulltext/S1097-2765(20)30888-1)

[![DOI](https://zenodo.org/badge/461916697.svg)](https://zenodo.org/badge/latestdoi/461916697)

# High-throughput sequencing SELEX for the determination of DNA-binding protein specificities *in vitro*

> Pantier, R., Chhatbar, K., Alston, G., Lee, H.Y. and Bird, A., 2022. High-throughput sequencing SELEX for the determination of DNA-binding protein specificities in vitro. [*STAR protocols, 3*(3), p.101490.](https://star-protocols.cell.com/protocols/1750)

> Pantier, R., Chhatbar, K., Quante, T., Skourti-Stathaki, K., Cholewa-Waclaw, J., Alston, G., Alexander-Howden, B., Lee, H.Y., Cook, A.G., Spruijt, C.G. and Vermeulen, M., 2021. SALL4 controls cell fate in response to DNA base composition. [*Molecular cell, 81*(4), pp.845-858.](https://www.cell.com/molecular-cell/fulltext/S1097-2765(20)30888-1)

High-throughput sequencing SELEX (HT-SELEX) is a powerful technique for unbiased determination of preferred target motifs of DNA-binding proteins in vitro. The procedure depends upon selection of DNA binding sites from a random library of oligonucleotides by purifying protein-DNA complexes and amplifying bound DNA using the polymerase chain reaction. Here, we describe an optimized step-by-step protocol for HT-SELEX compatible with Illumina sequencing. We also introduce a bioinformatic pipeline (eme_selex) facilitating the detection of promiscuous DNA binding by analyzing the enrichment of all possible k-mers.

## eme_selex

eme_selex (***Every Motif Ever for SELEX Analysis***) is a Python package to perform k-mer abundance analysis in DNA sequences. ***eme_selex*** is developed to perform fast and efficient analysis of short k-mers (tested with k-mers up to length 10). 

While ***eme_selex*** can be used for general purpose k-mer analysis, motivation to develop ***eme_selex*** is to perform [**S**ystemic **E**volution of **L**igands by **EX**ponential enrichment coupled with **H**igh **T**hroughput sequencing (HT-SELEX)](https://en.wikipedia.org/wiki/Systematic_evolution_of_ligands_by_exponential_enrichment) analysis in a Pythonic way. By default, for every k-mer, ***eme_selex*** quantifies the fraction of reads containing that k-mer *in a non-redundant manner*. After the quantification, a basic position frequency matrix (PFM) for the top 50 k-mers is generated. If the user wants to generate more PFMs, they can change the *top* keyword argument to a desired number.

### Installation

```bash
pip install eme_selex
```

## Tutorial for HT-SELEX analysis

Jupyter notebooks detailing the usage of ***eme_selex*** and extensive analysis for HT-SELEX are hosted here [https://eme-selex.readthedocs.io](https://eme-selex.readthedocs.io)
