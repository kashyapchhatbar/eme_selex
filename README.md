[![DOI](https://zenodo.org/badge/461916697.svg)](https://zenodo.org/badge/latestdoi/461916697)

# eme_selex

eme_selex (***Every Motif Ever for SELEX Analysis***) is a Python package to perform k-mer abundance analysis in DNA sequences. ***eme_selex*** is developed to perform fast and efficient analysis of short k-mers (tested with k-mers up to length 10). 

While ***eme_selex*** can be used for general purpose k-mer analysis, motivation to develop ***eme_selex*** is to perform [**S**ystemic **E**volution of **L**igands by **EX**ponential enrichment coupled with **H**igh **T**hroughput sequencing (HT-SELEX)](https://en.wikipedia.org/wiki/Systematic_evolution_of_ligands_by_exponential_enrichment) analysis in a Pythonic way. By default, for every k-mer, ***eme_selex*** quantifies the fraction of reads containing that k-mer *in a non-redundant manner*. After the quantification, a basic position frequency matrix (PFM) for the top 50 k-mers is generated. If the user wants to generate more PFMs, they can change the *top* keyword argument to a desired number.

## Installation

```bash
pip install eme_selex
```

## Tutorial for HT-SELEX analysis

Jupyter notebooks detailing the usage of ***eme_selex*** and extensive analysis for HT-SELEX are hosted here [https://eme_selex.readthedocs.io](https://eme_selex.readthedocs.io)
