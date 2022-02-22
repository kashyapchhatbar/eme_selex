---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Bioinformatics Workflow

## Before you begin

- Please install the following software in order to replicate the analysis.

    - [`flexbar`](https://github.com/seqan/flexbar) for quality trimming of sequencing reads
    - [`snakemake`](https://snakemake.readthedocs.io/en/stable/) for workflow management
    - [`eme_selex`](https://github.com/kashyapchhatbar/every-motif-ever) for kmer abundance
    - [`jupyterlab`](https://jupyter.org/install) for web-based interactive development environment 
    - [`tqdm`](https://tqdm.github.io/) for displaying progress
    - [`pandas`](https://pandas.pydata.org/) for data frame creation and manipulation
    - [`seaborn`](https://seaborn.pydata.org/) for data visualisation
    - [`plotly`](https://plotly.com/python/getting-started/) for data visualisation
    - [`logomaker`](https://logomaker.readthedocs.io/) for making motif logos
    - [`upsetplot`](https://upsetplot.readthedocs.io/) for visualising set overlaps

    ```{tip} 
    Use [anaconda](https://anaconda.org/) to create a separate environment and efficient installation. The following command will install all the relevant packages
    ```

    ```bash
    conda create -n eme_selex -c bioconda flexbar snakemake pip
        jupyterlab tqdm pandas seaborn

    conda activate eme_selex

    pip install eme_selex logomaker upsetplot

    conda install -c plotly plotly==5.6.0
    ```

## Preprocessing of fastq files

### Quality trim sequencing reads

- Remove low quality sequencing reads and trim them down to the number of bases in your random library. In our case, that number is 20. Use [`flexbar`](https://github.com/seqan/flexbar) with the following arguments below.

    ```bash
    flexbar --reads {input} --post-trim-length 20 --min-read-length 20 
      --qtrim-threshold 30 --output-reads {output} --fasta-output --number-tags 
      --stdout-log > {log}
    ```

    ```{tip} 
    Use a workflow manager like [`snakemake`](https://snakemake.readthedocs.io/en/stable/) to automate this step for all samples.
    ```
### Automate quality trimming using a workflow manager

- Using [`snakemake`](https://snakemake.readthedocs.io/en/stable/), automate preprocessing of all samples. Create `Snakefile` in the parent directory, all the samples should be present in `fastq` sub-directory.

    ```python
    samples = [f"RV{s:02d}" for s in range(1, 40)]

    rule preprocess:
        input: expand("fasta/{sample}.fasta.gz", sample=samples)

    rule flexbar:
        input: "fastq/{sample}.fastq.gz"
        output: "fasta/{sample}.fasta.gz"
        log: "logs/flexbar.{sample}.log"
        shell: "flexbar --reads {input} --post-trim-length 20 --min-read-length 20 
                  --qtrim-threshold 30 --output-reads {output} --fasta-output --number-tags 
                  --stdout-log > {log}"
    ```

- Execute the rule `flexbar` using the following command below.

    ```bash
    snakemake flexbar --cores 1
    ```

    ```{tip}
    Increase the number of simultaneous executions by changing the parameter ``--cores``
    ```
(kmer_fraction_from_file)=
## Calculate k-mer abundance

- Initialize [`defaultdict`](https://docs.python.org/3/library/collections.html#collections.defaultdict) to store the data.

    ```python
    from collections import defaultdict
    counts, fractions, models = defaultdict(dict), defaultdict(dict), defaultdict(dict)
    ```

- Calculate kmer abundance, fraction and top position frequency matrix (PFM) models using [`eme_selex`](https://github.com/kashyapchhatbar/every-motif-ever). For example, to calculate fraction of 5-mers, set `k=5` and execute the function over all samples.

    ```python
    from eme_selex.eme import kmer_fraction_from_file as kf
    from tqdm.auto import tqdm

    samples = [f"RV{s:02d}" for s in range(1, 40)]

    k = 5

    for sample in tqdm(df["SampleName"].values, leave=False):
        c, f, m = kf(f"fasta/{sample}.fasta.gz", k=k)
        counts[sample] = c
        fractions[sample] = f
        models[sample] = m
    ```

### Saving the results

- Store the results in binary format for easy access and comparative analysis later.

    ```python
    import pickle
    import bz2

    with bz2.open(f"results/results_{k}.bz2", "wb") as wH:
        pickle.dump(counts, wH)
        pickle.dump(fractions, wH)
        pickle.dump(models, wH)
    ```
