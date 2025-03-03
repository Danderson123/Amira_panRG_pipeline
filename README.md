# Amira_panRG_pipeline

## Overview

Here are some brief instructions to run the Amira panRG construction pipeline on a directory of reference assemblies. The pipeline is implemented in snakemake which is a workflow manager built on top of Python. The steps involved in the pipeline are defined by rules in the Snakefile in this directory. Snakemake automatically detects what files need to exist and what rules need to be run to create the output file that is specified by `rule all`. You will need to use a MacOS or Linux terminal to run this pipeline. The pipeline is intended to be used on a cluster computing system. It should work on a local machine but there are several rules that operate on thousands of files so it will take a long time to run if you do not have many CPUs available.

## Installation

You first need to clone and then cd into the repository by running:
```
git clone https://github.com/Danderson123/Amira_panRG_pipeline && cd Amira_panRG_pipeline
```
We use poetry to manage Amira_panRG_pipeline's python dependencies. You can setup a virtual conda environment and install poetry by running:
```
conda create -n amira_panRG_pipeline python=3.9 && conda activate amira_panRG_pipeline && pip install poetry
```
The python dependencies can then be installed by running:
```
poetry install
```
You also need to have [Conda](https://docs.anaconda.com/miniconda/miniconda-install/) and [MAFFT](https://mafft.cbrc.jp/alignment/software/source.html) in your PATH.

## Inputs

`E_coli_config.yaml` is the file telling the snakemake rules which parameters and file paths to use for the different tools in the pipeline.

You will need to add various files and databases to the subdirectories in `data` in order to use this pipeline. See below for details.
* `data/reference_assemblies`: This is where you put all of the reference assemblies for the species you want to build a panRG for. You will need at least 2 assemblies for the pipeline to work. The extension of these files must be `.fa`.
* `data/test_assemblies`: It is not essential to put anything in here. Gene annotation tools often can give the same gene multiple different names. Panaroo is designed to unify the names and make sure that the same gene is named the same thing across all of the input samples. The point of this directory is so that you can get the gene names of any test samples relative to the names output by Panaroo but the pipeline does not include the sequences of these samples in the final panRG. I used this feature to estimate the accuracy of Pandora's gene calling when you apply it to samples not in the panRG.
* `data/bakta_db`: This is the database needed by bakta to annotate you reference assemblies. I used the full database. It is quite a large download and it can be downloaded from here https://zenodo.org/records/10522951.
* `data/reference_genes`: This is where I have stored a FASTA of some additional genes I want Pandora to find (namely the AMR gene reference sequences and plasmid-specific genes). The FASTAs are included as a sample when running Panaroo so the sequences are included in the final panRG. If you want to supplement the panRG with additional sequences let me know, it is not too hard to do but it requires modifying the snakemake workflow.
* `data/poppunk_databases`: This is not essential. This directory is where I put the poppunk database for the species I am interested in. PopPUNK is a tool that clusters bacterial isolates into lineages. By default the pipeline skips running PopPUNK but if you choose not to skip it, then one sample per PopPUNK lineage will be chosen to include in the panRG. This is useful if you want to build a panRG from hundreds or thousands of samples. The database for *E. coli* can be downloaded from here https://ftp.ebi.ac.uk/pub/databases/pp_dbs/escherichia_coli_v2_refs.tar.bz2.

## Usage

To run the pipeline on your local machine use:
`snakemake --cores <NUMBER OF CPUs> --configfile E_coli_config.yaml --nolock --rerun-incomplete --use-conda`
