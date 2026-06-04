# Align and Trees bash script

A bash script designed to align multiple genes at a time, filter low sample size alignments, produce individual gene trees, and run a final consensus tree based on all individual gene trees.

## Pipeline Architecture

1. **Alignments:** MAFFT to align individual genes
2. **Filtering:** For loop to check number of samples in each gene alignment and filter out genes with few samples
3. **Gene Trees:** IQ-TREE 2 to produce individual gene trees
4. **Filtering 2:** nw_ed to collapse low support branches in gene trees before generating consensus tree
5. **Consensus tree:** ASTRAL to generate a consensus tree

## Quick Start

```
align_and_trees.sh -j path/to/astral.jar -d path/to/base_directory/ -samples min_number_of_samples -support min_support_value
```

## Dependencies
This pipeline is executed via a Bash script and requires a Unix-like environment (Linux, macOS, or WSL). Ensure the following bioinformatics tools and dependencies are installed and accessible within your system's $PATH prior to execution.



* MAFFT: Required for multiple sequence alignment.
* IQ-TREE 2: Required for maximum likelihood gene tree inference and bootstrapping.
* Newick Utilities (specifically nw_ed): Required to programmatically collapse low-support branches in the concatenated gene trees.
* ASTRAL: Required for coalescent-based species tree inference. Note: You must provide the absolute path to the ASTRAL .jar file using the -j flag when executing the pipeline.
* Java (v1.8 or higher): Required to execute the ASTRAL .jar file.



Installation Note: With the exception of Java and ASTRAL, all dependencies can be easily managed and installed using the Bioconda package manager (e.g., `conda install bioconda::mafft bioconda::iqtree bioconda::newick_utils`).

## Requirements
A `merged_loci` folder must be set up containing individual unaligned multi-fasta files for each gene (where each file contains the sequences for all available samples for that specific gene).

## Parameters

| Parameter | Description |
|-----------|-------------|
| `-j` or `--jar` | Absolute path to the ASTRAL .jar file. |
| `-d` or `--dir` | (Optional) Path to the base directory containing the `merged_loci` folder. Default = Current working directory. |
| `-samples` | (Optional) Minimum number of samples needed in gene alignment before being filtered out. Default = 4 |
| `-support` | (Optional) Minimum branch bootstrap support value needed to avoid collapsing of branches in gene trees before generating consensus tree. Default = 50|


## Pipeline Outputs

* **loci_aligned/:** Fasta alignments of each gene with more than minimum number of samples
* **loci_discarded/:** Fasta alignments of each gene with less than minimum number of samples
* **gene_trees/:** Gene tree outputs from IQ-TREE 2
* **gene_trees/all_gene_trees.newick:** Concatenation of all gene trees
* **gene_trees/all_gene_trees_collapsed.tre:** Concatenation of all gene trees with branches collapsed if low support value
* **gene_trees/astral_all_genes.tre:** Final consensus tree



## Citation

**Align and Trees** is open for all. If you use this pipeline in your research, please cite this GitHub repository directly:

> Batalla, M. I. (2026). *Align and Trees*: A Bash Script to Automatically Align Genes and Produce a Final Consensus Tree. https://github.com/martinbatalla/align_and_trees