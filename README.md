# gsoc19-single-cell-pipeline

>Sean Nesdoly 2019-03-16  
>MEng student, McGill University  
>Montreal, QC, Canada

`Seurat` single cell RNA-seq analysis launched with a single step `GenPipes`
pipeline.

This is implemented as the first step of the application process for a Google
Summer of Code 2019 project that has been proposed by the [Canadian Centre for
Computational
Genomics](https://bitbucket.org/mugqic/gsoc_2019#markdown-header-develop-genpipe-single-cell-pipeline).
The goal of this selection test is to write a single step `GenPipes` pipeline
that executes an R script to run an analysis on single cell RNA-seq data using
the [Seurat R package](https://satijalab.org/seurat/). Specifically, the R
script performs the following:

1.  Imports 10xGenomics scRNA data
2.  Filters cells with more than 20% mitochondrial genes
3.  Runs a PCA on the most variable genes
4.  Runs a tSNE analysis
5.  Plots the tSNE result in a pdf
