## Single cell RNA-seq analysis with the Seurat package.
##
## This script performs the following:
##   1. Imports 10xGenomics scRNA-seq data
##   2. Filters cells with more than 20% mitochondrial gene expression
##   3. Runs a PCA on the most variable genes
##   4. Runs a tSNE analysis
##   5. Plots the tSNE result in a pdf
##
## Code is based off of the following tutorial:
##   https://satijalab.org/seurat/pbmc3k_tutorial_v3.html
##
## The raw 10xGenomics data was obtained from:
##   https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
##
## This is the selection test for a Google Summer of Code 2019 project proposed
## by the Canadian Center for Computational Genomics. Further details can be
## found at: https://bitbucket.org/mugqic/gsoc_2019
##
## Sean Nesdoly 2019-03-21
## MEng student, McGill University
## Montreal, QC, Canada

## Load packages ---------------------------------------------------------------
library(tidyverse)
library(Seurat)

## Import 10xGenomics scRNA-seq data -------------------------------------------
raw_counts <- Read10X(data.dir = "./filtered_gene_bc_matrices/hg19/")

# Initialize Seurat object with raw 10xGenomics RNA-seq data
seurat <- CreateSeuratObject(raw.data = raw_counts, # gene (rows) counts per cell (columns)
							 project = "10X_PBMC")  # project name

## Filter cells with more than 20% mitochondrial gene expression ---------------
# One indication for cellular stress or damage are cells that have a large
# proportion of reads mapping to mitochondrial genes.

# Obtain list of mitochondrial genes by subsetting for those starting with 'MT'
mito_genes <- grep("^MT-", rownames(seurat@data), value = TRUE)

# For each cell, compute fraction of expression that is by mitochondrial genes
fraction_mito <- Matrix::colSums(seurat@data[mito_genes, ]) / Matrix::colSums(seurat@data)

# Add to existing Seurat object as a metadata slot
seurat <- AddMetaData(seurat, fraction_mito, "fraction_mito_expr")

# Remove cells that have >20% mitochondrial gene expression
min_mito_expr <- 0
max_mito_expr <- 0.20
seurat <- FilterCells(seurat,
					  subset.names = "fraction_mito_expr", # subset by mitochondrial gene expression
					  low.thresholds = min_mito_expr,      # minimum threshold
					  high.thresholds = max_mito_expr)     # maximum threshold

## Principal Component Analysis (PCA) of the most variable genes ---------------
# Scale and log-transform gene expression for each cell
seurat <- NormalizeData(object = seurat,
						normalization.method = "LogNormalize",
						scale.factor = 10000)

# Find significantly variable genes across cells using default thresholds
seurat <- FindVariableGenes(object = seurat,
							do.plot = F)

# Regress out variation from coverage and mitochondrial gene expression
seurat <- ScaleData(object = seurat,
					vars.to.regress = c("nUMI", "fraction_mito_expr"))

# Run PCA on the most variable genes
seurat <- RunPCA(seurat,
				 pc.genes = seurat@var.genes, # only use the most variable genes
				 do.print = FALSE)

# Project the PCA (computed on the most variable genes) onto all genes
seurat <- ProjectPCA(seurat, do.print = FALSE)

# Visually select the number of principal components to use in downstream analyses
PCElbowPlot(seurat)
num_PCs <- 9 # here, I selected the top 9 PCs as being the most informative

## Cluster cells using only the most informative PCs ---------------------------
seurat <- FindClusters(seurat, 
					   reduction.type = "pca", 
					   dims.use = 1:num_PCs,   # only use informative PCs
					   print.output = FALSE, 
					   save.SNN = TRUE
)

## Run a tSNE analysis ---------------------------------------------------------
seurat <- RunTSNE(object = seurat, 
				  dims.use = 1:num_PCs) # only use informative PCs

## Plot the tSNE result in a pdf and save to file ------------------------------
pdf(file = "tSNE.pdf")
DimPlot(object = seurat, reduction = "tsne", 
		plot.title = paste(seurat@project.name,
						   "tSNE using the top",
						   num_PCs,
						   "principal components."))
dev.off()