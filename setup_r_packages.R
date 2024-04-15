
install.packages("optparse")
install.packages("Seurat")
install.packages("jsonlite")
install.packages("BiocManager")
install.packages("remotes")

library("BiocManager")
library("remotes")

BiocManager::install("SingleCellExperiment")
BiocManager::install("SC3")

remotes::install_github("mojaveazure/seurat-disk")


library("optparse")
library("Seurat")
library("jsonlite")
library("BiocManager")
library("remotes")
library("SingleCellExperiment")
library("SC3")
library("SeuratDisk")