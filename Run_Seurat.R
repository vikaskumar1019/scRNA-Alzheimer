library(Seurat)  
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S1_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S1_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S2_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S2_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S3_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S3_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S4_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S4_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S5_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S5_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S6_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S6_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S7_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S7_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S8_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S8_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S9_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S9_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S10_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S10_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S11_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S11_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S12_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S12_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S13_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S13_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S14_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S14_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S15_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S15_New.rds")

library(Seurat)
library("tidyseurat")
setwd('/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S16_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX')
solote_matrix <- ReadMtx("matrix.mtx","barcodes.tsv","features.tsv",feature.column=1)
solote_seuratobj <- CreateSeuratObject(count=solote_matrix,min.cells=3,project="SoloTE")
solote_seuratobj$percent_mt <- PercentageFeatureSet(solote_seuratobj,pattern="^chrM-|^ENSG00000198695|^ENSG00000198712|^ENSG00000198727|^ENSG00000198763|^ENSG00000198786|^ENSG00000198804|^ENSG00000198840|^ENSG00000198886|^ENSG00000198888|^ENSG00000198899|^ENSG00000198938|^ENSG00000209082|^ENSG00000210049|^ENSG00000210077|^ENSG00000210082|^ENSG00000210100|^ENSG00000210107|^ENSG00000210112|^ENSG00000210117|^ENSG00000210127|^ENSG00000210135|^ENSG00000210140|^ENSG00000210144|^ENSG00000210151|^ENSG00000210154|^ENSG00000210156|^ENSG00000210164|^ENSG00000210174|^ENSG00000210176|^ENSG00000210184|^ENSG00000210191|^ENSG00000210194|^ENSG00000210195|^ENSG00000210196|^ENSG00000211459|^ENSG00000212907|^ENSG00000228253",assay = "RNA")
seuratobj <- subset(solote_seuratobj,subset = nFeature_RNA>=200 & percent_mt<=5)
saveRDS(seuratobj,file="Final_S16_New.rds")

======

seurat_obj_A1 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S2_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S2_New.rds")
seurat_obj_A1$label <- "A"
seurat_obj_A1$replicate <- "Replicate1"
seurat_obj_A2 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S7_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S7_New.rds")
seurat_obj_A2$label <- "A"
seurat_obj_A2$replicate <- "Replicate2"
seurat_obj_A3 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S8_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S8_New.rds")
seurat_obj_A3$label <- "A"
seurat_obj_A3$replicate <- "Replicate3"
seurat_obj_A4 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S11_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S11_New.rds")
seurat_obj_A4$label <- "A"
seurat_obj_A4$replicate <- "Replicate4"
seurat_obj_A5 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S12_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S12_New.rds")
seurat_obj_A5$label <- "A"
seurat_obj_A5$replicate <- "Replicate5"
seurat_obj_A6 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S13_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S13_New.rds")
seurat_obj_A6$label <- "A"
seurat_obj_A6$replicate <- "Replicate6"
seurat_obj_A7 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S15_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S15_New.rds")
seurat_obj_A7$label <- "A"
seurat_obj_A7$replicate <- "Replicate7"
seurat_obj_A8 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S16_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S16_New.rds")
seurat_obj_A8$label <- "A"

seurat_obj_B1 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S1_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S1_New.rds")
seurat_obj_B1$label <- "B"
seurat_obj_B1$replicate <- "Replicate1"
seurat_obj_B2 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S3_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S3_New.rds")
seurat_obj_B2$label <- "B"
seurat_obj_B2$replicate <- "Replicate2"
seurat_obj_B3 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S4_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S4_New.rds")
seurat_obj_B3$label <- "B"
seurat_obj_B3$replicate <- "Replicate3"
seurat_obj_B4 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S5_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S5_New.rds")
seurat_obj_B4$label <- "B"
seurat_obj_B4$replicate <- "Replicate4"
seurat_obj_B5 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S6_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S6_New.rds")
seurat_obj_B5$label <- "B"
seurat_obj_B5$replicate <- "Replicate5"
seurat_obj_B6 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S9_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S9_New.rds")
seurat_obj_B6$label <- "B"
seurat_obj_B6$replicate <- "Replicate6"
seurat_obj_B7 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S10_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S10_New.rds")
seurat_obj_B7$label <- "B"
seurat_obj_B7$replicate <- "Replicate7"
seurat_obj_B8 <-readRDS("/home/vikku19/software/SoloTE-main/New_Run1/GSE163577_S14_Out/Alz1_SoloTE_output/Alz1_locustes_MATRIX/Final_S14_New.rds")
seurat_obj_B8$label <- "B"
seurat_obj_B8$replicate <- "Replicate8"

#With new method https://satijalab.org/seurat/articles/integration_introduction

combined_seurat_obj <- merge(seurat_obj_A1, y = list(seurat_obj_A2,seurat_obj_A3,seurat_obj_A4,seurat_obj_B1,seurat_obj_B2,seurat_obj_B3,seurat_obj_B4), add.cell.ids = c("A_Rep1", "A_Rep2","A_Rep3","A_Rep4","B_Rep1", "B_Rep2","B_Rep3","B_Rep4"),project = "CombinedConditions")
New<- JoinLayers(combined_seurat_obj, overwrite = TRUE)
New[["RNA"]] <- split(New[["RNA"]], f = New$label) # split only works after join other wise it is throwing an error # @ layers    :List of 2
New <- NormalizeData(New,assay="RNA")
New <- FindVariableFeatures(New,assay="RNA")
New <- ScaleData(New,assay="RNA")
New <- RunPCA(New,assay="RNA")
New <- IntegrateLayers(object = New, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)

New[["RNA"]] <- JoinLayers(New[["RNA"]])
New <- FindNeighbors(New, reduction = "integrated.cca", dims = 1:30)
New <- FindClusters(New, resolution = 1)
New<- RunUMAP(New, dims = 1:30, reduction = "integrated.cca")

saveRDS(New,file="Final_S2_New.rds")


# Plotting
jpeg("Integrate_plot4.jpg", width = 800, height = 600, quality = 100)
new_plot<-DimPlot(New, reduction = "umap", group.by = c("label","seurat_clusters"))
#new_plot<-DimPlot(New, reduction = "umap", split.by = "label")
print(new_plot)
dev.off()

jpeg("Integrate_plot2.jpg", width = 800, height = 600, quality = 100)
new_plot<-DimPlot(New, reduction = "umap", group.by = c("label","seurat_clusters"))
#new_plot<-DimPlot(New, reduction = "umap", split.by = "label")
print(new_plot)
dev.off()


library(Libra)
colnames(New@meta.data)[colnames(New@meta.data) == "seurat_clusters"] <- "cell_type"
DE = run_de(New)
write.csv(DE,"DE_result.csv")




seurat_obj_A2 <- RenameIdents(seurat_obj_A1, `0` = "Oligodendro", `1` = "Endothelial", `2` = "Pericyte", `3` = "Astrocyte", `4` = "Endothelial", `5` = "SMC", `6` = "Endothelial", `7` = "Oligodendro", `8` = "Astrocyte", `9` = "Endothelial", `10` = "Microglia", `11` = "Astrocyte", `12` = "Oligodendro", `13` = "Pericyte" , `14` = "Fibroblast", `15` = "Pericyte", `16` = "Oligodendro", `17` = "Endothelial", `18` = "SMC", `19` = "Oligodendro", `20` = "Endothelial" , `21` = "Endothelial", `22` = "Oligodendro", `23` = "Oligodendro", `24` = "Ependymal", `25` = "Pericyte", `26` = "Oligodendro", `27` = "Oligodendro", `28` = "Oligodendro", `29` = "Oligodendro")
DimPlot(seurat_obj_A2,label = T) + NoLegend()

new_cluster_names <- c(`0` = "Oligodendro", `1` = "Endothelial", `2` = "Pericyte", `3` = "Astrocyte", `4` = "Endothelial", `5` = "SMC", `6` = "Endothelial", `7` = "Oligodendro", `8` = "Astrocyte", `9` = "Endothelial", `10` = "Microglia", `11` = "Astrocyte", `12` = "Oligodendro", `13` = "Pericyte" , `14` = "Fibroblast", `15` = "Pericyte", `16` = "Oligodendro", `17` = "Endothelial", `18` = "SMC", `19` = "Oligodendro", `20` = "Endothelial" , `21` = "Endothelial", `22` = "Oligodendro", `23` = "Oligodendro", `24` = "Ependymal", `25` = "Pericyte", `26` = "Oligodendro", `27` = "Oligodendro", `28` = "Oligodendro", `29` = "Oligodendro")

seurat_obj_A2 <- RenameIdents(seurat_obj_A1, new_cluster_names)
seurat_obj_A2$renamed_clusters <- Idents(seurat_obj_A2)
head(seurat_obj_A2@meta.data)


library(Libra)
str(seurat_obj_A2)
head(seurat_obj_A2@meta.data)
DE = run_de(seurat_obj_A1)
DE = run_de(seurat_obj_A2, cell_type_col = "renamed_clusters", label_col = "label")
write.csv(DE,"DE_result_Type.csv")





===============

combined_seurat_obj_A <- merge(seurat_obj_A1, y = list(seurat_obj_A2,seurat_obj_A3,seurat_obj_A4), add.cell.ids = c("A_Rep1", "A_Rep2", "A_Rep3", "A_Rep4"), project = "Obj_A")
combined_seurat_obj_B <- merge(seurat_obj_B1, y = list(seurat_obj_B2,seurat_obj_B3,seurat_obj_B4), add.cell.ids = c("B_Rep1", "B_Rep2", "B_Rep3", "B_Rep4"), project = "Obj_B")

combined_seurat_obj_A <- merge(seurat_obj_A1, y = list(seurat_obj_A2), add.cell.ids = c("A_Rep1", "A_Rep2"), project = "Obj_A")
combined_seurat_obj_B <- merge(seurat_obj_B1, y = list(seurat_obj_B2), add.cell.ids = c("B_Rep1", "B_Rep2"), project = "Obj_B")
az.list <- list(combined_seurat_obj_A,combined_seurat_obj_B)
ifnb <- merge(x = combined_seurat_obj_A, y = combined_seurat_obj_B)
#Normalize and FindVariableFeatures in all datasets
for (i in 1:length(az.list)) {
  az.list[[i]] <- NormalizeData(az.list[[i]], verbose = TRUE)
  az.list[[i]] <- FindVariableFeatures(az.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = TRUE)
  #az.list[[i]] <- ScaleData(az.list[[i]], verbose = TRUE)
  #az.list[[i]] <- RunPCA(az.list[[i]],verbose = TRUE)
}
# New<- JoinLayers(az.list, overwrite = TRUE)
az.anchors <- FindIntegrationAnchors(az.list, anchor.features = 2000, dims = 1:30
data.input <- GetAssayData(New, assay = "RNA", slot = "data")
az.integrated <- IntegrateData(anchorset =  az.anchors, dims = 1:30)

az.integrated <- ScaleData(az.integrated, verbose = TRUE)
az.integrated <- RunPCA(az.integrated, npcs = 30, verbose = TRUE)

Error:
! GetAssayData doesn't work for multiple layers in v5 assay.

Hi,
Based on the error, it seems that you are using objects with v5 assays. If so, we recommend running integration directly on your 
full object using our new IntegrateLayers() workflow as this is intended to be used on v5 objects. Alternatively, 
if you would like to use your previous code, you can create v3 objects in Seurat v5 by 
setting options(Seurat.object.assay.version = "v3") before running CreateSeuratObject.


elbowPlot(az.integrated, ndims=50)
az.integrated <- RunUMAP(az.integrated, reduction = "pca", dims = 1:30)
az.integrated <- RunTSNE(az.integrated, reduction = "pca", check_duplicates = FALSE, dims = 1:30)
az.integrated <- FindNeighbors(az.integrated, reduction = "pca", dims = 1:30)
az.integrated <- FindClusters(az.integrated, resolution = 0.35)

#save RDS for downstream analysis
saveRDS(az.integrated, "skin_integrated.rds")
	

================================================================
# Add the RNA as assay 
DefaultAssay(combined_seurat_obj) <- "RNA"
New<- JoinLayers(combined_seurat_obj, overwrite = TRUE)
New <- NormalizeData(New,assay="RNA")
New <- FindVariableFeatures(New,assay="RNA")
# After this run the integration step
New <- ScaleData(New,assay="RNA")
New <- RunPCA(New,assay="RNA")

New <- FindNeighbors(New,assay="RNA")
New <- FindClusters(New,assay="RNA")
New <- RunUMAP(New, dims = 1:10)
====================================================================

New[["RNA"]] <- split(New[["RNA"]], f = New$label)
New <- IntegrateLayers(object = New, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
New[["RNA"]] <- JoinLayers(New[["RNA"]])
New <- FindNeighbors(New, reduction = "integrated.cca", dims = 1:30)
New <- FindClusters(New, resolution = 1)
New<- RunUMAP(New, dims = 1:30, reduction = "integrated.cca")


# Plotting
jpeg("Integrate_plot.jpg", width = 800, height = 600, quality = 100)
new_plot<-DimPlot(New, reduction = "umap", group.by = c("label","seurat_clusters"))
new_plot<-DimPlot(New, reduction = "umap", split.by = "label")
print(new_plot)
dev.off()

jpeg("Integrate_plot_Sep.jpg", width = 800, height = 600, quality = 100)
new_plot<-DimPlot(New, reduction = "umap", split.by = "label")
print(new_plot)
dev.off()



# Start Libra
library(Libra)
colnames(New@meta.data)[colnames(New@meta.data) == "seurat_clusters"] <- "cell_type"
DE = run_de(New)
write.csv(DE,"DE_result.csv")

jpeg("TE_plot.jpg", width = 800, height = 600, quality = 100)
new_plot<-FeaturePlot(New, features = c("SoloTE-hAT-Blackjack"))
print(new_plot)
dev.off()


jpeg("Ctrl_plot_3.jpg", width = 800, height = 600, quality = 100)
new_plot<-FeaturePlot(seurat_obj_B3, features = c("SoloTE-hAT-Blackjack","SoloTE-hAT-Tag1","SoloTE-MULE-MuDR","SoloTE-hAT-Tag1","SoloTE-MULE-MuDR","SoloTE-PiggyBac","SoloTE-Dong-R4","SoloTE-Helitron","SoloTE-TcMar","SoloTE-tRNA-Deu","SoloTE-hAT-Ac"))
print(new_plot)
dev.off()














