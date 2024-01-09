#./configure --prefix=$HOME/local/R --enable-R-shlib --without-x --with-pcre1
export PATH="$HOME/local/R/bin:$PATH"
export RSTUDIO_WHICH_R="$HOME/local/R/lib/R/bin/R"

library(clustree)
library(dplyr)
library(Seurat)
library(DESeq2)
library(cluster)
suppressPackageStartupMessages({
    library(rlang)
})
library(ggplot2)
library(dplot)
library(grid)
library(SignacX)
library(future)

#load data

mj001_data <- Read10X(data.dir = "~/Documents/workstation/data/MJ001") #1
mj002_data <- Read10X(data.dir = "~/Documents/workstation/data//MJ002") #2
mj003_data <- Read10X(data.dir = "~/Documents/workstation/data//MJ003") #3
mj005_data <- Read10X(data.dir = "~/Documents/workstation/data//MJ005") #4
mj006_data <- Read10X(data.dir = "~/Documents/workstation/data//MJ006") #5
mj007_data <- Read10X_h5("~/Documents/workstation/data//MJ007/filtered_feature_bc_matrix.h5") #12
mj008_data <- Read10X_h5("~/Documents/workstation/data//MJ008/filtered_feature_bc_matrix.h5") #6
mj009_data <- Read10X(data.dir = "~/Documents/workstation/data//MJ009") #7

mj016_data <- Read10X_h5("~/Documents/workstation/data//MJ016/filtered_feature_bc_matrix.h5") #8
mj017_data <- Read10X_h5("~/Documents/workstation/data//MJ017/filtered_feature_bc_matrix.h5") #9
mj018_data <- Read10X_h5("~/Documents/workstation/data//MJ018/filtered_feature_bc_matrix.h5") #10
mj019_data <- Read10X_h5("~/Documents/workstation/data//MJ019/filtered_feature_bc_matrix.h5") #11

#downsample size
downsize  = 6400


#processing individual data


drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj001_data,), value = FALSE)
mj001_count.data <- mj001_data[-drop.index, ]
mj001_me <- CreateSeuratObject(counts = mj001_count.data, project = "Pt15_POD1194", min.cells = 5)
mj001_me$pos <- "Pt15_POD1194"
mj001_me$mt <- PercentageFeatureSet(mj001_me, pattern = "^MT-")
mj001_me_Q1 = quantile(mj001_me$nFeature_RNA)[2]
mj001_me_Q3 = quantile(mj001_me$nFeature_RNA)[4]
mj001_me_interQ = (mj001_me_Q3-mj001_me_Q1)*1.5
mj001_me_upperb = mj001_me_Q3 + mj001_me_interQ
mj001_me_lowerb = mj001_me_Q1 - mj001_me_interQ
#write.csv(mj001_me$orig.ident,'~/Desktop/int/cell/mj001_cells.csv')

##########################################Mj002##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj002_data,), value = FALSE)
mj002_count.data <- mj002_data[-drop.index, ]
mj002_me <- CreateSeuratObject(counts = mj002_count.data, project = "Pt13_POD1032_IEL", min.cells = 5)
mj002_me$pos <- "Pt13_POD1032_IEL"
mj002_me$mt <- PercentageFeatureSet(mj002_me, pattern = "^MT-")
mj002_me_Q1 = quantile(mj002_me$nFeature_RNA)[2]
mj002_me_Q3 = quantile(mj002_me$nFeature_RNA)[4]
mj002_me_interQ = (mj002_me_Q3-mj002_me_Q1)*1.5
mj002_me_upperb = mj002_me_Q3 + mj002_me_interQ
mj002_me_lowerb = mj002_me_Q1 - mj002_me_interQ
#write.csv(mj002_me$orig.ident,'~/Desktop/int/cell/mj002_cells.csv')

##########################################Mj003##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj003_data,), value = FALSE)
mj003_count.data <- mj003_data[-drop.index, ]
mj003_me <- CreateSeuratObject(counts = mj003_count.data, project = "Pt13_POD1032_LPL", min.cells = 5)
mj003_me$pos <- "Pt13_POD1032_LPL"
mj003_me$mt <- PercentageFeatureSet(mj003_me, pattern = "^MT-")
mj003_me_Q1 = quantile(mj003_me$nFeature_RNA)[2]
mj003_me_Q3 = quantile(mj003_me$nFeature_RNA)[4]
mj003_me_interQ = (mj003_me_Q3-mj003_me_Q1)*1.5
mj003_me_upperb = mj003_me_Q3 + mj003_me_interQ
mj003_me_lowerb = mj003_me_Q1 - mj003_me_interQ
#write.csv(mj003_me$orig.ident,'~/Desktop/int/cell/mj003_cells.csv')

##########################################Mj005##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj005_data,), value = FALSE)
mj005_count.data <- mj005_data[-drop.index, ]
mj005_me <- CreateSeuratObject(counts = mj005_count.data, project = "Pt14_POD1764", min.cells = 5)
mj005_me$pos <- "Pt14_POD1764"
mj005_me$mt <- PercentageFeatureSet(mj005_me, pattern = "^MT-")
mj005_me_Q1 = quantile(mj005_me$nFeature_RNA)[2]
mj005_me_Q3 = quantile(mj005_me$nFeature_RNA)[4]
mj005_me_interQ = (mj005_me_Q3-mj005_me_Q1)*1.5
mj005_me_upperb = mj005_me_Q3 + mj005_me_interQ
mj005_me_lowerb = mj005_me_Q1 - mj005_me_interQ
#write.csv(mj005_me$orig.ident,'~/Desktop/int/cell/mj005_cells.csv')

##########################################Mj006##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj006_data,), value = FALSE)
mj006_count.data <- mj006_data[-drop.index, ]
mj006_me <- CreateSeuratObject(counts = mj006_count.data, project = "Pt21_POD626", min.cells = 5)
mj006_me$pos <- "Pt21_POD626"
mj006_me$mt <- PercentageFeatureSet(mj006_me, pattern = "^MT-")
mj006_me_Q1 = quantile(mj006_me$nFeature_RNA)[2]
mj006_me_Q3 = quantile(mj006_me$nFeature_RNA)[4]
mj006_me_interQ = (mj006_me_Q3-mj006_me_Q1)*1.5
mj006_me_upperb = mj006_me_Q3 + mj006_me_interQ
mj006_me_lowerb = mj006_me_Q1 - mj006_me_interQ
#write.csv(mj006_me$orig.ident,'~/Desktop/int/cell/mj006_cells.csv')

##########################################Mj007##########################################

mj007_me <- CreateSeuratObject(counts = mj007_data, project = "D251", min.cells = 5)
mj007_me$clonotype <- 'Un'
mj007_me$pos <- "D251"
mj007_me$pre <- "Un"
mj007_me$post <- "Un"
mj007_me$pre_post <- "Un; Un"
mj007_me$mt <- PercentageFeatureSet(mj007_me, pattern = "^MT-")
mj007_me_Q1 = quantile(mj007_me$nFeature_RNA)[2]
mj007_me_Q3 = quantile(mj007_me$nFeature_RNA)[4]
mj007_me_interQ = (mj007_me_Q3-mj007_me_Q1)*1.5
mj007_me_upperb = mj007_me_Q3 + mj007_me_interQ
mj007_me_lowerb = mj007_me_Q1 - mj007_me_interQ
mj007_me <- subset(mj007_me, subset = nFeature_RNA > mj007_me_lowerb & nFeature_RNA < mj007_me_upperb  & mt < 15, downsample = downsize)
#mj007_me <- NormalizeData(mj007_me, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj007_me <- FindVariableFeatures(mj007_me, selection.method = "vst", nfeatures = 5000)
#write.csv(mj007_me$orig.ident,'~/Desktop/int/cell/mj007_cells.csv'))
#mj007_me <- SCTransform(mj007_me, verbose = FALSE)

##########################################Mj008##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj008_data,), value = FALSE)
mj008_count.data <- mj008_data[-drop.index, ]
mj008_me <- CreateSeuratObject(counts = mj008_count.data, project = "Pt04_POD1606_IEL", min.cells = 5)
mj008_me$pos <- "Pt04_POD1606_IEL"
mj008_me$mt <- PercentageFeatureSet(mj008_me, pattern = "^MT-")
mj008_me_Q1 = quantile(mj008_me$nFeature_RNA)[2]
mj008_me_Q3 = quantile(mj008_me$nFeature_RNA)[4]
mj008_me_interQ = (mj008_me_Q3-mj008_me_Q1)*1.5
mj008_me_upperb = mj008_me_Q3 + mj008_me_interQ
mj008_me_lowerb = mj008_me_Q1 - mj008_me_interQ
#write.csv(mj008_me$orig.ident,'~/Desktop/int/cell/mj008_cells.csv')

##########################################Mj009##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj009_data,), value = FALSE)
mj009_count.data <- mj009_data[-drop.index, ]
mj009_me <- CreateSeuratObject(counts = mj009_count.data, project = "Pt04_POD1606_LPL", min.cells = 5)
mj009_me$pos <- "Pt04_POD1606_LPL"
mj009_me$mt <- PercentageFeatureSet(mj009_me, pattern = "^MT-")
mj009_me_Q1 = quantile(mj009_me$nFeature_RNA)[2]
mj009_me_Q3 = quantile(mj009_me$nFeature_RNA)[4]
mj009_me_interQ = (mj009_me_Q3-mj009_me_Q1)*1.5
mj009_me_upperb = mj009_me_Q3 + mj009_me_interQ
mj009_me_lowerb = mj009_me_Q1 - mj009_me_interQ
#write.csv(mj009_me$orig.ident,'~/Desktop/int/cell/mj009_cells.csv')

##########################################Mj016##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj016_data,), value = FALSE)
mj016_count.data <- mj016_data[-drop.index, ]
mj016_me <- CreateSeuratObject(counts = mj016_count.data, project = "Pt16_POD1004_IEL", min.cells = 5)
mj016_me$pos <- "Pt16_POD1004_IEL"
mj016_me$mt <- PercentageFeatureSet(mj016_me, pattern = "^MT-")
mj016_me_Q1 = quantile(mj016_me$nFeature_RNA)[2]
mj016_me_Q3 = quantile(mj016_me$nFeature_RNA)[4]
mj016_me_interQ = (mj016_me_Q3-mj016_me_Q1)*1.5
mj016_me_upperb = mj016_me_Q3 + mj016_me_interQ
mj016_me_lowerb = mj016_me_Q1 - mj016_me_interQ
#write.csv(mj016_me$orig.ident,'~/Desktop/int/cell/mj016_cells.csv')


##########################################Mj017##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj017_data,), value = FALSE)
mj017_count.data <- mj017_data[-drop.index, ]
mj017_me <- CreateSeuratObject(counts = mj017_count.data, project = "Pt16_POD1004_LPL", min.cells = 5)
mj017_me$pos <- "Pt16_POD1004_LPL"
mj017_me$mt <- PercentageFeatureSet(mj017_me, pattern = "^MT-")
mj017_me_Q1 = quantile(mj017_me$nFeature_RNA)[2]
mj017_me_Q3 = quantile(mj017_me$nFeature_RNA)[4]
mj017_me_interQ = (mj017_me_Q3-mj017_me_Q1)*1.5
mj017_me_upperb = mj017_me_Q3 + mj017_me_interQ
mj017_me_lowerb = mj017_me_Q1 - mj017_me_interQ
#write.csv(mj017_me$orig.ident,'~/Desktop/int/cell/mj017_cells.csv')

##########################################Mj018##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj018_data,), value = FALSE)
mj018_count.data <- mj018_data[-drop.index, ]
mj018_me <- CreateSeuratObject(counts = mj018_count.data, project = "Pt21_POD1145_IEL", min.cells = 5)
mj018_me$pos <- "Pt21_POD1145_IEL"
mj018_me$mt <- PercentageFeatureSet(mj018_me, pattern = "^MT-")
mj018_me_Q1 = quantile(mj018_me$nFeature_RNA)[2]
mj018_me_Q3 = quantile(mj018_me$nFeature_RNA)[4]
mj018_me_interQ = (mj018_me_Q3-mj018_me_Q1)*1.5
mj018_me_upperb = mj018_me_Q3 + mj018_me_interQ
mj018_me_lowerb = mj018_me_Q1 - mj018_me_interQ
#write.csv(mj018_me$orig.ident,'~/Desktop/int/cell/mj018_cells.csv')

##########################################Mj019##########################################

drop.index <- grep(pattern = "^TR[AB]|^IG[HLK][A-Z][1-9]{0,1}", x = rownames(mj019_data,), value = FALSE)
mj019_count.data <- mj019_data[-drop.index, ]
mj019_me <- CreateSeuratObject(counts = mj019_count.data, project = "Pt21_POD1145_LPL", min.cells = 5)
mj019_me$pos <- "Pt21_POD1145_LPL"
mj019_me$mt <- PercentageFeatureSet(mj019_me, pattern = "^MT-")
mj019_me_Q1 = quantile(mj019_me$nFeature_RNA)[2]
mj019_me_Q3 = quantile(mj019_me$nFeature_RNA)[4]
mj019_me_interQ = (mj019_me_Q3-mj019_me_Q1)*1.5
mj019_me_upperb = mj019_me_Q3 + mj019_me_interQ
mj019_me_lowerb = mj019_me_Q1 - mj019_me_interQ
#write.csv(mj019_me$orig.ident,'~/Desktop/int/cell/mj019_cells.csv')


#category meta data function
cat <- function(info) {
    categories <- info[,7]
    categories <- as.character(categories)
    categories[info[,5] == 'CD4_HvG' & info[,6] == "CD4_H'vG"] <- "Persistent HvG"
    categories[info[,5] == 'CD8_HvG' & info[,6] == "CD8_H'vG"] <- "Persistent HvG"
    categories[info[,5] == 'CD4_HvG' & info[,6] == "CD4_nonH'vG"] <- "Tolerant HvG"
    categories[info[,5] == 'CD8_HvG' & info[,6] == "CD8_nonH'vG"] <- "Tolerant HvG"
    categories[info[,5] == 'CD4_HvG' & info[,6] == "Un"] <- "Missing HvG"
    categories[info[,5] == 'CD8_HvG' & info[,6] == "Un"] <- "Missing HvG"
    categories[info[,5] == 'CD4_nonHvG' & info[,6] == "CD4_H'vG"] <- "Acquired H'vG"
    categories[info[,5] == 'CD8_nonHvG' & info[,6] == "CD8_H'vG"] <- "Acquired H'vG"
    categories[info[,5] == 'Un' & info[,6] == "CD4_H'vG"] <- "De novo H'vG"
    categories[info[,5] == 'Un' & info[,6] == "CD8_H'vG"] <- "De novo H'vG"
    categories[info[,5] == 'CD4_nonHvG' & info[,6] == "CD4_nonH'vG"] <- "Persistent nonHvG"
    categories[info[,5] == 'CD8_nonHvG' & info[,6] == "CD8_nonH'vG"] <- "Persistent nonHvG"
    categories[categories == 'Un; Un'] <- "Others"
    categories <- as.factor(categories)
    return(categories)
}



##########################################DownSample##########################################

info = read.csv('~/Desktop/int/cell/mj005.csv')
mj005_me$clonotype <- info[,4]
mj005_me$pre <- info[,5]
mj005_me$post <- info[,6]
mj005_me$pre_post <- info[,7]
mj005_me$categories <- cat(info)

mj005_me <- subset(mj005_me, subset = nFeature_RNA > mj005_me_lowerb & nFeature_RNA < mj005_me_upperb  & mt < 15)
Idents(mj005_me) <- 'categories'
mj005_me_hvg <- subset(mj005_me,idents = c("Persistent HvG","Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG"))
mj005_me_un <- subset(mj005_me,idents = 'Others')
mj005_me_un_sub <- subset(mj005_me_un, downsample = downsize-length(mj005_me_hvg$orig.ident))
mj005_me <- merge(mj005_me_hvg, mj005_me_un_sub,  project = "Pt14_POD1764")
Idents(mj005_me) <- 'orig.ident'



#Idents(mj005_me) <- 'categories'
#mj005_me <- subset(mj005_me, subset = nFeature_RNA > mj005_me_lowerb & nFeature_RNA < mj005_me_upperb  & mt < 15)
#mj005_me_hvg <- subset(mj005_me,idents = c("Persistent HvG","Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG"))
#mj005_me_un <- subset(mj005_me,idents = 'Others')
#mj005_me_un_sub <- subset(mj005_me_un, downsample = downsize-length(mj005_me_hvg$orig.ident))
#mj005_me <- merge(mj005_me_hvg, mj005_me_un_sub,  project = "Pt14_POD1764")
#Idents(mj005_me) <- 'orig.ident'
#mj005_me <- NormalizeData(mj005_me, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj005_me <- FindVariableFeatures(mj005_me, selection.method = "vst", nfeatures = 5000)
#mj005_me <- SCTransform(mj005_me, verbose = FALSE)



info = read.csv('~/Desktop/int/cell/mj006.csv')
mj006_me$clonotype <- info[,4]
mj006_me$pre <- info[,5]
mj006_me$post <- info[,6]
mj006_me$pre_post <- info[,7]
mj006_me$categories <- cat(info)
mj006_me <- subset(mj006_me, subset = nFeature_RNA > mj006_me_lowerb & nFeature_RNA < mj006_me_upperb  & mt < 15)

Idents(mj006_me) <- 'categories'
mj006_me_hvg <- subset(mj006_me,idents = c("Persistent HvG","Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG"))
mj006_me_un <- subset(mj006_me,idents = 'Others')
mj006_me_un_sub <- subset(mj006_me_un, downsample = downsize-length(mj006_me_hvg$orig.ident))
mj006_me <- merge(mj006_me_hvg, mj006_me_un_sub,  project = "Pt21_POD626")
Idents(mj006_me) <- 'orig.ident'



#mj006_me <- subset(mj006_me, subset = nFeature_RNA > mj006_me_lowerb & nFeature_RNA < mj006_me_upperb  & mt < 15)
#mj006_me_hvg <- subset(mj006_me,idents = c("Persistent HvG","Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG"))
#mj006_me_un <- subset(mj006_me,idents = 'Others')
#mj006_me_un_sub <- subset(mj006_me_un, downsample = downsize-length(mj006_me_hvg$orig.ident))
#mj006_me <- merge(mj006_me_hvg, mj006_me_un_sub,  project = "Pt21_POD626")
#Idents(mj006_me) <- integrated8.rpca
#mj006_me <- NormalizeData(mj006_me, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj006_me <- FindVariableFeatures(mj006_me, selection.method = "vst", nfeatures = 5000)
#mj006_me <- SCTransform(mj006_me, verbose = FALSE)



info = read.csv('~/Desktop/int/cell/mj001.csv')
mj001_me$clonotype <- info[,4]
mj001_me$pre <- info[,5]
mj001_me$post <- info[,6]
mj001_me$pre_post <- info[,7]
mj001_me$categories <- cat(info)
mj001_me <- subset(mj001_me, subset = nFeature_RNA > mj001_me_lowerb & nFeature_RNA < mj001_me_upperb  & mt < 15,downsample = downsize)

#mj001_me <- NormalizeData(mj001_me,normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj001_me <- FindVariableFeatures(mj001_me, selection.method = "vst", nfeatures = 5000)
#mj001_me <- SCTransform(mj001_me, verbose = TRUE)

info = read.csv('~/Desktop/int/cell/mj002.csv')
mj002_me$clonotype <- info[,4]
mj002_me$pre <- info[,5]
mj002_me$post <- info[,6]
mj002_me$pre_post <- info[,7]
mj002_me$categories <- cat(info)
mj002_me <- subset(mj002_me, subset = nFeature_RNA > mj002_me_lowerb & nFeature_RNA < mj002_me_upperb  & mt < 15, downsample = downsize)

#mj002_me <- NormalizeData(mj002_me,normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj002_me <- FindVariableFeatures(mj002_me, selection.method = "vst", nfeatures = 5000)
#mj002_me <- SCTransform(mj002_me, verbose = FALSE)
info = read.csv('~/Desktop/int/cell/mj003.csv')
mj003_me$clonotype <- info[,4]
mj003_me$pre <- info[,5]
mj003_me$post <- info[,6]
mj003_me$pre_post <- info[,7]
mj003_me$categories <- cat(info)
mj003_me <- subset(mj003_me, subset = nFeature_RNA > mj003_me_lowerb & nFeature_RNA < mj003_me_upperb  & mt < 15, downsample = downsize)


#mj003_me <- NormalizeData(mj003_me, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj003_me <- FindVariableFeatures(mj003_me, selection.method = "vst", nfeatures = 5000)
#mj003_me <- SCTransform(mj003_me, verbose = FALSE)

info = read.csv('~/Desktop/int/cell/mj008.csv')
mj008_me$clonotype <- info[,4]
mj008_me$pre <- info[,5]
mj008_me$post <- info[,6]
mj008_me$pre_post <- info[,7]
mj008_me$categories <- cat(info)
mj008_me <- subset(mj008_me, subset = nFeature_RNA > mj008_me_lowerb & nFeature_RNA < mj008_me_upperb  & mt < 15, downsample = downsize)


#mj008_me <- NormalizeData(mj008_me,normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj008_me <- FindVariableFeatures(mj008_me, selection.method = "vst", nfeatures = 5000)

#mj008_me <- SCTransform(mj008_me, verbose = FALSE)
info = read.csv('~/Desktop/int/cell/mj009.csv')
mj009_me$clonotype <- info[,4]
mj009_me$pre <- info[,5]
mj009_me$post <- info[,6]
mj009_me$pre_post <- info[,7]
mj009_me$categories <- cat(info)
mj009_me <- subset(mj009_me, subset = nFeature_RNA > mj009_me_lowerb & nFeature_RNA < mj009_me_upperb  & mt < 15, downsample = downsize)




#Idents(mj009_me) <- 'categories'
#mj009_me <- subset(mj009_me, subset = nFeature_RNA > mj009_me_lowerb & nFeature_RNA < mj009_me_upperb  & mt < 15)
#mj009_me_hvg <- subset(mj009_me,idents = c("Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG"))
#mj009_me_un <- subset(mj009_me,idents = 'Others')
#mj009_me_un_sub <- subset(mj009_me_un, downsample = downsize-length(mj009_me_hvg$orig.ident))
#mj009_me <- merge(mj009_me_hvg, mj009_me_un_sub,  project = "Pt04_POD1606_LPL")
#Idents(mj009_me) <- 'orig.ident'
#mj009_me <- NormalizeData(mj009_me, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
#mj009_me <- FindVariableFeatures(mj009_me, selection.method = "vst", nfeatures = 5000)
#mj009_me <- SCTransform(mj009_me, verbose = FALSE)



info = read.csv('~/Desktop/int/cell/mj016.csv')
mj016_me$clonotype <- info[,4]
mj016_me$pre <- info[,5]
mj016_me$post <- info[,6]
mj016_me$pre_post <- info[,7]
mj016_me$categories <- cat(info)
mj016_me <- subset(mj016_me, subset = nFeature_RNA > mj016_me_lowerb & nFeature_RNA < mj016_me_upperb  & mt < 15, downsample = downsize)


info = read.csv('~/Desktop/int/cell/mj017.csv')
mj017_me$clonotype <- info[,4]
mj017_me$pre <- info[,5]
mj017_me$post <- info[,6]
mj017_me$pre_post <- info[,7]
mj017_me$categories <- cat(info)
mj017_me <- subset(mj017_me, subset = nFeature_RNA > mj017_me_lowerb & nFeature_RNA < mj017_me_upperb  & mt < 15, downsample = downsize)


info = read.csv('~/Desktop/int/cell/mj018.csv')
mj018_me$clonotype <- info[,4]
mj018_me$pre <- info[,5]
mj018_me$post <- info[,6]
mj018_me$pre_post <- info[,7]
mj018_me$categories <- cat(info)
mj018_me <- subset(mj018_me, subset = nFeature_RNA > mj018_me_lowerb & nFeature_RNA < mj018_me_upperb  & mt < 15, downsample = downsize)



info = read.csv('~/Desktop/int/cell/mj019.csv')
mj019_me$clonotype <- info[,4]
mj019_me$pre <- info[,5]
mj019_me$post <- info[,6]
mj019_me$pre_post <- info[,7]
mj019_me$categories <- cat(info)
mj019_me <- subset(mj019_me, subset = nFeature_RNA > mj019_me_lowerb & nFeature_RNA < mj019_me_upperb  & mt < 15, downsample = downsize)


##########################################Integrating##########################################

samples <- list( mj001_me,mj002_me,mj003_me,mj005_me,mj006_me,mj008_me,mj009_me,mj016_me,mj017_me,mj018_me,mj019_me,mj007_me)

mj.list <- lapply(X = samples, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 6400)
})

features <- SelectIntegrationFeatures(object.list = samples, nfeatures = 6400)

mj.list <- lapply(X = mj.list, FUN = function(x) {

    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors.rpca <- FindIntegrationAnchors(object.list = mj.list, reduction = "rpca", anchor.features =features)


integrated8.rpca <- IntegrateData(anchorset = anchors.rpca,dims = 1:64)
rm(anchors.rpca)
saveRDS(integrated8.rpca,'~/data/integrated8.rpca.rds')



integrated8.rpca <- ScaleData(integrated8.rpca, verbose = TRUE)
integrated8.rpca <- ScaleData(integrated8.rpca, assay='RNA', verbose = TRUE)
integrated8.rpca <- RunPCA(integrated8.rpca, verbose = TRUE)
integrated8.rpca <- RunUMAP(integrated8.rpca, dims = 1:32)
integrated8.rpca <- FindNeighbors(integrated8.rpca, dims = 1:32)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.64)


integrated8.rpca <- SCTransform(integrated8.rpca, conserve.memory = TRUE)
integrated8.rpca <- RunPCA(integrated8.rpca, verbose = FALSE)
integrated8.rpca <- RunUMAP(integrated8.rpca, dims = 1:30, verbose = FALSE)
integrated8.rpca <- FindNeighbors(integrated8.rpca, dims = 1:30, verbose = FALSE)


################################################################################
#Downstream analysis
#Including specific cell shown etc.
#If you have any questions, please contact the authors

colors <- c('#E91E63','#3F51B5','#00BCD4','#8BC34A','#FFC107','#795548',
            '#9C27B0','#2196F3','#009688','#CDDC39','#FF9800','#9E9E9E',
            '#673AB7','#03A9F4','#4CAF50','#FFEB3B','#FF5722','#607D8B',
            '#F44336','#000000',
            '#F48FB1','#9FA8DA','#80DEEA','#C5E1A5','#FFE082','#BCAAA4',
            '#CE93D8','#90CAF9','#80CBC4','#E6EE9C','#FFCC80','#EEEEEE'
            )

orders <- c('CD4_HvG','CD8_HvG','CD4_nonHvG','CD8_nonHvG','Un')
orders_post <- c("CD4_H'vG","CD8_H'vG","CD4_nonH'vG","CD8_nonH'vG",'Un')
color5 <- c('#DDDDDD','#FF9800','#4CAF50','#3F51B5','#E91E63')
genelist = c("CD79A", "CD74", "CD3E", "CD3D", "CD3G", "CD4", "CD8A", "TRDC", "CD69", "ITGAE", "ITGA1", "CXCR6", "RUNX3", "PRDM1", "ZNF683", "CCR7", "KLF2", "S1PR1", "SELL", "MKI67", "CD28","TMIGD2", "TIGIT","TCF7", "RORA", "RORC", "AHR", "IL23R", "IL23A", "IL22RA1", "IL22", "IL17RA","IL17A","CCL20", "TNF", "IFNG", "IFNGR1","TBX21", "EOMES", "IL4", "IL5", "IL13", "GATA3", "BCL6", "CXCR5", "PDCD1", "ICOS", "IL2RA", "CXCL13", "IL21", "IL6", "STAT3", "IRF4", "IL21R", "IL6R", "BATF", "STAT5A", "MAF", "GZMB", "GZMA","PRF1","GNLY", "KLRD1", "KLRC1", "KLRK1","FOXP3", "CTLA4","IL10","TGFB1", "IL2","NR4A1", "SATB1", "EGR1", "IL4I1", "NFATC1","JUN", "FOS","NFKBIA", "NFKBID","CD83", "PTPN22", "CXCR4", "BTG2", "HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1", "HLA-E", "HLA-F", "LYN", "CCL4")

delist = c("TRDC", "TRDV1", "TRDV2", "TRDV3", "CD69", "ITGAE", "ITGA1", "RUNX3", "CCR7", "KLF2", "S1PR1", "SELL", "MKI67", "CD28", "GZMB", "GZMA", "PRF1", "GNLY", "KLRD1")

color7 <- c('#9E9E9E','#EEEEEE','#FFC107','#CDDC39','#009688','#03A9F4','#673AB7','#E91E63')
order_cat <- c("Persistent HvG","Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG","Others","Control")
color7 <- c('#BBBBBB','#EEEEEE',"#9C27B0","#3F51B5","#00BCD4","#4CAF50","#FF9800","#E91E63")

sel_labels <- c("CD4_HvG; CD4_H'vG","Un; CD4_H'vG","CD4_HvG; Un","CD4_HvG; CD4_nonH'vG","CD4_nonHvG; CD4_nonH'vG","CD8_HvG; CD8_H'vG","Un; CD8_H'vG","CD8_HvG; Un","CD8_HvG; CD8_nonH'vG","CD8_nonHvG; CD8_nonH'vG","Un; Un")
sel_labels_order <- c("CD4_HvG; CD4_H'vG","Un; CD4_H'vG","CD4_HvG; Un","CD4_HvG; CD4_nonH'vG","CD4_nonHvG; CD4_nonH'vG","CD8_HvG; CD8_H'vG","Un; CD8_H'vG","CD8_HvG; Un","CD8_HvG; CD8_nonH'vG","CD8_nonHvG; CD8_nonH'vG","Others","Un; Un")
labels_colors <- c('#EEEEEE','#BBBBBB','#FF9800','#FFEB3B','#00BCD4','#2196F3','#3F51B5','#4CAF50','#009688','#673AB7','#9C27B0','#E91E63')
labels_colors <- c('#EEEEEE',
    '#BBBBBB',

    '#FF9800',
    '#FFEB3B',
    '#00BCD4',
    #'#74b9ff',
    '#2196F3',
    '#3F51B5',

    '#4CAF50',
    #'#81ecec',
    '#009688',
    '#673AB7',
    #'#fd79a8',
    '#9C27B0',
    '#E91E63')

lc = length(table(integrated8.rpca$pre_post))
sel_lab = c()
sel_col = c()
for (i in seq(12)){
    if (sel_labels_order[i]%in% rownames(table(integrated8.rpca$pre_post)))
        {
        print (sel_labels_order[i])
        print (rev(labels_colors)[i])
        sel_lab = append(sel_lab, sel_labels_order[i])
        sel_col = append(sel_col, rev(labels_colors)[i])
        }
}

sel_col = rev(sel_col)


sel_labels_order2 <- c("Persistent HvG","Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG","Others","Unmappable")

labels_colors2 <- c('#EEEEEE','#BBBBBB','#9C27B0','#3F51B5','#00BCD4','#4CAF50','#FF9800','#E91E63')
lc2 = length(table(integrated8.rpca$categories))
sel_lab2 = c()
sel_col2 = c()
for (i in seq(8)){
    if (sel_labels_order2[i]%in% rownames(table(integrated8.rpca$categories)))
        {
        print (sel_labels_order2[i])
        print (rev(labels_colors2)[i])
        sel_lab2 = append(sel_lab2, sel_labels_order2[i])
        sel_col2 = append(sel_col2, rev(labels_colors2)[i])
        }
}
sel_col2 = rev(sel_col2)


ncluster <- integrated8.rpca$integrated_snn_res.0.5
recluster <- function(ncluster){
    ncluster <- as.character(ncluster)
    for(i in seq(length(ncluster))){
        ncluster[i] = paste('c',sprintf("%02d", as.numeric(ncluster[i])+1),sep='')
    }
    ncluster <- as.factor(ncluster)
    return(ncluster)
}

integrated8.rpca$groups <- recluster(ncluster)


integrated8.rpca$repos <- factor(x = integrated8.rpca$pos, levels =c(
'D251',            'Pt14_POD1764',     'Pt15_POD1194',     'Pt21_POD626',
'Pt04_POD1606_IEL','Pt13_POD1032_IEL', 'Pt16_POD1004_IEL', 'Pt21_POD1145_IEL',
'Pt04_POD1606_LPL','Pt13_POD1032_LPL', 'Pt16_POD1004_LPL', 'Pt21_POD1145_LPL'
))



####################Figures Plotting & tables generating########################


pdf("~/Desktop/int/UMAP.pdf", width = 8, height = 4.5)
DimPlot(integrated8.rpca ,label.size = 6,label = TRUE, cols = colors, group.by = 'groups') +NoLegend()
dev.off()

pdf("~/Desktop/int/UMAP_HvG.pdf", width = 8, height = 4.5)
DimPlot(integrated8.rpca , pt.size = 1, order = orders, cols = color5, group.by='pre')
dev.off()

pdf("~/Desktop/int/UMAP_split.pdf", width = 16, height = 9)
DimPlot(integrated8.rpca ,label.size = 3,label = FALSE,cols = colors,split.by = 'repos',ncol = 4, group.by = 'groups')
dev.off()


pdf("~/Desktop/int/UMAP_HvG_split.pdf", width = 16, height = 9)
DimPlot(integrated8.rpca , pt.size = 1, order = orders,cols = color5, group.by='pre',split.by = 'repos',ncol = 4)
dev.off()


pdf("~/Desktop/int/UMAP_H'vG.pdf", width = 8, height = 4.5)
DimPlot(integrated8.rpca , pt.size = 1, order = orders_post, cols = color5, group.by='post')
dev.off()


pdf("~/Desktop/int/UMAP_H'vG_split.pdf", width = 16, height = 9)
DimPlot(integrated8.rpca , pt.size = 1, order = orders_post,cols = color5, group.by='post',split.by = 'repos',ncol = 4)
dev.off()


pdf("~/Desktop/int/UMAP_Pre_Post.pdf", width = 8, height = 4.5)
DimPlot(integrated8.rpca , pt.size = 2, order = sel_lab, cols =sel_col, group.by='pre_post')
dev.off()


pdf("~/Desktop/int/UMAP_Pre_Post_split.pdf", width = 16, height = 9)
DimPlot(integrated8.rpca , pt.size = 1, order = sel_lab, cols =sel_col, group.by='pre_post',split.by = 'repos',ncol = 4)
dev.off()


pdf("~/Desktop/int/UMAP_L18.pdf", width = 8, height = 4.5)
DimPlot(integrated8.rpca ,label.size = 6,label = TRUE, cols = colors, group.by = 'L18') +NoLegend()
dev.off()

pdf("~/Desktop/int/UMAP_L18_split.pdf", width = 16, height = 9)
DimPlot(integrated8.rpca ,label.size = 3,label = FALSE,cols = colors,split.by = 'repos',ncol = 4, group.by = 'L18')
dev.off()


cate <- integrated8.rpca$categories
cate[is.na(cate)] <- 'Control'
integrated8.rpca$cate <- cate
pdf("~/Desktop/int/UMAP_CAT.pdf", width = 8, height = 4.5)
DimPlot(integrated8.rpca , pt.size = 1, order = order_cat , cols = color7, group.by='cate')
dev.off()
#sel_labels_order2#sel_col2
pdf("~/Desktop/int/UMAP_CAT_split.pdf", width = 16, height = 9)
DimPlot(integrated8.rpca , pt.size = 1, order = order_cat ,cols =color7, group.by='cate',split.by = 'repos',ncol = 4)
dev.off()


Idents(integrated8.rpca) <- 'groups'
jpeg("~/Desktop/int/HEAT.jpeg",width = 6400, height = 3600,res=300)
DoHeatmap(integrated8.rpca, features = genelist,assay='RNA', group.colors =colors)
dev.off()



Idents(integrated8.rpca) <- 'groups'
integrated8.rpca.c618 <- subset(integrated8.rpca, idents=c('c06','c18'))
jpeg("~/Desktop/int/HEAT_C618.jpeg",width = 3200, height = 1800,res=300)
DoHeatmap(integrated8.rpca.c618, features = delist,assay='RNA', group.colors =c('#795548','#607D8B'))
dev.off()


Idents(integrated8.rpca) <- 'groups'
genename = 'CD69'
pdf(paste("~/Desktop/int/FEAT_",genename,".pdf",sep = ''), width = 8, height = 4.5)
FeaturePlot(integrated8.rpca,features = genename,min.cutoff = 'q9',order=TRUE)
dev.off()


integrated8.rpca.markers <- FindAllMarkers(integrated8.rpca, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
integrated8_top30 <- integrated8.rpca.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

jpeg(paste("~/Desktop/int/HEAT_markers0.5.jpeg",sep=''), width = 9600, height = length(integrated8_top30$gene)*40,res=300)
DoHeatmap(integrated8.rpca , features = integrated8_top30$gene, group.colors =colors)
dev.off()
write.csv(integrated8.rpca.markers,'~/Desktop/int/markers.csv')


integrated8.rpca$cat <- integrated8.rpca$categories
integrated8.rpca$cat[integrated8.rpca$cat == "Persistent HvG"] <- "01.Persistent HvG"
integrated8.rpca$cat[integrated8.rpca$cat == "Tolerant HvG"] <- "02.Tolerant HvG"
integrated8.rpca$cat[integrated8.rpca$cat == "Missing HvG"] <- "03.Missing HvG"
integrated8.rpca$cat[integrated8.rpca$cat == "Acquired H'vG"] <- "04.Acquired H'vG"
integrated8.rpca$cat[integrated8.rpca$cat == "De novo H'vG"] <- "05.De novo H'vG"
integrated8.rpca$cat[integrated8.rpca$cat == "Persistent nonHvG"] <- "06.Persistent nonHvG"

catc <- c("Persistent HvG","Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG")
Idents(integrated8.rpca) <- 'categories'
integrated8.rpca.cat <- subset(integrated8.rpca, idents=catc)

"#E91E63" "#FF9800" "#4CAF50" "#00BCD4" "#3F51B5" "#9C27B0"
c("#00BCD4","#3F51B5","#4CAF50","#E91E63","#9C27B0","#FF9800")


jpeg("~/Desktop/int/HEAT_cat_genelist.jpeg", width = 1600*3, height = 1600*3,res=300)
DoHeatmap(integrated8.rpca.cat, features = genelist, group.colors = c("#E91E63","#FF9800","#4CAF50","#00BCD4","#3F51B5","#9C27B0"),assay ='RNA',group.by='cat')
dev.off()


DefaultAssay(integrated8.rpca.cat) <- 'RNA'
integrated8.rpca.cat.markers <- FindAllMarkers(integrated8.rpca.cat, only.pos = TRUE)
write.csv(integrated8.rpca.cat.markers,'~/Desktop/int/DE_list.csv')
write.csv(integrated8.rpca@assays$RNA@data,'~/Desktop/int/data.csv')

mj018,'~/Desktop/int/MJ018_data.csv')

integrated8_top20 <- integrated8.rpca.cat.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
jpeg(paste("~/Desktop/int/HEAT_cat_makers.jpeg",sep=''), width = 1600*3, height = 2000*3,res=300)
DoHeatmap(integrated8.rpca.cat, features = integrated8_top20$gene, group.colors =c("#00BCD4","#3F51B5","#4CAF50","#E91E63","#9C27B0","#FF9800"),group.by='categories',assay='RNA')
dev.off()


mj001_marker_hvg <- FindAllMarkers(mj001_me)
mj002_marker_hvg <- FindAllMarkers(mj002_me)
mj003_marker_hvg <- FindAllMarkers(mj003_me)
mj005_marker_hvg <- FindAllMarkers(mj005_me)
mj006_marker_hvg <- FindAllMarkers(mj006_me)

mj007_marker_hvg <- FindAllMarkers(mj007_me)
mj008_marker_hvg <- FindAllMarkers(mj008_me)
mj009_marker_hvg <- FindAllMarkers(mj009_me)


top20 <- mj001_marker_hvg %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
mj001_me <- ScaleData(mj001_me, verbose = TRUE)
jpeg("~/Desktop/DE_list/mj001_marker_hvg.jpeg",width = 1000, height = 1000)
DoHeatmap(mj001_me, features = top20$gene) + NoLegend()
dev.off()
top20 <- mj002_marker_hvg %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
mj002_me <- ScaleData(mj002_me, verbose = TRUE)
jpeg("~/Desktop/DE_list/mj002_marker_hvg.jpeg",width = 1000, height = 1000)
DoHeatmap(mj002_me, features = top20$gene) + NoLegend()
dev.off()
top20 <- mj003_marker_hvg %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
mj003_me <- ScaleData(mj003_me, verbose = TRUE)
jpeg("~/Desktop/DE_list/mj003_marker_hvg.jpeg",width = 1000, height = 1000)
DoHeatmap(mj003_me, features = top20$gene) + NoLegend()
dev.off()
top20 <- mj005_marker_hvg %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
mj005_me <- ScaleData(mj005_me, verbose = TRUE)
jpeg("~/Desktop/DE_list/mj005_marker_hvg.jpeg",width = 1000, height = 1000)
DoHeatmap(mj005_me, features = top20$gene) + NoLegend()
dev.off()
top20 <- mj006_marker_hvg %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
mj006_me <- ScaleData(mj006_me, verbose = TRUE)
jpeg("~/Desktop/DE_list/mj006_marker_hvg.jpeg",width = 1000, height = 1000)
DoHeatmap(mj006_me, features = top20$gene) + NoLegend()
dev.off()

write.csv(mj001_marker_hvg,'~/Desktop/DE_list/mj001_marker_hvg.csv')
write.csv(mj002_marker_hvg,'~/Desktop/DE_list/mj002_marker_hvg.csv')
write.csv(mj003_marker_hvg,'~/Desktop/DE_list/mj003_marker_hvg.csv')
write.csv(mj005_marker_hvg,'~/Desktop/DE_list/mj005_marker_hvg.csv')
write.csv(mj006_marker_hvg,'~/Desktop/DE_list/mj006_marker_hvg.csv')

mj001_marker_cd4hvg <- FindAllMarkers(subset(mj001_me, only.pos=TRUE, idents=c('CD4_HvG','CD4_nonHvG')))
mj002_marker_cd4hvg <- FindAllMarkers(subset(mj002_me, only.pos=TRUE, idents=c('CD4_HvG','CD4_nonHvG')))
mj003_marker_cd4hvg <- FindAllMarkers(subset(mj003_me, only.pos=TRUE, idents=c('CD4_HvG','CD4_nonHvG')))
mj005_marker_cd4hvg <- FindAllMarkers(subset(mj005_me, only.pos=TRUE, idents=c('CD4_HvG','CD4_nonHvG')))
mj006_marker_cd4hvg <- FindAllMarkers(subset(mj006_me, only.pos=TRUE, idents=c('CD4_HvG','CD4_nonHvG')))

write.csv(mj001_marker_cd4hvg,'~/Desktop/DE_list/mj001_marker_cd4hvg.csv')
write.csv(mj002_marker_cd4hvg,'~/Desktop/DE_list/mj002_marker_cd4hvg.csv')
write.csv(mj003_marker_cd4hvg,'~/Desktop/DE_list/mj003_marker_cd4hvg.csv')
write.csv(mj005_marker_cd4hvg,'~/Desktop/DE_list/mj005_marker_cd4hvg.csv')
write.csv(mj006_marker_cd4hvg,'~/Desktop/DE_list/mj006_marker_cd4hvg.csv')


top20 <- mj001_marker_cd8hvg %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
jpeg("~/Desktop/DE_list/mj001_marker_cd8hvg.jpeg",width = 1000, height = 1000)
DoHeatmap(subset(mj001_me, only.pos=TRUE, idents=c('CD8_HvG','CD8_nonHvG')), features = top20$gene) + NoLegend()
dev.off()
top20 <- mj002_marker_cd8hvg %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
jpeg("~/Desktop/DE_list/mj002_marker_cd8hvg.jpeg",width = 1000, height = 1000)
DoHeatmap(subset(mj002_me, only.pos=TRUE, idents=c('CD8_HvG','CD8_nonHvG')), features = top20$gene) + NoLegend()
dev.off()
top20 <- mj003_marker_cd8hvg %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
jpeg("~/Desktop/DE_list/mj003_marker_cd8hvg.jpeg",width = 1000, height = 1000)
DoHeatmap(subset(mj003_me, only.pos=TRUE, idents=c('CD8_HvG','CD8_nonHvG')), features = top20$gene) + NoLegend()
dev.off()
top20 <- mj005_marker_cd8hvg %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
jpeg("~/Desktop/DE_list/mj005_marker_cd8hvg.jpeg",width = 1000, height = 1000)
DoHeatmap(subset(mj005_me, only.pos=TRUE, idents=c('CD8_HvG','CD8_nonHvG')), features = top20$gene) + NoLegend()
dev.off()
top20 <- mj006_marker_cd8hvg %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
jpeg("~/Desktop/DE_list/mj006_marker_cd8hvg.jpeg",width = 1000, height = 1000)
DoHeatmap(subset(mj006_me, only.pos=TRUE, idents=c('CD8_HvG','CD8_nonHvG')), features = top20$gene) + NoLegend()
dev.off()


DefaultAssay(integrated8.rpca) <- 'RNA'
test_gene = c('GZMK', 'ITGB2', 'CCL4','KLRG1','CRTAM','ITGAE','ITGA1','CD7','IL7R','KLRB1','CCR6',
    'CD69', 'IL23R', 'PDCD1', 'DUSP6', 'RGS1', 'KCNK5', 'CXCR6', 'IL10', 'IL2', 'CA10',
    'CD3E','CD3D','CD3G', 'TRDC', 'TRBC1', 'CD4', 'CD8A')

jpeg("~/Desktop/int/Int_split.jpeg", width = 3600, height = 2400, res=300)
DotPlot(integrated8.rpca, features = rev(test_gene) , group.by = 'groups', cols = c('#DDDDDD','#3F51B5'),assay='RNA') + RotatedAxis()
dev.off()

gamma = c("TRGV2", "TRGV3", "TRGV4", "TRGV5", "TRGV8", "TRGV9", "TRGV11", "TRGV1", "TRGV5P", "TRGV6", "TRGV7", "TRGV10", "TRGVA", "TRGVB, TRGJ1", "TRGJ2", "TRGJP", "TRGJP1", "TRGJP2")
delta = c("TRDV1", "TRDV2", "TRDV3", "TRAV14DV4", "TRAV29DV5", "TRAV23DV6", "TRAV36DV7", "TRAV38-2DV8 TRDJ1", "TRDJ2", "TRDJ3", "TRDJ4")


jpeg("~/Desktop/int/Int_gamma.jpeg", width = 2400, height = 1200, res=300)
DotPlot(integrated8.rpca, features = rev(gamma) , group.by = 'repos', cols = c('#DDDDDD','#3F51B5'),assay='RNA') + RotatedAxis()
dev.off()

jpeg("~/Desktop/int/Int_delta.jpeg", width = 1800, height = 1200, res=300)
DotPlot(integrated8.rpca, features = rev(delta) , group.by = 'repos', cols = c('#DDDDDD','#3F51B5'),assay='RNA') + RotatedAxis()
dev.off()


mj001_marker_cd8hvg <- FindAllMarkers(subset(mj001_me, only.pos=TRUE, idents=c('CD8_HvG','CD8_nonHvG')))
mj002_marker_cd8hvg <- FindAllMarkers(subset(mj002_me, only.pos=TRUE, idents=c('CD8_HvG','CD8_nonHvG')))
mj003_marker_cd8hvg <- FindAllMarkers(subset(mj003_me, only.pos=TRUE, idents=c('CD8_HvG','CD8_nonHvG')))
mj005_marker_cd8hvg <- FindAllMarkers(subset(mj005_me, only.pos=TRUE, idents=c('CD8_HvG','CD8_nonHvG')))
mj006_marker_cd8hvg <- FindAllMarkers(subset(mj006_me, only.pos=TRUE, idents=c('CD8_HvG','CD8_nonHvG')))

write.csv(mj001_marker_cd8hvg,'~/Desktop/DE_list/mj001_marker_cd8hvg.csv')
write.csv(mj002_marker_cd8hvg,'~/Desktop/DE_list/mj002_marker_cd8hvg.csv')
write.csv(mj003_marker_cd8hvg,'~/Desktop/DE_list/mj003_marker_cd8hvg.csv')
write.csv(mj005_marker_cd8hvg,'~/Desktop/DE_list/mj005_marker_cd8hvg.csv')
write.csv(mj006_marker_cd8hvg,'~/Desktop/DE_list/mj006_marker_cd8hvg.csv')


integrated.hvg.qui.cd4 <- subset(integrated.hvg.qui,  idents=c('CD4_HvG','CD4_nonHvG'))
ma <- FindAllMarkers(integrated.hvg.qui.cd8,assay = 'RNA',test.use = 'DESeq2')
top30 <- ma %>% group_by(cluster) %>% top_n(n = 30, wt = avg_logFC)


cell_file = 'int/cell/mj001_cells.csv'
cl1 = makelabel(cell_file,pt15_post_output_2,labels,pt15_tcr_file, True)
count_clone(cl1)
cl1.to_csv('int/cell/mj001.csv')

cell_file = 'int/cell/mj002_cells.csv'
cl2 = makelabel(cell_file,pt13_iel_post_output_2,labels,pt13_iel_tcr_file,True)
count_clone(cl2)
cl2.to_csv('int/cell/mj002.csv')

cell_file = 'int/cell/mj003_cells.csv'
cl3 = makelabel(cell_file,pt13_lpl_post_output_2,labels,pt13_lpl_tcr_file,True)
count_clone(cl3)
cl3.to_csv('int/cell/mj003.csv')

cell_file = 'int/cell/mj005_cells.csv'
cl5 = makelabel(cell_file,pt14_post_output,labels,pt14_tcr_file,True)
count_clone(cl5)
cl5.to_csv('int/cell/mj005.csv')

cell_file = 'int/cell/mj006_cells.csv'
cl6 = makelabel(cell_file,pt21_post_output_2,labels,pt21_tcr_file,True)
count_clone(cl6)
cl6.to_csv('int/cell/mj006.csv')

cell_file = 'int/cell/mj008_cells.csv'
cl8 = makelabel(cell_file,pt4_iel_post_output,labels,pt4_iel_tcr_file,True)
count_clone(cl8)
cl8.to_csv('int/cell/mj008.csv')

cell_file = 'int/cell/mj009_cells.csv'
cl9 = makelabel(cell_file,pt4_lpl_post_output,labels,pt4_lpl_tcr_file,True)
count_clone(cl9)
cl9.to_csv('int/cell/mj009.csv')



jpeg("~/Desktop/int/MJ010-11_DOT_gamma.jpeg", width = 2400, height = 800, res=300)
DotPlot(integrated8.bm, features = rev(gamma) , group.by = 'name', cols = c('#DDDDDD','#3F51B5'),assay ='RNA') + RotatedAxis()
dev.off()

jpeg("~/Desktop/int/MJ010-11_DOT_delta.jpeg", width = 1800, height = 1200, res=300)
DotPlot(integrated8.bm, features = rev(delta) , group.by = 'name', cols = c('#DDDDDD','#3F51B5'),assay ='RNA') + RotatedAxis()
dev.off()


Idents(integrated8.rpca) <- 'orig.ident'
integrated8.rpca <- FindNeighbors(integrated8.rpca, dims = 1:32)

integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.05,reduction = "umap")
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.1)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.15)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.2)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.25)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.3)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.35)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.4)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.45)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.5)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.55)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.6)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.7)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.8)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.85)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.9)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.95)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 1)

integrated8.rpca <- RunUMAP(integrated8.rpca, dims = 1:30)



suppressPackageStartupMessages({
  library(rlang)
})
library(ggplot2)
DoMultiBarHeatmap(integrated8.rpca, assay = 'RNA', features = integrated8_top30$gene, group.by='groups', additional.group.by = 'pos',cols.use = clos.use)
DoMultiBarHeatmap(integrated8.rpca.cat, assay = 'RNA', features = genelist, group.by='categories', additional.group.by = c('pos','integrated_snn_res.0.75'))

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}

library(grid)

library(cluster, quietly = TRUE)
dist.matrix <- dist(x = Embeddings(object = integrated8.rpca)[, dims])
clusters <- dataset$celltype
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
dataset$sil <- sil[, 3]

# mixing metric
max.k <- 300
mm <- max.k - MixingMetric(object = dataset, grouping.var = "replicate", reduction = reduction, dims = dims, max.k = max.k)

DefaultAssay(object = dataset) <- "RNA"
# Local structure preservation
ls <- LocalStruct(object = dataset, grouping.var = "replicate", reduction = reduction, reduced.dims = dims, orig.dims = 1:30)
ls <- unname(obj = unlist(x = ls))

all.metrics <- list(
  silhouette = dataset$sil, 
  mixing.metric = mm,
  local.struct = ls
)


saveRDS(object = all.metrics, file = args[2])]

da <- dist(t(int_sub@assays$integrated@data))


plot(sil, main ="Silhouette plot - K-means")

install.packages("FactoMineR")

install.packages("factoextra")
library("FactoMineR")
library("factoextra")


res = 'integrated_snn_res.0.65'
silh <- function(res)
{
    Idents(integrated8.rpca) <- 'integrated_snn_res.0.75'
    int_sub = subset(integrated8.rpca, downsample = 100)
    da <- as.matrix(dist(t(int_sub@assays$integrated@data), diag = TRUE, upper =TRUE))
    int_sub$clustering <- Idents(int_sub)
    si <- silhouette(as.integer(int_sub$clustering), dmatrix = da)
    return (mean(si[,3]))
}


for (gene in c('RGS1', 'ITGAE', 'CD69', 'CXCR6', 'CCR7', 'CXCR5', 'EOMES', 'TBX21', 'FOXP3', 'IL22', 'IL17A', 'PRF1', 'GZMA', 'CD28'))
{
    pdf(paste("~/Desktop/int/TREE_",gene,".pdf",sep = ''), width = 16, height = 9)
    clustree(integrated8.rpca, prefix = "integrated_snn_res.", node_colour = 'CD28',node_colour_aggr = 'median')
    dev.off()
}


pdf(paste("~/Desktop/int/TREE.pdf"), width = 16, height = 12)
clustree(integrated8.rpca,node_colour_aggr = 'median')
dev.off()


integrated8.rpca <- FindNeighbors(integrated8.rpca, dims = 1:30)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.1)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.2)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.3)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.4)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.5)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.6)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.7)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.8)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.9)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 1)

integrated8.rpca <- RunUMAP(integrated8.rpca, dims = 1:30)

clustree(integrated8.rpca,node_colour_aggr = 'median', prefix = "groups", node_colour = 'CD28')


Idents(integrated8.rpca) <- 'orig.ident'
write.csv(mj019@assays$RNA@data,'~/Desktop/int/data.csv')

write.csv(integrated8.rpca$clonotype,'~/Desktop/int/Clone.csv')
write.csv(integrated8.rpca$orig.ident,'~/Desktop/int/Cell.csv')
write.csv(integrated8.rpca$pre_post, '~/Desktop/int/Detail.csv')
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.5)
write.csv(integrated8.rpca$seurat_clusters, '~/Desktop/int/Cluster.csv')


mj018 <- subset(integrated8.rpca, idents='Pt21_POD1145_IEL')
write.csv(mj018$orig.ident,'~/Desktop/int/MJ018_cell.csv')
write.csv(mj018@assays$RNA@data,'~/Desktop/int/MJ018_data.csv')
write.csv(mj018$pre_post, '~/Desktop/int/MJ018_detail.csv')
write.csv(mj018$seurat_clusters, '~/Desktop/int/MJ18_cluster.csv')


mj019 <- subset(integrated8.rpca, idents='Pt21_POD1145_LPL')
write.csv(mj019$orig.ident,'~/Desktop/int/MJ019_cell.csv')
write.csv(mj019@assays$RNA@data,'~/Desktop/int/MJ019_data.csv')
write.csv(mj019$pre_post, '~/Desktop/int/MJ019_detail.csv')
write.csv(mj019$seurat_clusters, '~/Desktop/int/MJ19_cluster.csv')


anchors<- FindTransferAnchors(reference = bm,query = integrated8.rpca, k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:64
  )


integrated8.rpca<- MapQuery(
    anchorset = anchors, 
    query = integrated8.rpca,
    reference = bm, 
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2", 
      predicted_ADT = "ADT"),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )


################################################################################



round(table(immune.combined2$predicted.celltype[immune.combined2$orig.ident == tis])/sum(table(immune.combined2$predicted.celltype))*100,2)



table(integrated8.rpca$orig.ident)

filter(integrated8.rpca, integrated8.rpca$orig.ident == 'D251' & integrated8.rpca$groups == 'c01') 



table(integrated8.rpca$groups[integrated8.rpca$orig.ident == 'D251'])



for (sample in names(table(integrated8.rpca$orig.ident)))
{
    print (sample)
    print (table(integrated8.rpca$groups[integrated8.rpca$orig.ident == sample]))
}


'D251',            'Pt14_POD1764',     'Pt15_POD1194',     'Pt21_POD626',
'Pt04_POD1606_IEL','Pt13_POD1032_IEL', 'Pt16_POD1004_IEL', 'Pt21_POD1145_IEL',
'Pt04_POD1606_LPL','Pt13_POD1032_LPL', 'Pt16_POD1004_LPL', 'Pt21_POD1145_LPL'


'Pt14_POD1764','Pt04_POD1606_IEL','Pt04_POD1606_LPL','Pt21_POD1145_IEL','Pt21_POD1145_LPL'

'Pt13_POD1032_IEL','Pt16_POD1004_LPL','Pt15_POD1194','Pt16_POD1004_IEL','Pt16_POD1004_LPL','Pt21_POD626'

mj019 <- subset(integrated8.rpca, idents='Pt21_POD1145_LPL')

Idents(integrated8.rpca) <- 'orig.ident'


chr <- subset(integrated8.rpca, idents=c('Pt04_POD1606_IEL','Pt04_POD1606_LPL','Pt14_POD1764','Pt21_POD1145_IEL','Pt21_POD1145_LPL'))
qui <- subset(integrated8.rpca, idents=c('Pt13_POD1032_IEL','Pt13_POD1032_LPL','Pt15_POD1194','Pt16_POD1004_IEL','Pt16_POD1004_LPL','Pt21_POD626'))


Idents(chr) <- 'cate'
for (sample in c("Persistent HvG","De novo H'vG","Acquired H'vG","Missing HvG","Persistent nonHvG"))
{
    print(sample)
    mark <- FindMarkers(chr, ident.1 = "Tolerant HvG", ident.2 = sample, assay = 'SCT', slot = "data",logfc.threshold = 0)
    write.csv(mark, paste('Dropbox/Jianing/scRNA-seq TRM paper/306/chronic rejection/SCT/Tolerant HvG vs. ',sample,'.csv',sep=''))
}



Idents(qui) <- 'cate'
for (sample in c("Persistent HvG","De novo H'vG","Missing HvG","Persistent nonHvG"))
{
    print(sample)
    mark <- FindMarkers(qui, ident.1 = "Tolerant HvG", ident.2 = sample, assay = 'SCT', slot = "data",logfc.threshold = 0)
    write.csv(mark, paste('Dropbox/Jianing/scRNA-seq TRM paper/306/quiescent/SCT/Tolerant HvG vs. ',sample,'.csv',sep=''))
}



integrated8.rpca$catcq <- integrated8.rpca$categories
for (sample in c('Pt04_POD1606_IEL','Pt04_POD1606_LPL','Pt14_POD1764','Pt21_POD1145_IEL','Pt21_POD1145_LPL'))
{
    print (sample)
    integrated8.rpca$catcq[(integrated8.rpca$catcq == "Persistent HvG") & (integrated8.rpca$orig.ident == sample)] <- "+ Persistent HvG"
    integrated8.rpca$catcq[(integrated8.rpca$catcq == "Tolerant HvG")& (integrated8.rpca$orig.ident == sample)] <- "+ Tolerant HvG"
    integrated8.rpca$catcq[(integrated8.rpca$catcq == "Missing HvG")& (integrated8.rpca$orig.ident == sample)] <- "+ Missing HvG"
    integrated8.rpca$catcq[(integrated8.rpca$catcq == "Acquired H'vG")& (integrated8.rpca$orig.ident == sample)] <- "+ Acquired H'vG"
    integrated8.rpca$catcq[(integrated8.rpca$catcq == "De novo H'vG")& (integrated8.rpca$orig.ident == sample)] <- "+ De novo H'vG"
    integrated8.rpca$catcq[(integrated8.rpca$catcq == "Persistent nonHvG")& (integrated8.rpca$orig.ident == sample)] <- "+ Persistent nonHvG"
}

#integrated8.rpca$catcq <- integrated8.rpca$categories
for (sample in c('Pt13_POD1032_IEL','Pt13_POD1032_LPL','Pt15_POD1194','Pt16_POD1004_IEL','Pt16_POD1004_LPL','Pt21_POD626'))
{
    print (sample)
    integrated8.rpca$catcq[(integrated8.rpca$catcq == "Persistent HvG")& (integrated8.rpca$orig.ident == sample)] <- "- Persistent HvG"
    integrated8.rpca$catcq[(integrated8.rpca$catcq == "Tolerant HvG")& (integrated8.rpca$orig.ident == sample)] <- "- Tolerant HvG"
    integrated8.rpca$catcq[(integrated8.rpca$catcq == "Missing HvG")& (integrated8.rpca$orig.ident == sample)] <- "- Missing HvG"
    integrated8.rpca$catcq[(integrated8.rpca$catcq == "Acquired H'vG")& (integrated8.rpca$orig.ident == sample)] <- "- Acquired H'vG"
    integrated8.rpca$catcq[(integrated8.rpca$catcq == "De novo H'vG")& (integrated8.rpca$orig.ident == sample)] <- "- De novo H'vG"
    integrated8.rpca$catcq[(integrated8.rpca$catcq == "Persistent nonHvG")& (integrated8.rpca$orig.ident == sample)] <- "- Persistent nonHvG"
}


Idents(integrated8.rpca) <- 'catcq'
for (sample in c("Persistent HvG","De novo H'vG","Missing HvG","Persistent nonHvG","Tolerant HvG"))
{
    s1 = paste("+",sample)
    s2 = paste("-",sample)
    print (s1)
    mark <- FindMarkers(integrated8.rpca, ident.1 = s1, ident.2 = s2, assay = 'SCT', slot = "data",logfc.threshold = 0)
    write.csv(mark, paste('Dropbox/Jianing/scRNA-seq TRM paper/306/cvq/SCT/',sample,'.csv',sep=''))
}



################################################################################
mj019 <- subset(integrated8.rpca, idents='Pt21_POD1145_LPL')





Idents(integrated8.rpca) <- 'pre'

for (samples in names(table(integrated8.rpca$orig.ident)))
{
    output <- list()
    for (pres in c("CD4_HvG",  "CD8_HvG", "CD4_nonHvG", "CD8_nonHvG", "Un"))
    {
        temp <- table(integrated8.rpca$groups[(integrated8.rpca$pre == pres)& (integrated8.rpca$orig.ident == samples)])
        output <- append(output, temp)
    }
    write.csv(output, paste('~/Dropbox/Jianing/scRNA-seq TRM paper/306/cell/',samples,'_pre.csv',sep=''))
}


for (samples in names(table(integrated8.rpca$orig.ident)))
{
    output <- list()
    for (pres in c("CD4_H'vG", "CD8_H'vG", "CD4_nonH'vG", "CD8_nonH'vG", "Un"))
    {
        temp <- table(integrated8.rpca$groups[(integrated8.rpca$post == pres)& (integrated8.rpca$orig.ident == samples)])
        output <- append(output, temp)
    }
    write.csv(output, paste('~/Dropbox/Jianing/scRNA-seq TRM paper/306/cell/',samples,'_post.csv',sep=''))
}


for (samples in names(table(integrated8.rpca$orig.ident)))
{
    output <- list()
    for (pres in c("Persistent HvG", "Tolerant HvG", "Missing HvG", "Acquired H'vG", "De novo H'vG","Persistent nonHvG","Others"))
    {
        temp <- table(integrated8.rpca$groups[(integrated8.rpca$cate == pres)& (integrated8.rpca$orig.ident == samples)])
        output <- append(output, temp)
    }
    write.csv(output, paste('~/Dropbox/Jianing/scRNA-seq TRM paper/306/cell/',samples,'_com.csv',sep=''))
}



Idents(integrated8.rpca) <- 'post'
output <- list()
for (samples in c("Persistent nonHvG", "Tolerant HvG", "Missing HvG", "Acquired H'vG", "De novo H'vG","Persistent nonHvG","Others"))
{
    output <- append(output, table(integrated8.rpca$groups[integrated8.rpca$cate == samples]))
}



################################################################################
tolerant HVG cells
cluster group: 
c01, c02: TRM group; c03,c04,c07: Teff/TRM group; c05: nonTRM group; c06,c08: Tfh group; c10: Treg group.
quiescent vs chronic rejection

which(data$V1 > 2 | data$V2 < 4) , ]


thc <- subset(integrated8.rpca, idents='Tolerant HvG')


thc$catcq <- thc$categories
for (sample in c('Pt04_POD1606_IEL','Pt04_POD1606_LPL','Pt14_POD1764','Pt21_POD1145_IEL','Pt21_POD1145_LPL'))
{
    print (sample)
    thc$catcq[(thc$catcq == "Tolerant HvG")& (thc$orig.ident == sample)] <- "chronic rejection Tolerant HvG"
}


for (sample in c('Pt13_POD1032_IEL','Pt13_POD1032_LPL','Pt15_POD1194','Pt16_POD1004_IEL','Pt16_POD1004_LPL','Pt21_POD626'))
{
    thc$catcq[(thc$catcq == "Tolerant HvG")& (thc$orig.ident == sample)] <- "quiescent Tolerant HvG"
}


Idents(thc) <- 'groups'
trmg <-  subset(thc, idents= c('c01', 'c02'))
tefg <-  subset(thc, idents= c('c03', 'c04', 'c07'))
ntrm <-  subset(thc, idents= c('c05'))
tfhg <-  subset(thc, idents= c('c06', 'c08'))
treg <-  subset(thc, idents= c('c10'))


deparse(substitute(sample))
deparse(substitute(sample))

Idents(integrated8.rpca) <- 'catcq'

for (sample in list(trmg, tefg, ntrm, tfhg, treg))
{
    print (deparse(substitute(sample)))
    mark <- FindMarkers(sample, ident.1 = 'chronic rejection Tolerant HvG', ident.2 = 'quiescent Tolerant HvG', assay = 'SCT', slot = "data",logfc.threshold = 0)
    write.csv(mark, paste('Dropbox/Jianing/scRNA-seq TRM paper/306/cvq/SCT/',sample,'.csv',sep=''))
}


Idents(trmg) <- 'catcq'
mark <- FindMarkers(trmg, ident.1 = 'chronic rejection Tolerant HvG', ident.2 = 'quiescent Tolerant HvG', assay = 'SCT', slot = "data",logfc.threshold = 0)
write.csv(mark, 'Dropbox/Jianing/scRNA-seq TRM paper/306/cvq/SCT/trmg.csv')

Idents(tefg) <- 'catcq'
mark <- FindMarkers(tefg, ident.1 = 'chronic rejection Tolerant HvG', ident.2 = 'quiescent Tolerant HvG', assay = 'SCT', slot = "data",logfc.threshold = 0)
write.csv(mark, 'Dropbox/Jianing/scRNA-seq TRM paper/306/cvq/SCT/tefg.csv')


Idents(ntrm) <- 'catcq'
mark <- FindMarkers(ntrm, ident.1 = 'chronic rejection Tolerant HvG', ident.2 = 'quiescent Tolerant HvG', assay = 'SCT', slot = "data",logfc.threshold = 0)
write.csv(mark, 'Dropbox/Jianing/scRNA-seq TRM paper/306/cvq/SCT/ntrm.csv')

Idents(tfhg) <- 'catcq'
mark <- FindMarkers(tfhg, ident.1 = 'chronic rejection Tolerant HvG', ident.2 = 'quiescent Tolerant HvG', assay = 'SCT', slot = "data",logfc.threshold = 0)
write.csv(mark, 'Dropbox/Jianing/scRNA-seq TRM paper/306/cvq/SCT/tfhg.csv')

Idents(treg) <- 'catcq'
mark <- FindMarkers(treg, ident.1 = 'chronic rejection Tolerant HvG', ident.2 = 'quiescent Tolerant HvG', assay = 'SCT', slot = "data",logfc.threshold = 0)
write.csv(mark, 'Dropbox/Jianing/scRNA-seq TRM paper/306/cvq/SCT/treg.csv')



################################################################################




thc$regroup <- thc$groups
thc$regroup[which(thc$groups == 'c01' | thc$groups == 'c02')] <- 'TRM group'
thc$catcq <- thc$categories
thc$catcq[(thc$catcq == "Persistent HvG")& (thc$orig.ident == sample)] <- "- Persistent HvG"




Idents(integrated8.rpca) <- 'orig.ident'
integrated8.rpca <- FindNeighbors(integrated8.rpca, dims = 1:32)


DefaultAssay(integrated8.rpca) <- 'integrated'
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0)
#integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.05,reduction = "umap")
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.1)
#integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.15)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.2)
#integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.25)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.3)
#integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.35)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.48)
#integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.45)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.5)
#integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.55)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.6)
#integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.65)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.7)
#integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.75)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.8)
#integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.85)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.9)
#integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 0.95)
integrated8.rpca <- FindClusters(integrated8.rpca, resolution = 1)


ncluster <- integrated8.rpca$integrated_snn_res.0.48
recluster <- function(ncluster){
    ncluster <- as.character(ncluster)
    for(i in seq(length(ncluster))){
        ncluster[i] = paste('c',sprintf("%02d", as.numeric(ncluster[i])+1),sep='')
    }
    ncluster <- as.factor(ncluster)
    return(ncluster)
}
integrated8.rpca$groups.0.55 <- recluster(ncluster)


clustree(integrated8.rpca, prefix = "groups.", node_colour = "purple", node_size = 10,
         node_alpha = 0.8)

'FOXP3', 'CCR7', 
for (gene in c('RGS1', 'ITGAE', 'CD69', 'CXCR6', 'CCR7', 'CXCR5', 'EOMES', 'TBX21', 'FOXP3', 'IL22', 'IL17A', 'PRF1', 'GZMA', 'CD28'))
for (gene in c('CD69', 'RGS1', 'TBX21', 'EOMES', 'CD28', 'SELL', 'KLF2', 'MKI67'))
{
    pdf(paste("~/Desktop/int/clustree/TREE_",gene,".pdf",sep = ''), width = 16, height = 9)
    clustree(integrated8.rpca, prefix = "groups.", node_colour = gene,node_colour_aggr = 'median')
    dev.off()
}



pdf("~/Desktop/int/clustree/TREE_FOXP3.pdf", width = 16, height = 9)
clustree(integrated8.rpca, prefix = "groups.", node_colour = 'FOXP3',node_colour_aggr = 'median',node_text_colour = "red",node_alpha = 0.8)+ scale_fill_distiller(palette = "RdPu")
dev.off()





small_subset<- RandomSubsetData(integrated8.rpca, 0.8)

small_subset <- RunPCA(small_subset, verbose = TRUE)
small_subset <- RunUMAP(small_subset, dims = 1:30)
small_subset <- FindNeighbors(small_subset,  dims = 1:30)
small_subset <- FindClusters(small_subset, resolution = 0.5)


setNames(data.frame(matrix(ncol = 15, nrow = 0)), )


l <- data.frame(matrix(ncol=15,nrow=0, dimnames=list(NULL, names(table(integrated8.rpca$groups)))))
ljc <- data.frame(matrix(ncol=2,nrow=0, dimnames=list(NULL, c('Jaccard', 'Cluster'))))

for (k in seq(10)){
    print (k)
    for (i in seq(0.25,0.75,0.05)){
        print (i)
        small_subset<- RandomSubsetData(integrated8.rpca, 0.8)
        small_subset <- RunPCA(small_subset, verbose = FALSE)
        small_subset <- RunUMAP(small_subset, dims = 1:30, verbose = FALSE)
        small_subset <- FindNeighbors(small_subset,  dims = 1:30, verbose = FALSE)
        small_subset <- FindClusters(small_subset, resolution = i, verbose = TRUE)
        jd <- PairWiseJaccardSets(integrated8.rpca$groups, small_subset$seurat_clusters)
        small_subset <- NULL
        gc()
        for (j in seq(15)){
            ljc[nrow(ljc) + 1,] <- c(max(jd[j,]), names(table(integrated8.rpca$groups))[j])}
    }
}


ljc$Cluster <- as.factor(ljc$Cluster)
ljc$Jaccard <- as.numeric(ljc$Jaccard)

head(ljc)


pdf("~/Desktop/int/clustree/violin.pdf", width = 8, height = 6)
ggplot(ljc, aes(x=Cluster, y=Jaccard, fill = Cluster)) + geom_violin(scale = "width",draw_quantiles = c(0.25, 0.5, 0.75),alpha = 0.8) + scale_fill_manual(values=colors)
dev.off()



################################################################################



integrated8_top20 <- integrated8.rpca.cat.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
jpeg(paste("~/Desktop/int/HEAT_cat_makers.jpeg",sep=''), width = 1600*3, height = 2000*3,res=300)
DoHeatmap(integrated8.rpca.cat, features = integrated8_top20$gene, group.colors =c("#00BCD4","#3F51B5","#4CAF50","#E91E63","#9C27B0","#FF9800"),group.by='categories',assay='RNA')
dev.off()



gene2 <- c('CD69','PRDM1','RGS1','IL22','ICOS','IL21','IFNG','TNF','GZMA','CXCR6','PRDM1','IL22','RORC','BATF','EOMES','TNF','GZMB','ITGA1','RUNX3','ZNF683','TBX21','TRDC','GNLY','GZMA','GZMB','NKG7','RUNX3','EOMES','KLRK1','TRDC','CD8A','RUNX3','ZNF683','TBX21','EOMES','KLF2','GNLY','GZMB','NKG7','PRF1','CD8A','PDCD1','CCR7','KLF2','S1PR1','SELL','TCF7','CTLA4','ICOS','BCL6','CXCR5','MAF','EOMES','CD28','KLRK1','TCF7','CTLA4','PDCD1','IL21','BATF','CXCL13','CXCR5','MAF','CCR7','TNF','CD4','CXCR6','RGS1','CTLA4','FOXP3','ICOS','BATF','CD28')
gene3 <- c('CD4','CD8A','TRDC','CD69','RGS1','CXCR6','ITGA1','RORC','RUNX3','IL22','IFNG','TNF','GZMA','GZMB','PRF1','NKG7','KLRK1','EOMES','TBX21','ZNF683','KLF2','S1PR1','SELL','CCR7','TCF7','BCL6','MAF','CXCR5','PDCD1','CXCL13','CTLA4','ICOS','CD28','IL21','PRDM1','BATF','FOXP3')

integrated8.rpca$regroup <- as.character(integrated8.rpca$groups)
integrated8.rpca$regroup[integrated8.rpca$regroup %in% c('c01', 'c02')] <- "01.TRM"
integrated8.rpca$regroup[integrated8.rpca$regroup %in% c('c03', 'c04', 'c07')] <- "02.Teff/TRM"
integrated8.rpca$regroup[integrated8.rpca$regroup %in% c('c05')] <- "03.non TRM"
integrated8.rpca$regroup[integrated8.rpca$regroup %in% c('c06', 'c08')] <- "04.Tfh"
integrated8.rpca$regroup[integrated8.rpca$regroup %in% c('c10')] <- "05.Treg"

Idents(integrated8.rpca) <- 'groups'
integrated8.rpca <- subset(integrated8.rpca, idents=c('c01', 'c02','c03', 'c04', 'c07','c05','c06', 'c08','c10'))

integrated8.rpca$regroups <- as.factor(integrated8.rpca$regroup)
Idents(integrated8.rpca) <- 'regroups'

integrated8.rpca <- AverageExpression(integrated8.rpca, return.seurat=TRUE) 

jpeg(paste("~/Desktop/int/HEAT_gene3.jpeg",sep=''), width = 600*2, height = 900*2,res=300)
DoHeatmap(integrated8.rpca, features = gene3, group.colors =c('#E91E63','#3F51B5','#00BCD4','#8BC34A','#FFC107'),draw.lines = FALSE,angle = 30)
dev.off()


DefaultAssay
DoHeatmap(integrated8.rpca, features = gene2, group.colors =c("#00BCD4","#3F51B5","#4CAF50","#E91E63","#9C27B0"),draw.lines = FALSE,angle = 270)

################################################################################

catc <- c("Persistent HvG","Tolerant HvG","Missing HvG","Acquired H'vG","De novo H'vG","Persistent nonHvG")
integrated8.rpca.cat <- subset(integrated8.rpca, idents=catc)

data1 <- cbind(integrated8.rpca.cat$orig.ident, integrated8.rpca.cat$groups, integrated8.rpca.cat$cate) 

colnames(data1) <- c("col1", "col2") 


integrated8.rpca.cat$regroup <- as.character(integrated8.rpca.cat$groups)
integrated8.rpca.cat$regroup[integrated8.rpca.cat$regroup %in% c('c01', 'c02')] <- "01.TRM"
integrated8.rpca.cat$regroup[integrated8.rpca.cat$regroup %in% c('c03', 'c04', 'c07')] <- "02.Teff/TRM"
integrated8.rpca.cat$regroup[integrated8.rpca.cat$regroup %in% c('c05')] <- "03.non TRM"
integrated8.rpca.cat$regroup[integrated8.rpca.cat$regroup %in% c('c06', 'c08')] <- "04.Tfh"
integrated8.rpca.cat$regroup[integrated8.rpca.cat$regroup %in% c('c10')] <- "05.Treg"

data2 <- cbind(sample_ID = as.character(integrated8.rpca.cat$repos), cluster_ID = as.character(integrated8.rpca.cat$groups)
, cell_annotation = as.character(integrated8.rpca.cat$regroup), categories = as.character(integrated8.rpca.cat$cate))
write.csv(data2,'~/Desktop/categories.csv')


for (samples in c("Persistent nonHvG", "Tolerant HvG", "Missing HvG", "Acquired H'vG", "De novo H'vG","Persistent nonHvG","Others"))
{
    Idents(integrated8.rpca.cat) <- 'cate'
    subdata <- subset(integrated8.rpca.cat, idents = samples)
    data3 <- cbind(sample_ID = as.character(subdata$repos), cluster_ID = as.character(subdata$groups), cell_annotation = as.character(subdata$regroup), categories = as.character(subdata$cate))
    write.csv(data3,paste('~/Desktop/',samples,'data.csv'))
}



################################################################################
#treg cluster's acquired HVG gene expression by DE gene analysis
#persistent nonHVG in Treg cluster vs acquired H’vG in Treg cluster




Idents(integrated8.rpca) <- 'regroup'
subdata <- subset(integrated8.rpca, idents = '05.Treg')


Idents(subdata) <- 'categories'
mark <- FindMarkers(subdata, ident.1 = "Persistent nonHvG", ident.2 = "Acquired H'vG", assay = 'SCT', slot = "data",logfc.threshold = 0)
write.csv(mark, 'Dropbox/Jianing/scRNA-seq TRM paper/1003/Per_Acq_treg.csv')



################################################################################

#violin plot shows CD28, NKD2D expression in quiescent vs rejecting samples within CD8A positive cells in combined Teff/Trm cluster (c03, c04, c07)


tolerant HVG cells
cluster group: 
c01, c02: TRM group; c03,c04,c07: Teff/TRM group; c05: nonTRM group; c06,c08: Tfh group; c10: Treg group.
quiescent vs chronic rejection

which(data$V1 > 2 | data$V2 < 4) , ]


thc <- subset(integrated8.rpca, idents='Tolerant HvG')


thc$catcq <- thc$categories
for (sample in c('Pt04_POD1606_IEL','Pt04_POD1606_LPL','Pt14_POD1764','Pt21_POD1145_IEL','Pt21_POD1145_LPL'))
{
    print (sample)
    thc$catcq[(thc$catcq == "Tolerant HvG")& (thc$orig.ident == sample)] <- "chronic rejection Tolerant HvG"
}


for (sample in c('Pt13_POD1032_IEL','Pt13_POD1032_LPL','Pt15_POD1194','Pt16_POD1004_IEL','Pt16_POD1004_LPL','Pt21_POD626'))
{
    thc$catcq[(thc$catcq == "Tolerant HvG")& (thc$orig.ident == sample)] <- "quiescent Tolerant HvG"
}


Idents(thc) <- 'groups'
trmg <-  subset(thc, idents= c('c01', 'c02'))
tefg <-  subset(thc, idents= c('c03', 'c04', 'c07'))
ntrm <-  subset(thc, idents= c('c05'))
tfhg <-  subset(thc, idents= c('c06', 'c08'))
treg <-  subset(thc, idents= c('c10'))


deparse(substitute(sample))
deparse(substitute(sample))

Idents(integrated8.rpca) <- 'catcq'

for (sample in list(trmg, tefg, ntrm, tfhg, treg))
{
    print (deparse(substitute(sample)))
    mark <- FindMarkers(sample, ident.1 = 'chronic rejection Tolerant HvG', ident.2 = 'quiescent Tolerant HvG', assay = 'SCT', slot = "data",logfc.threshold = 0)
    write.csv(mark, paste('Dropbox/Jianing/scRNA-seq TRM paper/306/cvq/SCT/',sample,'.csv',sep=''))
}


chr <- subset(integrated8.rpca, idents=c('Pt04_POD1606_IEL','Pt04_POD1606_LPL','Pt14_POD1764','Pt21_POD1145_IEL','Pt21_POD1145_LPL'))
qui <- subset(integrated8.rpca, idents=c('Pt13_POD1032_IEL','Pt13_POD1032_LPL','Pt15_POD1194','Pt16_POD1004_IEL','Pt16_POD1004_LPL','Pt21_POD626'))


Idents(integrated8.rpca) <- 'groups'
tefg <-  subset(integrated8.rpca, idents= c('c03', 'c04', 'c07'))
tefg_cd8 <- subset(tefg , cells = WhichCells(tefg, slot = 'counts', expression =CD8A>0))


tefg_cd8$qvr <-'Others'

for (sample in c('Pt04_POD1606_IEL','Pt04_POD1606_LPL','Pt14_POD1764','Pt21_POD1145_IEL','Pt21_POD1145_LPL'))
{
    print (sample)
    tefg_cd8$qvr[tefg_cd8$orig.ident == sample] <- "Rejecting"
}


for (sample in c('Pt13_POD1032_IEL','Pt13_POD1032_LPL','Pt15_POD1194','Pt16_POD1004_IEL','Pt16_POD1004_LPL','Pt21_POD626'))
{
    tefg_cd8$qvr[tefg_cd8$orig.ident == sample] <- "Quiescent"
}


Idents(tefg_cd8) <- 'qvr'
tefg_cd8 <-  subset(tefg_cd8, idents= c('Quiescent','Rejecting'))

VlnPlot(object = tefg_cd8, features = 'CD28',slot = "counts", log = TRUE)

 

data2 <- cbind(sample_ID = as.character(tefg_cd8_gene$repos), cluster_ID = as.character(tefg_cd8_gene$groups)
, cell_annotation = as.character(tefg_cd8_gene$qvr), CD28 = as.character(tefg_cd8_cd28@assays$integrated@data), 
NKG2D = as.character(tefg_cd8_klrk1@assays$integrated@data))
write.csv(data2,'~/Desktop/cd28_nkg2d.csv')

jpeg("~/Desktop/cd28_nkg2d_vln.jpeg",width = 800, height = 600)
VlnPlot(object = tefg_cd8, features = c('CD28','KLRK1'),assay='integrated',slot = "data")
dev.off()

tefg_cd8_cd28 <- subset(tefg_cd8, features = "CD28")
tefg_cd8_klrk1 <- subset(tefg_cd8,features = "KLRK1")



as.character(tefg_cd8_gene@assays$integrated@data)



################################Trajectory######################################

library(monocle3)
library(SeuratWrappers)
library(patchwork)
integrated8.rpca <- readRDS('/Volumes/backup/integrated8.rpca0919.rds')
cds <- as.cell_data_set(integrated8.rpca)
rm(integrated8.rpca)
cds <- cluster_cells(cds, resolution=1e-4)
p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)
integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)
rm(integrated.sub)


cds <- learn_graph(cds, use_partition = TRUE, verbose = FALSE)

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 4]))

pdf("~/Desktop/int/UMAP.pdf", width = 8, height = 4.5)
cds <- cluster_cells(cds, resolution=1e-5)


p0 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 1]))
p1 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey40")

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 2]))
p2 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey40")

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 3]))
p3 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey40")


cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 4]))
p4 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey40")


colors <- c('#E91E63','#3F51B5','#00BCD4','#8BC34A','#FFC107','#795548',
            '#9C27B0','#2196F3','#009688','#CDDC39','#FF9800','#9E9E9E',
            '#673AB7','#03A9F4','#4CAF50','#FFEB3B','#FF5722','#607D8B',
            '#F44336','#000000',
            '#F48FB1','#9FA8DA','#80DEEA','#C5E1A5','#FFE082','#BCAAA4',
            '#CE93D8','#90CAF9','#80CBC4','#E6EE9C','#FFCC80','#EEEEEE'
            )

pm1 <- DimPlot(integrated.sub ,label.size = 4,label = TRUE, cols = colors, group.by = 'groups')+NoLegend()


jpeg("~/Desktop/UMAP_Trajectory1.jpeg", width = 1000, height = 400)
wrap_plots(pm1,p0, p1, p2, p3, p4)
dev.off()





p4_0 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 2]))
p4_2 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey40")

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 6]))
p4_6 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey40")

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 7]))
p4_7 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey40")


cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 9]))
p4_9 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey40")


cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 14]))
p4_14 <- plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,å
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey40")

jpeg("~/Desktop/UMAP_Trajectory_4.jpeg", width = 1000, height = 400)
wrap_plots(p4_0,p4_2,p4_6,p4_7,p4_9,p4_14)
dev.off()



cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)
rowData(cds)$gene_short_name <- row.names(rowData(cds))
head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)
deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))


jpeg("~/Desktop/genes_change.jpeg", width = 800, height = 400)
plot_cells(cds,
           genes=head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)
dev.off()


gene_modules <- find_gene_modules(cds[deg_ids,], resolution=c(10^seq(-6,-1)))
table(gene_modules$module)


################################Velocity######################################

integrated.loom <- as.loom(integrated.sub, filename = "~/Desktop/integrated.loom", verbose = TRUE)
integrated <- as.Seurat(x = integrated.loom)
DefaultAssay(object = integrated) <- "spliced"
integrated <- RunVelocity(object = integrated, deltaT = 1, kCells = 25, fit.quantile = 0.02)




library(SCP)


integrated.velo <- RunSCVELO(
  srt = integrated.sub, group_by = "groups",
  linear_reduction = "PCA", nonlinear_reduction = "UMAP"
)
VelocityPlot(srt = integrated.velo, reduction = "UMAP", group_by = "SubCellType")


################################################################################
integrated8.rpca <- readRDS('/Volumes/backup/integrated8.rpca0919.rds')
integrated8.rpca$barcode <- colnames(integrated8.rpca)
integrated8.rpca$UMAP_1 <- integrated8.rpca@reductions$umap@cell.embeddings[,1]
integrated8.rpca$UMAP_2 <- integrated8.rpca@reductions$umap@cell.embeddings[,2]
write.csv(integrated8.rpca@meta.data, file='~/Desktop/traj/metadata.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(integrated8.rpca, assay='RNA', slot='counts')
writeMM(counts_matrix, file='~/Desktop/traj/counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(integrated8.rpca@reductions$pca@cell.embeddings, file='~/Desktop/traj/pca.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='~/Desktop/traj/gene_names.csv',
  quote=F,row.names=F,col.names=F
)

######

velocyto run ~/Documents/workstation/data/jian/MJ001/ ~/Desktop/traj/Homo_sapiens.GRCh38.108.gtf 


