#设置工作路径
getwd()
setwd("F:\\scRNA-seq")
rm(list=ls())
options(stringsAsFactors = F)


#单细胞转录组基础分析一：分析环境搭建
#install.packages("devtools", dependencies=T)
#install.packages("BiocManager", dependencies=T)
#install.packages("tidyverse", dependencies=T)
#install.packages('Seurat', dependencies=T)
#BiocManager::install(c("SingleR","monocle", "DESeq2"),ask = F,update = F) 
#BiocManager::install(c("clusterProfiler","DOSE","pheatmap"),ask = F,update = F)
#BiocManager::install(c("org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db"),ask = F,update = F)
#devtools::install_github('RGLab/MAST', upgrade=F, build_vignettes = T)

#单细胞转录组基础分析二：数据质控与标准化
library(Seurat)
library(tidyverse)
rm(list=ls())
dir.create("QC")
##单细胞测序生信分析1-Seurat
https://www.jianshu.com/p/e7303b4db9b8
http://www.bio-info-trainee.com/6392.html
##构建Seurat对象
#min.cells:基因最少在3个细胞中表达，min.features：细胞中最少有200个基因表达
GSM4143678<- Read10X(data.dir = "E:/Westlake_University/Bioinformatics/SingleCellAnalysis/Nature/GSM4143678")
pbmc1 <- CreateSeuratObject(counts = GSM4143678,
                            min.cells = 3, 
                            min.features = 200)
head(pbmc1@meta.data)
GSM4143679<- Read10X(data.dir = "E:\\Westlake_University\\Bioinformatics\\SingleCellAnalysis\\Nature\\GSM4143679")
pbmc2 <- CreateSeuratObject(counts = GSM4143679,
                            min.cells = 3, 
                            min.features = 200)
head(pbmc2@meta.data)
GSM4143681<- Read10X(data.dir = "E:\\Westlake_University\\Bioinformatics\\SingleCellAnalysis\\Nature\\GSM4143681")
pbmc3 <- CreateSeuratObject(counts =GSM4143681,
                            min.cells = 3, 
                            min.features = 200)
head(pbmc3@meta.data)
GSM4143682<- Read10X(data.dir = "E:\\Westlake_University\\Bioinformatics\\SingleCellAnalysis\\Nature\\GSM4143682")
pbmc4 <- CreateSeuratObject(counts = GSM4143682,
                            min.cells = 3, 
                            min.features = 200)
head(pbmc4@meta.data)
GSM4143684<- Read10X(data.dir = "E:\\Westlake_University\\Bioinformatics\\SingleCellAnalysis\\Nature\\GSM4143684")
pbmc5 <- CreateSeuratObject(counts = GSM4143684,
                            min.cells = 3, 
                            min.features = 200)
head(pbmc5@meta.data)
GSM4143685<- Read10X(data.dir = "E:\\Westlake_University\\Bioinformatics\\SingleCellAnalysis\\Nature\\GSM4143685")
pbmc6 <- CreateSeuratObject(counts = GSM4143685,
                            min.cells = 3, 
                            min.features = 200)
head(pbmc6@meta.data)
#根据每个pbmc中有多少个细胞进行标记
pbmc1@meta.data$group <- rep("RT1",7363)
pbmc2@meta.data$group <- rep("RN1",5879)
pbmc3@meta.data$group <- rep("RT2",5529)
pbmc4@meta.data$group <- rep("RN2",7976)
pbmc5@meta.data$group <- rep("RT3",6712)
pbmc6@meta.data$group <- rep("RN3",5538)
table(pbmc1@meta.data$orig.ident)         #查看样本的细胞数量
#可以随意合并其中的几个分组
pbmcX = merge(pbmc1, pbmc2, 
             add.cell.ids = c("GSM4143678", "GSM4143679"),
             merge.data = TRUE)
			 
pbmc = merge(pbmc1, y=c(pbmc2, pbmc3, pbmc4, pbmc5, pbmc6), 
             add.cell.ids = c("GSM4143678", "GSM4143679", "GSM4143681", "GSM4143682", "GSM4143684", "GSM4143685"),
             merge.data = TRUE)
as.data.frame(pbmc@assays$RNA@counts[1:10, 1:2])
head(pbmc@meta.data)
saveRDS(pbmc,file="pbmc_raw.rds")  #用于后续分析
save(pbmc,file = 'pbmc_raw.Rdata')