rt=read.table("merge.txt",sep="\t",header=T,check.names=F)   
rt=read.table("exp.txt",sep="\t",header=T,check.names=F)   
rt2=read.table("OverlappedGenes.txt",sep="\t",header=F,check.names=F)   
View(rt2)
View(rt)
rownames(rt) <- rt[,1]
rownames(rt2) <- rt2[,1]
rt3 <- rt[rownames(rt2),]
rt3 <- rt3[,2:ncol(rt3)]
rt4 <- t(rt3)

conNum=100                                                        #normal组样品数目
treatNum=303                                                 #tumor组样品数目
grade=c(rep("control",conNum),rep("treat",treatNum))



##进行K-means聚类并分析结果 参考
#https://www.datalearner.com/blog/1051493902550927 
#https://blog.csdn.net/weixin_47580081/article/details/115018232
#https://blog.csdn.net/hfutxiaoguozhi/article/details/78828047

kmeans(x, centers, iter.max = 10, nstart = 1,algorithm = c("Hartigan-Wong", "Lloyd", "Forgy","MacQueen"), trace=FALSE)
#第一个参数是 x，它是kmeans的输入数据。只要是矩阵的数值数据就可以了。首先，我们使用R语言自带的数据集iris数据
#第二个参数是中心点选择 centers，它是中心点选择，可以有两种参数，第一种直接写一个数值，表示聚成几类，这种情况下它会自动随机选择初始中心点，R语言会选择数据集中的随机行作为初试中心点。还有一种方式是可以选择输入一个起始的中心点，那么kmeans就会根据选择的中心点开始聚类迭代。
#第三个参数是迭代次数 iter.max，即迭代次数，不写的话默认是10次。否则就是最大迭代次数。
#第四个参数是算法选择 algorithm，输入值是缩写的字符，算法有”Hartigan-Wong”, “Lloyd”, “Forgy”,”MacQueen”四种，注意Lloyd和Forgy是同一种算法的不同称谓。
#示例：km <- kmeans(iris[,1:4], 3)
library(cluster)
library(factoextra)
kmeans_input <- pheatmap_list_count_new
#确定最佳聚类数目
fviz_nbclust(kmeans_input, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)
#可以发现聚为四类最合适，当然这个没有绝对的，从指标上看，选择坡度变化不明显的点最为最佳聚类数目。
#设置随机数种子，保证实验的可重复进行
set.seed(123)

km <- kmeans(kmeans_input, 4, iter.max = 10, algorithm = c("Hartigan-Wong"))
km



# 结果输出
type <- km$cluster
type# 查看类别分布
table(type)  # 查看类别统计

centervec <- km$center
centervec # 查看簇中心点

# centervec对列取最大值
max <- apply(centervec,2,max) 
max

# centervec对列取最小值
min <- apply(centervec,2,min)
min

# 构建frame类型数据
df = data.frame(rbind(max,min,centervec))
df
write.csv(df,file="Kmenas_df.csv") ##输出Kmeans结果



# 绘制聚类散点图
# 如果本地没有factoextra和cluster包，先用以下命令进行下载
# install.packages('factoextra')
# install.packages('cluster')
library(cluster)
library(factoextra)
fviz_cluster(km, data = kmeans_input, repel = TRUE)
fviz_cluster(km, data = kmeans_input,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
             ellipse.type = "euclid",
             star.plot = TRUE, 
             repel = TRUE,
             ggtheme = theme_minimal()
)




#提取某个cluster里的全部基因
type <- km$cluster
type# 查看类别分布
table(type)  # 查看类别统计
cluster1 <- as.data.frame(type[type %in% 1]) # 查看cluste2统计
write.table(rownames(cluster1),file="cluster1_gene.txt",sep="\t",quote=F,col.names=F)   
cluster2 <- as.data.frame(type[type %in% 2]) # 查看cluste2统计
write.table(rownames(cluster2),file="cluster2_gene.txt",sep="\t",quote=F,col.names=F)   
cluster3 <- as.data.frame(type[type %in% 3]) # 查看cluster4统计
write.table(rownames(cluster3),file="cluster3_gene.txt",sep="\t",quote=F,col.names=F)   
cluster4 <- as.data.frame(type[type %in% 4]) # 查看cluster4统计
write.table(rownames(cluster4),file="cluster4_gene.txt",sep="\t",quote=F,col.names=F)   
cluster5 <- as.data.frame(type[type %in% 5]) # 查看cluster5统计
write.table(rownames(cluster5),file="cluster5_gene.txt",sep="\t",quote=F,col.names=F)
cluster7 <- as.data.frame(type[type %in% 7]) # 查看cluster4统计
write.table(rownames(cluster7),file="cluster7_gene.txt",sep="\t",quote=F,col.names=F)

cluster1_name <- rownames(cluster1)
cluster2_name <- rownames(cluster2)
cluster3_name <- rownames(cluster3)
cluster4_name <- rownames(cluster4)
cluster5_name <- rownames(cluster5)
cluster7_name <- rownames(cluster7)
cluster_genes <- c(cluster1_name,cluster2_name,cluster5_name,cluster7_name)
cluster_genes <- c(cluster3_name)
write.table(cluster_genes,file="cluster_genes.txt",sep="\t",quote=F,col.names=F)   
rt4 <- t(pheatmap_list_count_new[cluster_genes,])





##PCA Analysis  参考 https://blog.csdn.net/nikang3148/article/details/85246555
#导入输入矩阵
rt4<-t(pheatmap_list_count_new)

#prcomp函数
rt4.pca<-prcomp(rt4,scale=T,rank=7,retx=T) #相关矩阵分解 rank为PC的数目
#retx表四返回score，scale表示要标准化
summary(rt4.pca) #方差解释度
rt4.pca$sdev #特征值的开方
rt4.pca$rotation #特征向量，回归系数
rt4.pca$x #样本得分score

#princomp函数
rt4.pca<-princomp(rt4,cor=T,scores=T) 
#默认方差矩阵(cor=F),改为cor=T则结果与prcomp相同
summary(rt4.pca) #各主成份的SVD值以及相对方差
rt4.pca$loading #特征向量，回归系数
rt4.pca$score
screenplot(rt4.pca) #方差分布图
biplot(rt4.pca,scale=F) #碎石图,直接把x与rotation绘图，而不标准化

library(ggbiplot)
ggscreeplot(rt4.pca) #碎石图


ggbiplot(rt4.pca, obs.scale = 1, var.scale = 1,
         groups = grade, ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

  
library(ggfortify)

autoplot(rt4.pca,data = rt4.pca,col= 'grade',size=2,
         loadings =T,loadings.label = TRUE,
         frame = TRUE,frame.type='norm',
         label = TRUE, label.size = 3
)+  theme_classic()



###############################################就用这个
#导入输入矩阵
rt4<-t(pheatmap_list_count_new)
library(factoextra)
library(FactoMineR)
rt4.pca2 <- PCA(rt4,scale.unit = T,ncp=10,graph = T) #ncp是确定PC的数目
print(rt4.pca2)
summary(rt4.pca2)
get_eigenvalue(rt4.pca2)    ##特征值大于1 #标准化数据中特征值>1的变量解释能力较强  方差整体解释80%以上
write.table(as.data.frame(get_eigenvalue(rt4.pca2)), file="PCA_eigenvalue_variance.percent.txt",sep="\t", quote=FALSE)
##这时重新确定PC的数目
rt4.pca2 <- PCA(rt4,scale.unit = T,ncp=5,graph = T) #ncp是确定PC的数目


fviz_eig(rt4.pca2,addlabels = TRUE) #碎石图,展示方差解释度

#变量提取主要有get_pca_var()函数,输出变量在主成分投影上的坐标，变量与主成分PC的相关系数，相关系数的平方，变量对某一PC的相关贡献
#get_pca_var()等同于get_pca(element="var")
get_pca_var(rt4.pca2)#coord cor cos2 contribution
rt4.pca2$var #输出上述数据
get_pca_var(rt4.pca2)$coord
get_pca_var(rt4.pca2)$cos2

#变量坐标(coord)与相关性(cor)可视化
#coord是坐标(实际的loading)，与cor数值相同
#coord=eigen vector * stdev
#相关图中，靠近的变量表示正相关；对向的是负相关。
#箭头越远离远原点、越靠经圆圈表明PC对其的代表性高（相关性强）

fviz_pca_var(rt4.pca2) #变量相关性可视化图

#cos2代表不同主成分对变量的代表性强弱，对特定变量，所有组成份上的cos2之和为1，因为cos2为cor的平方，所以也认为是相关性越强，其结果与cor类似。
#cos2是coord的平方，表征特定变量在所有PC上的代表性，某个变量的所有cos2总和为1
library("corrplot")
corrplot(get_pca_var(rt4.pca2)$cos2, is.corr=FALSE, tl.cex=0.2, tl.col="blue")


fviz_cos2(rt4.pca2, choice = "var", axes = 1:2) 
# cos2在主成分上的加和，并排序

#不同按照cos2大小设定颜色梯度，也可以设置alpha梯度, axesc(1,2)代表PC1和PC2，axesc(2,3)代表PC2和PC3
fviz_pca_var(rt4.pca2,axes=c(1,2),
 col.var = "cos2",
gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE) # Avoid text overlapping

#contrib是每个变量对某一特定PC的贡献
#contrib=(var.cos2 * 100) / (total cos2 of the PC component)
#多个变量的contrib = [(Contrb1 * Eig1) + (Contrib2 * Eig2)]/(Eig1 + Eig2)
get_pca_var(rt4.pca2)$contrib
corrplot(get_pca_var(rt4.pca2)$contrib, is.corr=FALSE) 
fviz_contrib(rt4.pca2, choice = "var", axes = 1:5)

#根据contribution将变量颜色分类
fviz_pca_var(rt4.pca2,col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

#人为分组
bb<-as.factor(c(rep("Positive Contribution",32),rep("Negative Contribution",21)))  ###分组数目需要和总的基因数一致，因此这里的数字进行修改
bb<-as.factor(c(rep("Cluster1",25),rep("Cluster2",91)))  ###分组数目需要和总的基因数一致，因此这里的数字进行修改
bb<-as.factor(c(rep("Cluster3",121))) 

names(bb)<-row.names(rt4.pca2$var$contrib)
fviz_pca_var(rt4.pca2, col.var = bb, palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
             legend.title = "Cluster")


#样本坐标可视化
get_pca_ind(rt4.pca2) #coord cor cos2 contribution
rt4.pca2$ind #coord cos2 contrib dist
fviz_pca_ind(rt4.pca2)#score 可视化coord

# pointshape=635,      The length of fill variableshould be the same as the number of rows in the data.
fviz_pca_ind(rt4.pca2, geom=c("point","text"),
          addEllipses = T,
          pointshape=5,col.ind="black",pointsize="cos2",
           fill.ind = as.factor(grade),palette = "npg",
           #col.ind="cos2",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
           #col.ind="contrib",  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
         # label=wine.class,
            repel = TRUE)


fviz_pca_ind(rt4.pca2, axes = c(1, 2),label="none", #habillage只能是分类变量,一定要是因子形式
             addEllipses = TRUE, ellipse.type="norm",ellipse.level=0.9,
             # one of "confidence","t","norm","euclid","convex"
             habillage = as.factor(grade),palette = "jco",
             mean.point=F
             )

library(ggplot2)
library(ggpubr)
fviz_pca_biplot(rt4.pca2, axes = c(1,2),repel = F,
                addEllipses = T,ellipse.alpha=0.15,
                geom=c("point"),geom.var=c("arrow","text"),
                arrowsize=0.5,labelsize=2, #arrow与text大小
                pointshape=21,pointsize=1.5,fill.ind = as.factor(grade),col.ind="black", #point
               col.var=factor(c(rep("Cluster3",121),rep("Cluster2",0)))              
)%>%ggpar(xlab="PC1",ylab="PC2",title="PCA-Biplot",
          font.x=10,font.y=10,font.tickslab = 14,ggtheme=theme_bw(),xlim=c(-15,15),ylim=c(-2,3.5),
          legend.title = list(color="Variable",fill="Class"),font.legend = 12,
          )+
     ggpubr::fill_palette("jco")+ggpubr::color_palette("npg")+
  theme(axis.ticks.length= unit(-0.25, 'cm'), #设置y轴的刻度长度为负数，即向内
        axis.text.y.left  = element_text(margin = unit(c(0.5, 0.5, 0.5, 0.05), 'cm')),
        axis.text.x.bottom   = element_text(margin = unit(c(0.5, 0.5, 0.05, 0.5), 'cm'))
        )

















#设置工作路径
getwd()
setwd("G:\\OneDrive - 西湖大学\\Ma Lab\\Analysis_for_lab\\Shiniang\\Others\\scRNA")
#下载并加载所需的R包
#install.packages("Seurat")
#install.packages("patchwork")
#BiocManager::install("GEOquery")
#BiocManager::install("dplyr")
#BiocManager::install("ggplot")
#BiocManager::install("ggplot2")
rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(patchwork)
library(data.table)
library(ggplot2)
library(Seurat)
library(GEOquery)



#参考网址
##单细胞测序生信分析1-Seurat
https://www.jianshu.com/p/e7303b4db9b8
http://www.bio-info-trainee.com/6392.html
setwd()
#1.有直接的标准10X数据（喜大普奔）
pbmc.data <- Read10X(data.dir = "C:/Users/fhche/Desktop/GSE106273")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc", min.cells = 3, min.features = 200)
head(pbmc@meta.data)
saveRDS(pbmc,file="pbmc_raw.rds")  #用于后续分析

#2.多个10X数据可以用merge函数合并    例如：GSE135927，只有一个raw data能下载 下载后整理成GSM4038043、GSM4038044两个文件夹，分别含有barcodes.tsv、genes.tsv、matrix.mtx三个文件
###min.cells:基因最少在3个细胞中表达，min.features：细胞中最少有200个基因表达
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

#可以随意合并其中的几个分组
pbmcX = merge(pbmc1, pbmc2, 
             add.cell.ids = c("GSM4143678", "GSM4143679"),
             merge.data = TRUE)

pbmc = merge(pbmc1, y=c(pbmc3, pbmc5), 
             add.cell.ids = c("GSM4143679", "GSM4143682", "GSM4143684"),
             merge.data = TRUE)

pbmc = merge(pbmc2, y=c(pbmc4, pbmc6), 
             add.cell.ids = c( "GSM4143679",  "GSM4143682","GSM4143685" ),
             merge.data = TRUE)

pbmc = merge(pbmc1, y=c(pbmc2, pbmc3, pbmc4, pbmc5, pbmc6), 
             add.cell.ids = c("GSM4143678", "GSM4143679", "GSM4143681", "GSM4143682", "GSM4143684", "GSM4143685"),
             merge.data = TRUE)

head(pbmc@meta.data)
saveRDS(pbmc,file="pbmc_raw.rds")  #用于后续分析
save(pbmc,file = 'pbmc_raw.Rdata')
#3.GEO里只有矩阵数据    例如：GSE157703

library(data.table)
library("R.utils")
pca1 <- fread("GSM4773521_PCa1_gene_counts_matrix.txt.gz",
              data.table = F)
pca1[1:4,1:4]
d1=pca1[,-1]
rownames(d1)=pca1[,1]

pca2 <- fread("GSM4773522_PCa2_gene_counts_matrix.txt.gz",
              data.table = F)
pca2[1:4,1:4]
d2=pca2[,-1]
rownames(d2)=pca2[,1]

pbmc1 <- CreateSeuratObject(counts = d1,
                            min.cells = 3, 
                            min.features = 200,
                            project = "pca1")
pbmc2 <- CreateSeuratObject(counts = d2, 
                            min.cells = 3,
                            min.features = 200,
                            project = "pca2")
pbmc = merge(pbmc1, pbmc2, add.cell.ids = c("pca1", "pca2"),
             project = "PCA", merge.data = TRUE)
as.data.frame(pbmc@assays$RNA@counts[1:10, 1:2])
head(pbmc@meta.data)
saveRDS(pbmc,file="pbmc_raw.rds")  #用于后续分析
save(pbmc,file = 'pbmc_raw.Rdata')

#4.数据需要筛选后再构建Seurat矩阵   例如：GSE84465，需要筛选其中的Tumor细胞数据进一步分析
a<-read.table("GSE84465_GBM_All_data.csv.gz")
head(rownames(a))
tail(rownames(a),10)
a=a[1:(nrow(a)-5),]  #最后5行行名异常，剔除

b<-read.table("SraRunTable.txt",sep = ",",header = T)   #样本信息
table(b$patient_id)
table(b$Tissue)
table(b$Tissue, b$patient_id)

new.b <- b[,c("plate_id","Well","Tissue","patient_id")]
new.b$sample <- paste0("X",b.group$plate_id,".",b.group$Well)
head(new.b)
identical(colnames(a),new.b$sample)

#筛选肿瘤细胞
index<-which(new.b$Tissue=="Tumor")
group<-new.b[index,]   #筛选的是行
a.filt<-a[,index]        #筛选的是列
dim(a.filt)
identical(colnames(a.filt),group$sample)

#构建Seurat对象
library("Seurat")
datagroup=data.frame(Patient_ID=group$patient_id,
                     group=new.b[index,3],
                     row.names = group$sample)
pbmc <- CreateSeuratObject(counts = a.filt,
                           meta.data = datagroup,
                           min.cells = 3, 
                           min.features = 50)
head(pbmc@meta.data)
saveRDS(pbmc,file="pbmc_raw.rds")  #用于后续分析

#################最后记得
saveRDS(pbmc,file="pbmc_raw.rds")  #用于后续分析

#5.多个csv/txt矩阵合并：
# 读取多个csv
data1 <- read.csv("******", header = T,row.names= 1)
data2 <- read.csv("******", header = T,row.names= 1)
data3 <- read.csv("******", header = T,row.names= 1)

#如果是txt文件，就用data = read.table("******", header = T,row.names= 1)

# 合并成数据框
datan = data.frame(data1,data2,data3)

# 数据框转换成稀疏矩阵matrix
dataan <- as(as.matrix(datan), "dgCMatrix")

#6.数据分组比较多，想要R实现批量处理
library(data.table)
samples=list.files('GSE162025_RAW/')
samples 
pbmcList = lapply(samples,function(x){ 
  x=samples[1]
  print(x)
  y=file.path('GSE162025_RAW',x )  
  a=fread(y,data.table = F)
  a[1:4,1:4] 
  rownames(a)=a[,1]
  a=a[,-1]
  pbmc=CreateSeuratObject(a)
  return(pbmc)
})
pbmcList

pbmc=merge(x=pbmcList[[1]],
           y=pbmcList[ -1 ])
as.data.frame(pbmc@assays$RNA@counts[1:10, 1:2])
head(pbmc@meta.data)

library(stringr)
phe=str_split(rownames(pbmc@meta.data),'_',simplify = T)
head(phe)
pbmc@meta.data$patient=phe[,2]
pbmc@meta.data$location=phe[,3]
pbmc@meta.data$orig.ident=paste(phe[,3],phe[,2],sep = '_')
table(pbmc@meta.data) 



#########################################################
https://mp.weixin.qq.com/mp/appmsgalbum?__biz=MzIyMzMwNDQ2MA==&action=getalbum&album_id=1453861481023504385&scene=173&from_msgid=2247483783&from_itemidx=1&count=3&nolastread=1#wechat_redirect

Following Methods
#对构建的Seurat对象进行质控
###计算细胞中核糖体和线粒体基因比例并过滤MT  
pbmc[['percent.mt']]=PercentageFeatureSet(pbmc,pattern = "^MT-")  
scRNA.filter=subset(pbmc,percent.mt<5)
pbmc


#计算红细胞比例
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA.filter@assays$RNA)) 
HB.genes <- rownames(scRNA.filter@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA.filter[["percent.HB"]]<-PercentageFeatureSet(scRNA.filter, features=HB.genes) 
p0=VlnPlot(scRNA.filter,features = "percent.HB",group.by = "group")
ggplot2::ggsave(filename = 'p0.pdf')
#各组检测到的基因数量, 这里的group.by = "Patient_ID"需要根据实际情况进行修改
p1=VlnPlot(scRNA.filter,features = "nFeature_RNA",group.by = "group")
ggplot2::ggsave(filename = 'p1.pdf')
#各组检测到的counts数
p2=VlnPlot(scRNA.filter,features = "nCount_RNA",group.by = "group")
ggplot2::ggsave(filename = 'p2.pdf')
#基因数量和counts数的相关性
p3=FeatureScatter(scRNA.filter,feature1='nCount_RNA',feature2='nFeature_RNA',group.by="group")
ggplot2::ggsave(filename = 'p3.pdf')
#核糖体和线粒体基因
plot4 <- FeatureScatter(scRNA.filter, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggplot2::ggsave(filename = 'p4.pdf')
#数据的归一化和标准化
##首先去除样本/细胞效应
sce <- NormalizeData(scRNA.filter, normalization.method =  "LogNormalize",
                     scale.factor = 10000)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce,
                            selection.method = "vst", nfeatures = 2000)
##前10变化的基因
top10 <- head(VariableFeatures(sce), 20)

## plot variable features with and without labels
plot <- VariableFeaturePlot(sce)
p5 <- LabelPoints(plot=plot, points = top10, repel = TRUE)
ggplot2::ggsave(filename = 'p5.pdf',width =12,height = 8)

#########################################################
#首先需要降维聚类分群，代码比较简单，如下：

load(file = 'pbmc_raw.Rdata') 
library(Seurat)
# 步骤 ScaleData 的耗时取决于电脑系统配置（保守估计大于一分钟）
sce.big <- ScaleData(object = pbmc, 
 vars.to.regress = c('nCount_RNA'), 
 model.use = 'linear', 
 use.umi = FALSE)
 
sce.big <- FindVariableFeatures(object = sce.big, 
 mean.function = ExpMean, 
 dispersion.function = LogVMR, 
 x.low.cutoff = 0.0125, 
 x.high.cutoff = 4, 
 y.cutoff = 0.5)
 
 #PCA降维

length(VariableFeatures(sce.big)) 
sce.big <- RunPCA(object = sce.big, pc.genes = VariableFeatures(sce.big))
DimHeatmap(sce.big, dims = 1:12, cells = 100, balanced = TRUE)
#library("VizPCA")
#VizPCA( sce.big, pcs.use = 1:2)

#JackStrawPlot()函数提供可视化方法，用于比较每一个主成分的p-value的分布。
#虚线是均匀分布；显著的主成分富集有小p-Value基因，实线位于虚线左上方。下图表明保留10个pca主成分用于后续分析是比较合理的。
sce.big <- JackStraw(sce.big, num.replicate = 100)
sce.big <- ScoreJackStraw(sce.big, dims = 1:12)
E1 = JackStrawPlot(sce.big, dims = 1:12)
ggplot2::ggsave(filename = 'E1.pdf')

#ElbowPlot展示每个主成分对数据方差的解释情况
E2 <- ElbowPlot(sce.big)
ggplot2::ggsave(filename = 'E2.pdf')

#PCA中PC1和PC2做的因子得分图
E3 = DimPlot(sce.big, reduction = "pca",group.by  ="group")
ggplot2::ggsave(filename = 'E3.pdf')

#下面只是展现不同降维算法而已，并不要求都使用一次。
sce.big <- RunICA(sce.big )
sce.big <- RunTSNE(sce.big )
sce.big <- RunUMAP(sce.big,dims = 1:10)

DimPlot(object = sce.big, reduction = "pca") 
DimPlot(object = sce.big, reduction = "ica")
DimPlot(object = sce.big, reduction = "tsne")

# 针对PCA降维后的表达矩阵进行聚类 FindNeighbors+FindClusters 两个步骤,  注意调整resolution 参数来获得不同数目的细胞亚群,现在是分了16个cluster。
sce.big <- FindNeighbors(object = sce.big, dims = 1:20, verbose = FALSE) 
sce.big <- FindClusters(object = sce.big, resolution = 0.3,verbose = FALSE)
table(sce.big$orig.ident,sce.big@meta.data$RNA_snn_res.0.3)
table(sce.big@meta.data$group,sce.big@meta.data$RNA_snn_res.0.3)

library(gplots)
tab.1=table(sce.big@meta.data$RNA_snn_res.0.3,sce.big@meta.data$RNA_snn_res.0.3)
balloonplot(tab.1)


#tSNE降维与细胞分群
#install.packages("gplots")
##图片的输出可以通过以下代码进行修改
#ggplot2::ggsave(filename = "F1.pdf",width =12,height = 9)
set.seed(1234) ##这里为设置随机数
sce.big <- RunTSNE(sce.big, dims = 1:16, do.fast = TRUE )
F1 = DimPlot(object = sce.big, reduction = "tsne",group.by = 'RNA_snn_res.0.3',label=T)
ggplot2::ggsave(filename = 'F1.pdf',width =12,height = 9)
F2 = DimPlot(object = sce.big, reduction = "tsne",group.by = 'group',label=T)
ggplot2::ggsave(filename = 'F2.pdf',width =12,height = 9)


#Marker 在不同细胞类群中的表达情况
phe=data.frame(cell=rownames(sce.big@meta.data),
               cluster =sce.big@meta.data$seurat_clusters)
head(phe)
table(sce.big@meta.data$seurat_clusters)

for( i in unique(sce.big@meta.data$seurat_clusters) ){
  markers_df <- FindMarkers(object = sce.big, ident.1 = i, min.pct = 0.25)
  print(x = head(markers_df))
  markers_genes =  rownames(head(x = markers_df, n = 5))
  VlnPlot(object = sce.big, features =markers_genes,log =T )
  FeaturePlot(object = sce.big, features=markers_genes )
}

sce.markers <- FindAllMarkers(object = sce.big, only.pos = TRUE, min.pct = 0.25, 
thresh.use = 0.25)

###保存下R的Image，防止R studio aborted 数据丢失

##保存图片
FeaturePlot(sce.big, features = c("PLP1","CHI3L1","TF","HSPA6","APOD","DCN","LTF","MBP","STMN2","KCNQ2","LHFPL3"))
FeaturePlot(sce.big, features = c("CA9","SCNN1A","TNNC1","SLC37A4","SCARB2","CTR9","APEH","DPP8","CD93","LDHB","UEVLD","VWF","DHRS11","SRM","SUSD1","C1orf226","SLC12A1","ERMP1","PLA2G4F","DNAJC3"))
FeaturePlot(sce.big, features = c("VHL","PBRM1","IGJ","SPP1","TPSAB1","TPSB2","CLU","CXCL10","CCL17","CPA3","APOE","APOC1"))
FeaturePlot(sce.big, features = c("SYT7","PLA2G4F","GLDC","ABCA8","CCDC160","SCPEP1","CPVL"," FLRT3","ISLR","MAB21L4","SLC5A2","SLC5A3","ERMP1","CA10","CA4","CA8","CA2","PNP","HBA1","IRF6","CYP4F2","CYP4F3","CYP4A11","CYP26B1","XPNPEP2","BDKRB2","KLK1","TMPRSS4","KLK6","ST14"))
NegContribution=c("SYT7","PLA2G4F","GLDC","ABCA8","CCDC160","SCPEP1","CPVL"," FLRT3","ISLR","MAB21L4","SLC5A2","SLC5A3","ERMP1","CA10","CA4","CA8","CA2","PNP","HBA1","IRF6","CYP4F2","CYP4F3","CYP4A11","CYP26B1","XPNPEP2","BDKRB2","KLK1","TMPRSS4","KLK6","ST14")
rt2=read.table("OverlappedGenes.txt",sep="\t",header=F,check.names=F)
rt3=read.table("PosContribution.txt",sep="\t",header=F,check.names=F)
PosContribution=rt3[,1]
FeaturePlot(sce.big, features = PosContribution)

##目的基因的表达dotplot
genelist = c("VHL","PBRM1","CA9","SCNN1A","TNNC1","SLC37A4","SCARB2","CTR9","APEH","DPP8","CD93","LDHB","UEVLD","VWF","DHRS11","SRM","SUSD1","C1orf226","SLC12A1","ERMP1","PLA2G4F","DNAJC3")
p6 = DotPlot(sce.big, features = PosContribution)
ggplot2::ggsave(filename = 'p6.pdf',width =12,height = 7)
p7 = DotPlot(sce.big, features = NegContribution)
ggplot2::ggsave(filename = 'p7.pdf',width =12,height = 7)

FeaturePlot(sce.big, features = c("KRAS","NRAS","HRAS","MCRS1","TP53","MYC","FABP1"))
#每个cluster中的marker基因表达分析
## find markers for every cluster compared to all remaining cells, report only the positive ones
head(sce.markers)
##PS：看下是logFC还是log2FC
new_top10 <- sce.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
F3 = DoHeatmap(sce.big,new_top10$gene,size=3)
ggplot2::ggsave(filename = 'F3.pdf',width =14,height = 18)



####差异基因注释细胞

##以下方法三选一，建议第一种
#默认wilcox方法
diff.wilcox = FindAllMarkers(sce.big)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(all.markers, "cell_identify/diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "cell_identify/top10_diff_genes_wilcox.csv", row.names = F)
#专为单细胞设计的MAST
diff.mast = FindAllMarkers(scRNA, test.use = 'MAST')
all.markers = diff.mast %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(all.markers, "cell_identify/diff_genes_mast.csv", row.names = F)
write.csv(top10, "cell_identify/top10_diff_genes_mast.csv", row.names = F)
#bulkRNA经典方法DESeq2
diff.deseq2 = FindAllMarkers(scRNA, test.use = 'DESeq2', slot = 'counts')
all.markers = diff.deseq2 %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.csv(all.markers, "cell_identify/diff_genes_deseq2.csv", row.names = F)
write.csv(top10, "cell_identify/top10_diff_genes_deseq2.csv", row.names = F)
##top10基因绘制热图
top10_genes <- read.csv("cell_identify/top10_diff_genes_wilcox.csv")
top10_genes = CaseMatch(search = as.vector(top10_genes$gene), match = rownames(scRNA)) 
plot1 = DoHeatmap(scRNA, features = top10_genes, group.by = "seurat_clusters", group.bar = T, size = 4)
ggsave("cell_identify/top10_markers.pdf", plot=plot1, width=8, height=6) 
ggsave("cell_identify/top10_markers.png", plot=plot1, width=8, height=6)


#挑选部分基因
select_genes <- c("KRAS","NRAS","HRAS","MCRS1","TP53","MYC","WNT")
#vlnplot展示
d1 <- VlnPlot(scRNA, features = select_genes, pt.size=0, group.by="celltype", ncol=2)
ggsave("cell_identify/selectgenes_VlnPlot.png", d1, width=6 ,height=8)
#featureplot展示
d2 <- FeaturePlot(scRNA, features = select_genes, reduction = "tsne", label=T, ncol=2)
ggsave("cell_identify/selectgenes_FeaturePlot.png", d2, width=8 ,height=12)
d3=d1|d2
ggsave("cell_identify/selectgenes.png", p3, width=10 ,height=8)



##SingleR鉴定细胞类型
library(SingleR)
library(celldex)
refdata <- HumanPrimaryCellAtlasData()
testdata <- GetAssayData(sce.big, slot="data")
clusters <- sce.big@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                     method = "cluster", clusters = clusters, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"celltype_singleR.csv",row.names = F)
sce.big@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
sce.big@meta.data[which(sce.big@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

library(cowplot) 

G1 <- DimPlot(sce.big, reduction = "tsne", group.by = "celltype", pt.size=0.5)
ggplot2::ggsave(filename = 'G1.pdf',width =12,height = 9)

save(sce.big,file = 'sce.big.Rdata')


##细胞类型的再聚类
library(Seurat)
library(monocle)
library(tidyverse)
library(patchwork)
dir.create("subcluster")
#提取细胞子集
Cells.sub <- subset(sce.big@meta.data, celltype=="Endothelial_cells:lymphatic")
scRNAsub <- subset(sce.big, cells=row.names(Cells.sub))
#重新降维聚类：因为再聚类的细胞之间差异比较小，所以聚类函数FindClusters()控制分辨率的参数建议调高到resolution = 0.9。
##PCA降维
scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <-  rownames(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
ElbowPlot(scRNAsub, ndims=20, reduction="pca")
pc.num=1:10
##细胞聚类
scRNAsub <- FindNeighbors(scRNAsub, dims = pc.num) 
scRNAsub <- FindClusters(scRNAsub, resolution = 0.9)
table(scRNAsub@meta.data$seurat_clusters)
metadata <- scRNAsub@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'subcluster/cell_cluster.csv',row.names = F)
##非线性降维
#tSNE
scRNAsub = RunTSNE(scRNAsub, dims = pc.num)
embed_tsne <- Embeddings(scRNAsub, 'tsne')
write.csv(embed_tsne,'subcluster/embed_tsne.csv')
plot1 = DimPlot(scRNAsub, reduction = "tsne") 
ggsave("subcluster/tSNE.pdf", plot = plot1, width = 8, height = 7)
ggsave("subcluster/tSNE.png", plot = plot1, width = 8, height = 7)
#UMAP
scRNAsub <- RunUMAP(scRNAsub, dims = pc.num)
embed_umap <- Embeddings(scRNAsub, 'umap')
write.csv(embed_umap,'subcluster/embed_umap.csv') 
plot2 = DimPlot(scRNAsub, reduction = "umap") 
ggsave("subcluster/UMAP.pdf", plot = plot2, width = 8, height = 7)
ggsave("subcluster/UMAP.png", plot = plot2, width = 8, height = 7)
#合并tSNE与UMAP
plotc <- plot1+plot2+ plot_layout(guides = 'collect')
ggsave("subcluster/tSNE_UMAP.pdf", plot = plotc, width = 10, height = 5)
ggsave("subcluster/tSNE_UMAP.png", plot = plotc, width = 10, height = 5)


#Cluster差异分析
diff.wilcox = FindAllMarkers(scRNAsub)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
##要看是avg_log2FC还是avg_logFC
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "subcluster/diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "subcluster/top10_diff_genes_wilcox.csv", row.names = F)

#SingleR细胞鉴定
#Subcluster的细胞同样可以使用SingleR鉴定细胞类型。使用的时候注意调整参考数据库和分类标签，以便鉴定结果更有针对性。
#上节使用SingleR时使用的参考数据库是人类主要细胞图谱（HumanPrimaryCellAtlasData），采用分类标签是主分类标签（label.main）；
#这次建议使用人类免疫细胞数据（MonacoImmuneData），分类标签采用精细分类标签（label.fine）。

##精细细胞类型鉴定,如果是免疫细胞的话，就要选MonacoImmuneData
library(SingleR)
#refdata <- MonacoImmuneData()
refdata <- HumanPrimaryCellAtlasData()
testdata <- GetAssayData(scRNAsub, slot="data")
clusters <- scRNAsub@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine, 
                     method = "cluster", clusters = clusters, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)
write.csv(celltype,"subcluster/celltype_singleR.csv",row.names = F)
scRNAsub@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
scRNAsub@meta.data[which(scRNAsub@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
pp1 = DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='tsne')
pp2 = DimPlot(scRNAsub, group.by="celltype", label=T, label.size=5, reduction='umap')
pp3 = plotc <- pp1+pp2+ plot_layout(guides = 'collect')
ggsave("subcluster/tSNE_celltype.pdf", pp1, width=12 ,height=9)
ggsave("subcluster/UMAP_celltype.pdf", pp2, width=12 ,height=9)
ggsave("subcluster/celltype.pdf", pp3, width=12 ,height=9)
ggsave("subcluster/celltype.png", pp3, width=12 ,height=9)

save(scRNAsub,file = 'scRNAsub.Rdata')

#伪时间分析
library(monocle)
dir.create("pseudotime")
#scRNAsub是上一节保存的T细胞子集seurat对象
data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
#mycds <- detectGenes(mycds, min_expr = 2)  #很多教程不用

###选择代表性基因
##使用clusters差异表达基因
diff.genes <- read.csv('subcluster/diff_genes_wilcox.csv')
diff.genes <- subset(diff.genes,p_val_adj<0.01)$gene
mycds <- setOrderingFilter(mycds, diff.genes)
p1 <- plot_ordering_genes(mycds)
##使用seurat选择的高变基因
var.genes <- VariableFeatures(scRNAsub)
mycds <- setOrderingFilter(mycds, var.genes)
p2 <- plot_ordering_genes(mycds)
##使用monocle选择的高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds)
##结果对比
p1|p2|p3

###降维及细胞排序 使用disp.genes开展后续分析

#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
#State轨迹分布图
plot1 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("pseudotime/State.pdf", plot = plot1, width = 6, height = 5)
ggsave("pseudotime/State.png", plot = plot1, width = 6, height = 5)
##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
ggsave("pseudotime/Cluster.pdf", plot = plot2, width = 6, height = 5)
ggsave("pseudotime/Cluster.png", plot = plot2, width = 6, height = 5)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("pseudotime/Pseudotime.pdf", plot = plot3, width = 6, height = 5)
ggsave("pseudotime/Pseudotime.png", plot = plot3, width = 6, height = 5)
##合并作图
plotc <- plot1|plot2|plot3
ggsave("pseudotime/Combination.pdf", plot = plotc, width = 10, height = 3.5)
ggsave("pseudotime/Combination.png", plot = plotc, width = 10, height = 3.5)
##保存结果
write.csv(pData(mycds), "pseudotime/pseudotime.csv")

#轨迹图分面显示

p1 <- plot_cell_trajectory(mycds, color_by = "State") + facet_wrap(~State, nrow = 1)
p2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters") + facet_wrap(~seurat_clusters, nrow = 1)
plotc <- p1/p2
ggsave("pseudotime/trajectory_facet.png", plot = plotc, width = 6, height = 5)

#Monocle基因可视化
s.genes <- c("ITGB1","CCR7","KLRB1","GNLY")
p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
plotc <- p1|p2|p3
ggsave("pseudotime/genes_visual.png", plot = plotc, width = 8, height = 4.5)

#拟时相关基因聚类热图

#cluster差异基因
diff.genes <- read.csv('subcluster/diff_genes_wilcox.csv')
sig_diff.genes <- subset(diff.genes,p_val_adj<0.0001&abs(avg_logFC)>0.75)$gene
sig_diff.genes <- unique(as.character(sig_diff.genes))
diff_test <- differentialGeneTest(mycds[sig_diff.genes,], cores = 1, 
                              fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 0.01))
p1 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=3,
                             show_rownames=T, return_heatmap=T)
ggsave("pseudotime/pseudotime_heatmap1.png", plot = p1, width = 5, height = 8)
#高变基因
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
diff_test <- differentialGeneTest(mycds[disp.genes,], cores = 4, 
                              fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test, qval < 1e-04))
p2 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=5,
                             show_rownames=T, return_heatmap=T)
ggsave("pseudotime/pseudotime_heatmap2.png", plot = p2, width = 5, height = 10)

#BEAM分析
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
mycds_sub <- mycds[disp.genes,]
plot_cell_trajectory(mycds_sub, color_by = "State")
beam_res <- BEAM(mycds_sub, branch_point = 1, cores = 8)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]
mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 1e-4)),]
plot_genes_branched_heatmap(mycds_sub_beam,  branch_point = 1, num_clusters = 3, show_rownames = T)



###基因差异表达分析
library(Seurat)
library(tidyverse)
library(patchwork)
library(monocle)
library(clusterProfiler)
library(org.Hs.eg.db)
dir.create("enrich")

#比较cluster0和cluster1的差异表达基因
dge.cluster <- FindMarkers(sce.big,ident.1 = 0,ident.2 = 1)
sig_dge.cluster <- subset(dge.cluster, p_val_adj<0.01&abs(avg_log2FC)>1)
#比较B_cell和T_cells的差异表达基因
dge.celltype <- FindMarkers(sce.big, ident.1 = 'B_cell', ident.2 = 'T_cells', group.by = 'celltype')
sig_dge.celltype <- subset(dge.celltype, p_val_adj<0.01&abs(avg_log2FC)>1)
#比较拟时State1和State3的差异表达基因
p_data <- subset(pData(mycds),select='State')
scRNAsub <- subset(sce.big, cells=row.names(p_data))
scRNAsub <- AddMetaData(scRNAsub,p_data,col.name = 'State')
dge.State <- FindMarkers(scRNAsub, ident.1 = 1, ident.2 = 3, group.by = 'State')
sig_dge.State <- subset(dge.State, p_val_adj<0.01&abs(avg_logFC)>1)



#差异基因GO富集分析
ego_ALL <- enrichGO(gene          = row.names(sig_dge.cluster),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)
write.csv(ego_all,'enrich/enrichGO.csv')           
ego_CC <- enrichGO(gene          = row.names(sig_dge.cluster),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_MF <- enrichGO(gene          = row.names(sig_dge.cluster),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- enrichGO(gene          = row.names(sig_dge.cluster),
                   #universe     = row.names(dge.celltype),
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)           
ego_CC@result$Description <- substring(ego_CC@result$Description,1,70)
ego_MF@result$Description <- substring(ego_MF@result$Description,1,70)
ego_BP@result$Description <- substring(ego_BP@result$Description,1,70)
p_BP <- barplot(ego_BP,showCategory = 10) + ggtitle("barplot for Biological process")
p_CC <- barplot(ego_CC,showCategory = 10) + ggtitle("barplot for Cellular component")
p_MF <- barplot(ego_MF,showCategory = 10) + ggtitle("barplot for Molecular function")
plotc <- p_BP/p_CC/p_MF
ggsave('enrich/enrichGO.png', plotc, width = 12,height = 10)


#差异基因kegg富集分析
genelist <- bitr(row.names(sig_dge.cluster), fromType="SYMBOL",
                           toType="ENTREZID", OrgDb='org.Hs.eg.db')
genelist <- pull(genelist,ENTREZID)               
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa')
p1 <- barplot(ekegg, showCategory=20)
p2 <- dotplot(ekegg, showCategory=20)
plotc = p1/p2
ggsave("enrich/enrichKEGG.png", plot = plotc, width = 12, height = 10)


















##细胞间通讯CellChat
library(CellChat)
library(patchwork)

cellchat <- createCellChat(object = sce.big, meta = meta, group.by = "celltype")

##细胞间通讯iTALK
https://www.jianshu.com/p/dfe35a2a02d4

library(iTALK)
library(Seurat)
library(Matrix)
library(dplyr)

# iTALK 要求的矩阵: 行为细胞，列为基因
iTalk_data <- as.data.frame(t(sce.big@assays$RNA@counts))
# iTALK 要求包含cell_type列，我的细胞分群存储在seurat_cluster
iTalk_data$cell_type <- sce.big@meta.data$celltype
# iTALK 要求包含compare_group列（多样本），表示每个细胞的生物学分组/样本，我的细胞分组存放在Group
iTalk_data$compare_group <- sce.big@meta.data$Group

unique(iTalk_data$cell_type)
# "cd56_nk" "cd14_monocytes" "b_cells" "cytotoxic_t" "regulatory_t" "memory_t" "naive_t"
unique(iTalk_data$compare_group)
# "group1" "group2" "group3"

#通过所有细胞的高表达基因分析其中包含的配体-受体
my10colors <- my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282')

highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
# 通讯类型
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)

iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}

iTalk_res <- iTalk_res[order(iTalk_res$cell_from_mean_exprs*iTalk_res$cell_to_mean_exprs,decreasing=T),][1:20,]        
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)

LRPlot(iTalk_res[1:20,],datatype='mean count',cell_col=cell_col,link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:20],link.arr.width=iTalk_res$cell_to_mean_exprs[1:20])
























#细胞分群的自动化注释
#Warning message: package(s) not installed when version(s) same as current;  use `force = TRUE` to re-install: 'SingleR' 
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("tidyverse")
#BiocManager::install("viridis")
#BiocManager::install("mapproj")
#BiocManager::install("pheatmap")
library(SingleR)
library(celldex)
library(tidyverse)
library(viridis)
library(mapproj)
library(ggplot2)
library(pheatmap)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
save(hpca.se,file = 'hpca.se.Rdata')

##注意参考集里面是logcounts矩阵，后面对于单细胞数据集也要做类似的处理。 导入UMI count矩阵
#BiocManager::install("scater")
#BiocManager::install("SummarizedExperiment")
library(scater)
library(SummarizedExperiment)

##当合并了多个数据后，提取test.count容易引起内存占用爆炸，所以先查看有多少行数据，分批处理10000一批
nrow(sce.big[["RNA"]]@counts)
test.count1=as.data.frame(sce.big[["RNA"]]@counts[(1:5000),])
test.count2=as.data.frame(sce.big[["RNA"]]@counts[(5001:10000),])
test.count3=as.data.frame(sce.big[["RNA"]]@counts[(10001:15000),])
test.count4=as.data.frame(sce.big[["RNA"]]@counts[(15000:19051),])
test.count=merge(test.count1,test.count2)

test.count=as.data.frame(sce.big[["RNA"]]@counts)
load(file="hpca.se.RData")
common_hpca <- intersect(rownames(test.count), rownames(hpca.se))
hpca.se <- hpca.se[common_hpca,]
test.count_forhpca <- test.count[common_hpca,]
test.count_forhpca.se <- SummarizedExperiment(assays=list(counts=test.count_forhpca))
test.count_forhpca.se <- logNormCounts(test.count_forhpca.se)

##hpca的main大类注释，fine小类注释不是很准
pred.main.hpca <- SingleR(test = test.count_forhpca.se, ref = hpca.se, labels = hpca.se$label.main)
##保存热图
plotScoreHeatmap(results = pred.hesc)

result_main_hpca <- as.data.frame(pred.main.hpca$labels)
result_main_hpca$CB <- rownames(pred.main.hpca)
colnames(result_main_hpca) <- c('HPCA_Main', 'CB')
write.table(result_main_hpca, file = "HPCA_Main.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE) #保存下来，方便以后调用
head(result_main_hpca)

##把这个结果整合到test.seu@meta.data中，然后画tsne/umap展示一下
sce@meta.data$CB=rownames(sce.big@meta.data)
sce@meta.data=merge(sce.big@meta.data,result_main_hpca,by="CB")
rownames(sce.big@meta.data)=sce@meta.data$CB

##注释图
#install.packages("cowplot") 
library(cowplot) 

G1 <- DimPlot(sce.big, reduction = "tsne", group.by = "HPCA_Main", pt.size=0.5)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
G2 <- DimPlot(sce.big, reduction = "tsne", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
fig_tsne <- plot_grid(G2, G1, labels = c('ident','HPCA_Main'),rel_widths = c(2,3))
ggsave(filename = "tsne_final.pdf", plot = fig_tsne, device = 'pdf', width = 36, height = 12, units = 'cm')


#细胞间通讯的相关分析
https://mp.weixin.qq.com/mp/appmsgalbum?__biz=MzI1Njk4ODE0MQ==&action=getalbum&album_id=1644453005439614977&scene=173&from_msgid=2247498095&from_itemidx=1&count=3&nolastread=1#wechat_redirect


#CellChat三部曲1：使用CellChat对单个数据集进行细胞间通讯分析
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)


##CellChat具有不同细胞类型成分的多个数据集的细胞通讯比较分析
https://mp.weixin.qq.com/s/bEMky9vXRjTlc1Zb49340w
https://mp.weixin.qq.com/s?__biz=MzI1Njk4ODE0MQ==&mid=2247487845&idx=1&sn=7346517da8ce6f5b95fe7fb640a87403&chksm=ea1f17e7dd689ef1a0e176f38aac5cc8e70a67d2c0aeb11bd367e9c8b4dac620dbcc9a99cac2&scene=21#wechat_redirect


library(CellChat)
library(ggplot2)                  
library(patchwork)
library(igraph)















