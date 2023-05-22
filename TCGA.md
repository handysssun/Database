# TCGA数据库（Access TCGA Data）
## 1数据检索功能（部分需要账号）
1.1原始测序数据——fasta和fastq格式文件
1.2比对好的bam格式文件
1.3经过处理和标准化的文件（完全开放）
  
## 2 数据下载R包（支持断点续传）
```
library(TCGAbiolinks)
library(SummarizedExperiment) #加载R 包

work_dir <- "D:/colonadenocarcinoma" # 选择工程地点（也就是数据下载的位置）
project <- "TCGA-COAD" # 选择工程
data_category <- "TranomeProfiling" # 类似于在网站上直接进行选择，选择转录数据
data_type <- "Gene ExpressionQuantification" # 选择基因表达谱数据
workflow_type <- "HTSeq -Counts" # 选择counts 数据
legacy <- FALSE # 使用hg38

DataDirectory <-paste0(work_dir,"/GDC/",gsub("-","_",project))
FileNameData <- paste0(DataDirectory,"_","RNAseq_HTSeq_Counts",".rda")
query <- GDCquery(project = project, data.category = data_category, data.type = data_type, workflow.type = workflow_type, legacy = legacy) # 查询下载的数据情况

# 显示下载数据的总样本量
samplesDown <-getResults(query,cols=c("cases"))
cat("Total sample to download:",length(samplesDown))

# 显示下载数据的肿瘤样本量
dataSmTP <-TCGAquery_SampleTypes(barcode = samplesDown,typesample ="TP")
cat("Total TP samples to down:",length(dataSmTP))

# 显示下载数据的正常样本量
dataSmNT <-TCGAquery_SampleTypes(barcode = samplesDown,typesample ="NT")
cat("Total NT samples to down:",length(dataSmNT))

# 下载并整合数据
GDCdownload(query = query,directory = DataDirectory,files.per.chunk=6,method='client')
data <- GDCprepare(query = query,save = TRUE,directory = DataDirectory,save.filename =FileNameData)
data_expr <- assay(data)
dim(data_expr)
gene_expr_file <- paste0(DataDirectory, "_", "Gene_HTSeq_Counts", ".txt")
write.csv(data_expr,file ='raw_mRNAdata.csv')

```
另一种数据下载代码
```
#下载数据
rm(list=ls())
library(TCGAbiolinks)
library(stringr)
library(SummarizedExperiment)
projs=getGDCprojects()$project_id %>%
      str_subset("TCGA-LUAD")
projs
### 2.下载并整理表达矩阵
proj = projs
f1 = paste0(proj,"expf.Rdata")
if(!file.exists(f1)){
  query = GDCquery(project = proj, 
                   data.category = "Transcriptome Profiling",
                   data.type = "Gene Expression Quantification", 
                   workflow.type = "STAR - Counts",
  )
  
  GDCdownload(query)
  dat = GDCprepare(query)
  exp = assay(dat)
  #tpm = assay(dat,4)
  save(exp,file = f1)
}
load(f1)
### 3.下载并整理临床信息
f2 = paste0(proj,"clf.Rdata")
if(!file.exists(f2)){
  query = GDCquery(project = proj, 
                   data.category = "Clinical",
                   data.type = "Clinical Supplement",
                   file.type = "xml"
  )
  GDCdownload(query)
  dat = GDCprepare_clinic(query,clinical.info = "patient")
  k = apply(dat, 2, function(x){!all(is.na(x))});table(k) #选出空白列
  clinical = dat[,k] #去除空白列
  save(clinical,file = f2)
}
load(f2)
table(clinical$stage_event_pathologic_stage)
sum(table(clinical$stage_event_pathologic_stage))
### 4.表达矩阵行名ID转换
library(tinyarray)
exp = trans_exp_new(exp) #基因注释
exp[1:4,1:4]
### 5.基因过滤
exp = exp[apply(exp, 1, function(x) sum(x > 0) > 0.5*ncol(exp)), ]
nrow(exp)
### 6.分组信息获取
Group = make_tcga_group(exp) #根据列名区分肿瘤组织&正常组织
table(Group)
### 7.保存数据
save(exp,Group,proj,clinical,file = paste0(proj,".Rdata"))
```


## 3整理矩阵















