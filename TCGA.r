#TCGA数据库（Access TCGA Data）
#1 数据检索功能（部分需要账号）
#1.1原始测序数据——fasta和fastq格式文件
#1.2比对好的bam格式文件
#1.3经过处理和标准化的文件（完全开放）
  
#2 数据下载R包（支持断点续传）
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
query <- GDCquery(project = project,
data.category =data_category,
data.type = data_type,
workflow.type =workflow_type,
legacy = legacy) # 查询下载的数据情况

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

