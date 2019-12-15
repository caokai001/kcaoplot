#! bin/bash
#######################
# 目的： 提取每个TF非center区域突变富集情况，保存到txt。(按行存储)，python 依次读取
#######################   

# 加载包
library("tidyverse")
library("kcaoplot")
library("GenomicRanges")
library("parallel")

# 参数解析
singel_TF_mutation_dir <- "/public/home/kcao/Desktop/2018_NC_GS/Jasper_tf_somatic/"
# paired_TF_mutation_dir <- "/public/home/kcao/Desktop/2018_NC_GS/Jasper_paired_tf_somatic_norm/"
paired_TF_mutation_dir <- "/public/home/kcao/Desktop/2018_NC_GS/Jasper_paired_tf_somatic_norm_test/"
FIMO_TF_scan_dir <- "/public/home/kcao/genome_human/hg19_fimo/motif_dir/"

# 输出目录
tf_all_background <- "/public/home/kcao/Desktop/2018_NC_GS/tf_all_background.txt"
tf_all_center  <- "/public/home/kcao/Desktop/2018_NC_GS/tf_all_center.txt"



# 1.选择所有的TF 
tf_all <- list.files(FIMO_TF_scan_dir,pattern='*.meme')  # 747

# 2.开放区域合并
  dnase.peak <- bed2GRanges("/public/home/kcao/Desktop/2018_NC_GS/TESTDATA/E094-DNase.all.peaks.v2.bed")
  dnase.fdr <- bed2GRanges("/public/home/kcao/Desktop/2018_NC_GS/TESTDATA/E094-DNase.fdr0.01.peaks.v2.bed")
  dnase.mac <- bed2GRanges("/public/home/kcao/Desktop/2018_NC_GS/TESTDATA/E094-DNase.macs2.narrowPeak")
  dnase <- c(dnase.peak, dnase.fdr, dnase.mac)
  dnase <- GenomicRanges::reduce(dnase)
  # 3.2 读入snp信息
  somatic_gastric <- bed2GRanges("/public/home/kcao/Desktop/2018_NC_GS/TESTDATA/41467_2018_3828_MOESM6_ESM.txt", header = TRUE)
  names(mcols(somatic_gastric))[1:2] <- c("Ref", "Alt")
  somatic_gastric$id = c(1:length(somatic_gastric))

# 3.计算每个TF background 突变情况
## 命令行参数
args<-commandArgs(T)
i <-as.integer(args[1])
print(i)
# i=1
target_tf <- tf_all[i]




path_target_tf_fimo_dir <- paste0(FIMO_TF_scan_dir, str_replace(target_tf, ".meme", ""))
name_target.fimo <- list.files(path = path_target_tf_fimo_dir, pattern = "*.tsv")
# 转换格式GRanges ,cofactor.motif 换成全局变量
target.motif <- Fimo2GRanges(paste0(path_target_tf_fimo_dir, "/", name_target.fimo),"GRanges" ) 


# 4.过滤FIMO 
motif.ovl <- filter_fimo(dnase,target.motif)
# 拓展到1kb
site <- filter_fimo_to_1000bp_gr(motif.ovl)
# 与突变位置取交集
mut.table<-mut_pos_tb(site,somatic_gastric)

# 5.提取center 位置20bp 外作为background 区域
center_expend_each=20

mut.table.pos=mut.table[which(mut.table$dir=="+"),]
mut.table.neg=mut.table[which(mut.table$dir=="-"),]

## statistic postive mut number
  mut.table.pos$pos=mut.table.pos$start-mut.table.pos$motif.start
  mut.table.pos.df=aggregate(seqnames~pos,mut.table.pos,length)
  colnames(mut.table.pos.df)[2]="count"
  df=data.frame(pos=seq(-1000,1000,1),x=0)
  mut.table.pos.df=merge(df,mut.table.pos.df,all.x=TRUE,by="pos")
  mut.table.pos.df$count=ifelse(is.na(mut.table.pos.df$count),0,mut.table.pos.df$count)

######
  ## 判断 mut.table.neg 是否为零，如MYCN
  #####
  if(nrow(mut.table.neg)!=0){
    ## statistic negtive mut number
    mut.table.neg$pos=mut.table.neg$motif.start-mut.table.neg$start
    mut.table.neg.df=aggregate(seqnames~pos,mut.table.neg,length)
    colnames(mut.table.neg.df)[2]="count"
    df=data.frame(pos=seq(-1000,1000,1),x=0)
    mut.table.neg.df=merge(df,mut.table.neg.df,all.x=TRUE,by="pos")
    mut.table.neg.df$count=ifelse(is.na(mut.table.neg.df$count),0,mut.table.neg.df$count)

    # merge postive and negtive mut number
    mut.table.z=merge(mut.table.pos.df,mut.table.neg.df,by="pos")
    mut.table.z$mut=mut.table.z$count.x+mut.table.z$count.y
    # if we have many sample ；divided by the n sample
    mut.table.z$norm.mut=(mut.table.z$count.x+mut.table.z$count.y)
  }else{
    mut.table.z= mut.table.pos.df
    mut.table.z$mut = mut.table.z$count
    mut.table.z$norm.mut = mut.table.z$count
  }

#######################
# 5.2 提取 
require(tidyverse)
motif_len <- unique(nchar( head(mut.table$motif)))  # motif 长度
center_region_left <- 1001-ceiling(motif_len/2+center_expend_each)
center_region_right <- 1001+ceiling(motif_len/2+center_expend_each)

# calculate background
background <- mut.table.z$norm.mut[-c(center_region_left:center_region_right)]
# calculate center enrichment score
center_vector <- mut.table.z$norm.mut[c(center_region_left:center_region_right)]

# 6.将每个TF的一行数据保存到文件中
cat(paste0(c(target_tf ,background ),collapse = "\t"),file=tf_all_background,append = T,sep = "\n")

cat(paste0(c(target_tf ,center_vector ),collapse = "\t"),file=tf_all_center,append = T,sep = "\n")   
