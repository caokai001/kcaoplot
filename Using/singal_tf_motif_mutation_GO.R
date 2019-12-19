#!bin/bash

##################
# 对TF结合位点发生突变位置，附近基因进行GO富集
##################

# 加载包
library(kcaoplot)
library(GenomicRanges)
library(tidyverse)

# 接收参数
args=commandArgs(T)
TF_name=args[1]  # HIF1A_MA1106.1.meme

# 代码
# 1.读取FIMO.tsv
motif_dir="/public/home/kcao/genome_human/hg38_fimo/motif_dir"
target_TF_FIMO=paste0(motif_dir,"/",str_replace(TF_name,".meme",""),"/hg38_",TF_name,"_FIMO.tsv")
target.motif.gr=Fimo2GRanges(target_TF_FIMO,"GRanges")   # HIF1A.Granges

# 2.读取LNCAP.snp
LNCaP_SNP="/public/home/kcao/Desktop/TCGA_prostate/DATA_SNP/LNCaP.SNP.bed"
somatic_gastric<-bed2GRanges(LNCaP_SNP,header=TRUE)
names(mcols(somatic_gastric))[1:2]<-c("Ref","Alt")
somatic_gastric$id=c(1:length(somatic_gastric))

# 3.开放数据Dnase
DNase_macs2_narrowPeak="/public/home/kcao/Desktop/TCGA_prostate/DATA_SNP/PREC_cistrome_DNASE_rep1_45056.bed"
dnase<-bed2GRanges(DNase_macs2_narrowPeak)

# 4.过滤FIMO
convert2GR<-function(motif.ovl){
  motif_center_site<-GRanges(seqnames=motif.ovl$seqnames,
          IRanges(start=motif.ovl$start,end=motif.ovl$end),
          pval=motif.ovl$pval,
          dir=motif.ovl$dir,
          motif=motif.ovl$motif)
  return(motif_center_site)
}

# 
motif.ovl<-filter_fimo(dnase,target.motif.gr)
#
site=convert2GR(motif.ovl)

# 5.与突变位置取交集
mut.table<-mut_pos_tb(site,somatic_gastric)  # 270
write_delim(mut.table,"HINF1a_TF_binding_mutation.txt")
# 6.提取位置，注释基因  https://www.jianshu.com/p/4e2c8d7fce83
# 使用的snp 数据为hg38,下载注释信息


# 读取区间

## peak 读入
HINF1a <- read_delim("HINF1a_TF_binding_mutation.txt",delim = " ",col_names = TRUE)
convert2GR<-function(motif.ovl){
  motif_center_site<-GRanges(seqnames=motif.ovl$seqnames,
                             IRanges(start=motif.ovl$start,end=motif.ovl$end),
                             pval=motif.ovl$pval,
                             dir=motif.ovl$dir,
                             motif=motif.ovl$motif)
  return(motif_center_site)
}
region_GR = convert2GR(HINF1a)


# https://www.jianshu.com/p/4e2c8d7fce83
require(ChIPseeker)
library(GenomicFeatures)
library(clusterProfiler)
library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg38.knownGene

peakAnno <- annotatePeak(region_GR, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
View(as.data.frame(peakAnno))

gene.df=as.data.frame(peakAnno)

# 保存
write.table(drop_na(data.frame(gene.df$ENSEMBL)),file = "test",row.names = F)
tbl_df(as.data.frame(peakAnno)) %>%filter(annotation!="Distal Intergenic") %>% dplyr::select(ENSEMBL) %>%drop_na() %>% write.table("test_EMSEMBL.txt",row.names = F,col.names = F)


# GO
ego_BP <- enrichGO(gene = gene.df$ENSEMBL, 
                keyType = "ENSEMBL", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 1, 
                readable = TRUE)

head(ego_BP,2)
barplot(ego_BP,showCategory=20,drop=T)











/public/home/kcao/genome_human/hg38_fimo/motif_dir/HIF1A_MA1106.1/hg38_HIF1A_MA1106.1.meme_FIMO.tsv
