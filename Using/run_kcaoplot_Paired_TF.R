#! bin/R
## 加载包
library("tidyverse")
library("kcaoplot")
library("GenomicRanges")
library("parallel")

## 参数解析
singel_TF_mutation_dir <- "/public/home/kcao/Desktop/2018_NC_GS/Jasper_tf_somatic/"
paired_TF_mutation_dir <- "/public/home/kcao/Desktop/2018_NC_GS/paired_Jasper_tf_somatic/"
FIMO_TF_scan_dir <- "/public/home/kcao/genome_human/hg19_fimo/motif_dir/"

## 1.读取mutation_score,选取mutation score阈值
# 返回目录 2018_NC_GS
setwd(paste0(singel_TF_mutation_dir, "/.."))
data <- read_delim("tmp_merge_score.txt", delim = "\t", col_names = FALSE)
colnames(data) <- c("TF_name", "mutation_score")
cutoff <- as.vector(quantile(data$mutation_score, 0.8)) # 1.7
# TF filter vector
filter_tf <- data[data$mutation_score > cutoff,][["TF_name"]]
# 测试3个tf 
# filter_tf <- head(filter_tf, 3)


## 2.进行paired_TF_mutation_score calculate
# 排列组合
Permute <- function(TF_vector) {
  TF_combind_matrix <- combn(TF_vector, 2)
  TF_Permute_matrix <- cbind(TF_combind_matrix, apply(TF_combind_matrix, 2, rev))
  return(TF_Permute_matrix)
}
TF_Permute_matrix <- Permute(filter_tf)

# r$> Permute(filter_tf)                                                                          
#      [,1]                   [,2]                   [,3]                  
# [1,] "HOXC12_MA0906.1.meme" "HOXC12_MA0906.1.meme" "Arid5a_MA0602.1.meme"
# [2,] "Arid5a_MA0602.1.meme" "E2F3_MA0469.3.meme"   "E2F3_MA0469.3.meme"  
#      [,4]                   [,5]                   [,6]                  
# [1,] "Arid5a_MA0602.1.meme" "E2F3_MA0469.3.meme"   "E2F3_MA0469.3.meme"  
# [2,] "HOXC12_MA0906.1.meme" "HOXC12_MA0906.1.meme" "Arid5a_MA0602.1.meme"



## 3.循环进行计算
# 改写成函数，进行mclapply计算
#RunPairsTF <- function(i) {

 for (i in seq(ncol(TF_Permute_matrix))) {
  target_tf <- TF_Permute_matrix[1, i]
  co_factor <- TF_Permute_matrix[2, i]
  work_dir=paste0(paired_TF_mutation_dir, target_tf,"/",co_factor)
  if (!dir.exists(work_dir)) {
    dir.create(work_dir, recursive = TRUE)
    # system(paste0("mkdir -p ",work_dir))
  }

  setwd(work_dir)
  # print(c(target_tf,co_factor))

  # 4.1 读入数据target_tf,co_factor
  path_target_tf_fimo_dir <- paste0(FIMO_TF_scan_dir, str_replace(target_tf, ".meme", ""))
  path_co_factor_fimo_dir <- paste0(FIMO_TF_scan_dir, str_replace(co_factor, ".meme", ""))

  name_target.fimo <- list.files(path = path_target_tf_fimo_dir, pattern = "*.tsv")
  name_co_factor.fimo <- list.files(path = path_co_factor_fimo_dir, pattern = "*.tsv")

  # 转换格式GRanges
  target.motif <- Fimo2GRanges(paste0(path_target_tf_fimo_dir, "/", name_target.fimo),"GRanges" ) 
  cofactor.motif <- Fimo2GRanges(paste0( path_co_factor_fimo_dir,"/",name_co_factor.fimo), "GRanges")

  # warning(print(c(paste0(path_target_tf_fimo_dir, "/", name_target.fimo),"\n",paste0( path_co_factor_fimo_dir,"/",name_co_factor.fimo))))
  # warning(print(c(length(target.motif), "\t", length(cofactor.motif)))) # ok

  # 4.2 开放区域合并
  dnase.peak <- bed2GRanges("/public/home/kcao/Desktop/2018_NC_GS/TESTDATA/E094-DNase.all.peaks.v2.bed")
  dnase.fdr <- bed2GRanges("/public/home/kcao/Desktop/2018_NC_GS/TESTDATA/E094-DNase.fdr0.01.peaks.v2.bed")
  dnase.mac <- bed2GRanges("/public/home/kcao/Desktop/2018_NC_GS/TESTDATA/E094-DNase.macs2.narrowPeak")
  dnase <- c(dnase.peak, dnase.fdr, dnase.mac)
  dnase <- GenomicRanges::reduce(dnase)
  # 4.3 读入snp信息
  somatic_gastric <- bed2GRanges("/public/home/kcao/Desktop/2018_NC_GS/TESTDATA/41467_2018_3828_MOESM6_ESM.txt", header = TRUE)
  names(mcols(somatic_gastric))[1:2] <- c("Ref", "Alt")
  somatic_gastric$id = c(1:length(somatic_gastric))


  # 4.4 过滤target_fimo_motif
  motif.ovl <- filter_fimo(dnase, target.motif, cofactor_GRanges = "cofactor.motif" , extend_size = 1000)
  # warning(length(get("target.motif")))
  # warning(length(get("cofactor.motif")))
  # warning(nrow(get("motif.ovl")))
  #  拓展到1kb
  site <- filter_fimo_to_1000bp_gr(motif.ovl)  
  # 4.5 与突变位置取交集
  # 与突变位置取交集
  mut.table<-mut_pos_tb(site,somatic_gastric)
  # warning(nrow(mut.table))
  
  # 画1k附近突变分布图
  x.labs <- as.character(mapply(paste0, str_replace(target_tf, ".meme", ""), "(", str_replace(co_factor, ".meme", ""), ")"))

  
  p_1 <- plot_mutation_loci(mut.table, center_expend_each = 10,
                          tf_name = str_replace(target_tf, ".meme", "") ,
                          mutationScore_outdir="./")
  p_1 <- p_1 + xlab(paste0(x.labs," + flanki n g region (bp)"))
  p_1
  ggsave("mutation_distribution_1k.pdf",p_1)

  # 画flank_5bp 区间分布图
  
  p_2 <- plot_flank_5bp(motif.ovl,somatic_gastric)+  # labs(caption=str_replace(target_tf,".meme", ""))+
          xlab(paste0(x.labs," + flanking region (bp)"))
  ggsave("mutation_distribution_5bp.pdf",p_2)
  # 回到主目录
  # setwd(paste0(singel_TF_mutation_dir, "/.."))

  # 清理环境变量中“cofactor.motif”
  rm(list=c("cofactor.motif","motif.ovl"))  
#}
}


########################################### 运行
require(parallel)
# 单核计算
system.time({
res <- lapply(seq(ncol(TF_Permute_matrix)), RunPairsTF);
})
lapply(1:3,RunPairsTF)

# 多核并行计算
detectCores(logical = F)  # 8
mc <- getOption("mc.cores", 23)
system.time({
res <- mclapply(seq(ncol(TF_Permute_matrix)), RunPairsTF)
})
