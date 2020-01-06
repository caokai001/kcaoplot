#! bin/bash
mynew_FIMO="/public/home/kcao/Desktop/Cancer_cell_2019/motif_dir/AR_Moods/hg19_AR_motif_2fimo.tsv"
DNase_macs2_narrowPeak="/public/home/kcao/genome_human/chain_convert/h3k27ac_filter_AR_motif.bed"

MOESM6_ESM_SNP="/public/home/kcao/Desktop/TCGA_prostate/DATA_SNP/ICGC_within_SI_Data_1_filtered_variants_by_patient.bed"
TF_NAME="AR_motif"

## 1.modes_motif 转换成bed
(atac) [kcao@login Bedtools_result]$ cat /public/home/kcao/Desktop/Cancer_cell_2019/motif_dir/AR_Moods/hg19_AR_motif_moods.csv|
awk 'BEGIN{FS=",";OFS="\t"}{print $1,$3+1,$3+15,$2,$5,$4,$6}' >moods_AR_hg19.bed

## 2.h3k27ac_filter_AR_motif.bed 里面扫描AR motif
(atac) [kcao@login Bedtools_result]$ intersectBed -a /public/home/kcao/Desktop/Cancer_cell_2019/SCRIPT/Bedtools_result/moods_AR_hg19.bed -b \
/public/home/kcao/genome_human/chain_convert/h3k27ac_filter_AR_motif.bed -wa >AR_motif_filter_h3k27ac_ARpeak.bed
# 16096


## 3.对过滤后的AR motif 与snp 取交集
# 3.0 基因组size
[kcao@login hg19]$ cat hg19.fa.fai |grep -v "_"|cut -f1-2|sort -Vk 1 >hg19.genome.size

# 3.1 延伸motif 500bp，作为画图区域 (将bs 设置为10bp )
[kcao@login Bedtools_result]$ bedtools slop -i /public/home/kcao/Desktop/Cancer_cell_2019/SCRIPT/Bedtools_result/AR_motif_filter_h3k27ac_ARpeak.bed -g
 /public/home/kcao/genome_human/hg19/hg19.genome.size -b 500 >AR_motif_filtered_flank_50bp.bed
 
# 3.2 strands(+/-) 计算相对位置
[kcao@login Bedtools_result]$ sed '1d' /public/home/kcao/Desktop/TCGA_prostate/DATA_SNP/ICGC_within_SI_Data_1_filtered_variants_by_patient.bed |intersectBed \
-a stdin -b AR_motif_filtered_flank_50bp.bed -wa -wb|awk '{print $0"\t"($7+$8)/2}'

# 3.3 突变相对位置
[kcao@login Bedtools_result]$ sed '1d' /public/home/kcao/Desktop/TCGA_prostate/DATA_SNP/ICGC_within_SI_Data_1_filtered_variants_by_patient.bed |
intersectBed -a stdin -b AR_motif_filtered_flank_50bp.bed -wa -wb|awk '{print $0"\t"($7+$8)/2}'|
awk '{if($11=="+"){print $0"\t"($2-$13)}else{print $0"\t"($13-$2)}}' >AR_mutation_relative_position.bed

# 3.4 输出画图数据
[kcao@login Bedtools_result]$ cat AR_mutation_relative_position.bed |cut -f14 >plot_pos_data_AR_500bp.txt

# 3.5 画图.R
>>> radian
#################### 
# bedtools
# Tips: 未添加突变方向

library(tidyverse)
library(ggthemes)
library(RColorBrewer)
# display.brewer.all(type = "all")

# 过滤区间
count_data <- read.table(file = "https://raw.githubusercontent.com/caokai001/Bioinformatic_code/master/%E8%AF%BE%E9%A2%98%E7%9B%B8%E5%85%B3/%E5%A4%87%E4%BB%BD%E6%95%B0%E6%8D%AE%E5%A4%B9/plot_pos_data_AR_500bp.txt",header= F)
head(count_data)
tmp <-data.frame(pos=as.numeric((names(table(count_data$V1)))) ,count=as.numeric((table(count_data$V1))))
tmp<-tmp[which(abs(tmp$pos) <=500),][-1,]
# 分bin或者分组10bp
head(tmp)
tmp$group <- ceiling(tmp$pos/10)
head(tmp)

T <- tbl_df(tmp) %>% group_by(group) %>% summarize(sumCT=sum(count)) %>% ungroup()

ggplot(T,aes(x=group,y=sumCT))+
  #geom_bar(stat = "identity")+
  geom_line(color=brewer.pal(7, "Set1")[2],size=1.5)+
  # 添加线条图层(dashed line)
  geom_vline(xintercept=c(-0.7,0.7),linetype="dashed")+
  # 添加坐标轴线条
  # theme(axis.line.x.bottom = element_blank() )+
  geom_hline(yintercept = -0.5)+
  ggtitle("AR+NR")+  # <----
  ylim(0,max(T$sumCT))+
  xlab("")+
  ylab("Number of Mutation Per Bin")+
  scale_x_continuous(breaks = c(-50,0,50),
                     labels = c("-500","Motif","500"))+
  # 修改主题
  theme_classic()+
  theme(axis.text.x = element_text(size=15,face = "bold"),
        plot.title = element_text(size=15,face = "bold"),
        axis.title.y = element_text(size=15,face = "bold"))

# 画图结果链接
https://upload-images.jianshu.io/upload_images/9589088-658f16bed556e4a7.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240

## 4.添加6种突变方向
# 4.0 基因组size
[kcao@login hg19]$ cat hg19.fa.fai |grep -v "_"|cut -f1-2|sort -Vk 1 >hg19.genome.size

# 4.1 延伸motif 500bp，作为画图区域 (将bs 设置为10bp )
[kcao@login Bedtools_result]$ bedtools slop -i /public/home/kcao/Desktop/Cancer_cell_2019/SCRIPT/Bedtools_result/AR_motif_filter_h3k27ac_ARpeak.bed -g
 /public/home/kcao/genome_human/hg19/hg19.genome.size -b 500 >AR_motif_filtered_flank_50bp.bed
 
# 4.2 strands(+/-) 计算相对位置
pass

# 4.3 突变相对位置(***)
# 考虑正负链影响相对位置
# 突变方向需要单独校正,不考虑方向
##### 考虑碱基大小写，awk:toupper

[kcao@login Bedtools_result]$ sed '1d' /public/home/kcao/Desktop/TCGA_prostate/DATA_SNP/ICGC_within_SI_Data_1_filtered_variants_by_patient.bed |
intersectBed -a stdin -b AR_motif_filtered_flank_50bp.bed -wa -wb|
awk '{print $0"\t"($7+500)"\t"$8-500}' |
awk '{ if($11=="+" && ($2<$13)) print $0"\t"($2-$13)
		else if($11=="+" && ($2>$14)) print $0"\t"($2-$14)
		else if($2>= $13  && $2<=$14) print $0"\t0"
		else if($11=="-" && ($2<$13)) print $0"\t"($13-$2) 
		else print $0"\t"($14-$2)}' |
awk '{if(toupper($4)=="G"    && toupper($5)=="A") print $0"\tCT"
	else if(toupper($4)=="G" && toupper($5)=="C") print $0"\tCG"
	else if(toupper($4)=="G" && toupper($5)=="T") print $0"\tCA"
	else if(toupper($4)=="A" && toupper($5)=="C") print $0"\tTG"
	else if(toupper($4)=="A" && toupper($5)=="G") print $0"\tTC"
	else if(toupper($4)=="A" && toupper($5)=="T") print $0"\tTA"
	else print $0"\t"toupper($4)toupper($5)}' >AR_mutation_relative_position_direction.bed

# 4.4 输出画图数据
[kcao@login Bedtools_result]$ cat AR_mutation_relative_position_direction.bed |cut -f 15-16 >plot_pos_data_AR_500bp_direction.txt


# 4.5 画图.R 
>>> radian
library(tidyverse)
library(ggthemes)
library(RColorBrewer)

# 读取数据
count_data <- read.table(file ="https://raw.githubusercontent.com/caokai001/Bioinformatic_code/master/%E8%AF%BE%E9%A2%98%E7%9B%B8%E5%85%B3/%E5%A4%87%E4%BB%BD%E6%95%B0%E6%8D%AE%E5%A4%B9/plot_pos_data_AR_500bp_direction.txt",header=F)

colnames(count_data)=c("Pos","Mut")
# count_data$Pos=factor(count_data$Pos,levels=c(-500:500))
count_data$Mut=factor(count_data$Mut,levels=c("TG","TC","TA","CT","CG","CA"))
#View(count_data)

# 统计每个位置，突变情况
df <-count_data %>%tbl_df %>%group_by(Pos,Mut) %>%summarize(count=n()) %>%ungroup()

# View(df)

# 未分组画图(中心高度除以1.5)
ggplot(df,aes(x=Pos,y=count,fill=Mut))+
  geom_bar(stat="identity",colour="black")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black"))+
  scale_fill_manual(values = c("#33CC00","#33CCFF","#FF9933","#CC33FF","#FFFF33","#FF0000"))





### 分10bp组，画堆积图
# 按照bins=10bp 分组画图

# 分大于0 和小于0 ,注意数据处理不要是因子factor
df$group <- sapply(df$Pos, function(x){
  if(x>0){
    return(ceiling(x/10))}
  else {
	return(-ceiling(abs(x)/10))  }}
	)
	

Ts <- tbl_df(df) %>% group_by(group,Mut) %>% summarize(sumCT=sum(count)) %>% ungroup()
View(Ts)
# tips: 中间区域15bp,其他区域10bp,需要除以1.5倍差异
Ts$sumCT=as.numeric(Ts$sumCT)
Ts$sumCT= ifelse(Ts$group==0,round(Ts$sumCT/1.5,2),Ts$sumCT)

Ts$group=as.numeric(Ts$group)
Ts$Mut=factor(paste0(substr(Ts$Mut,1,1),">",substr(Ts$Mut,2,2)))
Ts$Mut =factor(Ts$Mut,levels = c("T>G","T>C","T>A","C>T","C>G","C>A"))


# 填补空缺值：避免异常的波动，缺失值按照0 填充.
Impute <- data.frame(
  group=as.numeric(rep(-50:50,each=6)),
  Mut=factor(rep(c("T>G","T>C","T>A","C>T","C>G","C>A"),times=101)),
  sumCT=0)
Ts <- merge(data.frame(Ts),Impute,by=c("group","Mut"),all.y=TRUE) %>%tbl_df() %>% 
  transmute(group,Mut,sumCT=ifelse(is.na(sumCT.x),sumCT.y,sumCT.x)) 



# 确保x,y 都是numeric 类型
ggplot(Ts,aes(x=group,y=sumCT,fill=Mut))+geom_area(position = "stack")+
  scale_fill_manual(values = c("#ffc1c1","#66cd00","#bebebe","#ee2c2c","#252525","#00b2ee"))+
  # 添加线条图层(dashed line)
  geom_vline(xintercept=c(-0.7,0.7),linetype="dashed",size=0.5)+
  ggtitle("AR+NR")+  # <---
  xlab("")+
  ylab("Number of Mutation Per Bin")+
  scale_x_continuous(breaks = c(-50,0,50),
                     labels = c("-500","Motif","500"))+
  # 修改图例顺序及其位置，去掉title
  guides(fill = guide_legend(reverse=TRUE,title=NULL)) +
  # scale_fill_discrete(labels=c("C>A","C>G","C>T","T>A","T>C","T>G"))+
  # 修改主题
  theme_classic()+
  theme(axis.text.x = element_text(size=15,face = "bold"),
        plot.title = element_text(size=15,face = "bold"),
        axis.title.y = element_text(size=15,face = "bold"))


####
图片效果：
https://upload-images.jianshu.io/upload_images/9589088-83e4ad6f5ef981b1.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240
