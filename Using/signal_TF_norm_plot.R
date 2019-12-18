#!/bin/R
library("tidyverse")
library(ggrepel)

# setwd("/public/home/kcao/Desktop/2018_NC_GS")
args=commandArgs(T)
file=args[1]


# 画图
# file=merge_mutation_Enrichment_Score
data=read_delim(file,col_names=FALSE,delim="\t")
colnames(data)[1:3]=c("TF_name","sep","score")
data <-data %>% mutate(color=case_when(score >4 ~ "high",2<score & score <4 ~"middle",TRUE ~"low" ))
# str_replace()
data <-data %>%mutate(label=ifelse(score >4,str_replace(data$TF_name,"_MA.*.meme",""),NA))
data$color=factor(data$color,levels=c("high","middle","low"))
P_1 <-ggplot(data,aes(x=TF_name,score))+
	geom_point(aes(col=data$color))+
	# theme(axis.text.x = element_text(angle = 90, hjust = 1))  # 旋转label 方向
	theme_classic()+
	ylim(0.5,6)+
	theme(#axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
	scale_colour_manual(values = c("red","orange", "black"),
						name="Mutation_score",
						labels=c("score >4", "2<score<4","1<score<2" ))+
	# geom_text(label=data$label,check_overlap = TRUE,size = 2,hjust = 0, nudge_x = 0.05)+
	#geom_label(label=data$label)+
	geom_label_repel(label=data$label)+
	ggtitle("2018 GS & 683 Jasper motif")
	
ggsave("merge_mutation_Enrichment_Score.pdf",P_1)
