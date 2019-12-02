#' 1000bp mutation loci
#'
#' @param mut.table mutation position
#' @param center_expend_each expand center region to statictic enrichment
#' @param tf_name Transcription factor name
#' @param mutationScore_outdir enrichment score output directory
#' @return
#' @import ggplot2
#' @export
#' @examples
#'
plot_mutation_loci<-function(mut.table,center_expend_each=5,tf_name="NA",mutationScore_outdir="./"){
  mut.table.pos=mut.table[which(mut.table$dir=="+"),]
  mut.table.neg=mut.table[which(mut.table$dir=="-"),]

  ## statistic postive mut number
  mut.table.pos$pos=mut.table.pos$start-mut.table.pos$motif.start
  mut.table.pos.df=aggregate(seqnames~pos,mut.table.pos,length)
  colnames(mut.table.pos.df)[2]="count"
  df=data.frame(pos=seq(-1000,1000,1),x=0)
  mut.table.pos.df=merge(df,mut.table.pos.df,all.x=TRUE,by="pos")
  mut.table.pos.df$count=ifelse(is.na(mut.table.pos.df$count),0,mut.table.pos.df$count)

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




  #### somatic enrichment score :FoldChange
  require(tidyverse)
  motif_len <- unique(nchar( head(mut.table$motif)))  # motif 长度
  center_region_left <- 1001-ceiling(motif_len/2+center_expend_each)
  center_region_right <- 1001+ceiling(motif_len/2+center_expend_each)
  # calculate background
  background=mut.table.z$norm.mut[-c(center_region_left:center_region_right)]%>%mean()
  # calculate center enrichment score

  center_vector <- mut.table.z$norm.mut[c(center_region_left:center_region_right)]
  top_score <- sapply(center_vector, function(x){return(x/background)}) %>%
    sort(decreasing = TRUE) %>%
    head(1)
  # save mutation score
  cat(paste0(tf_name," : ",top_score),
      file=paste0(mutationScore_outdir,"mutation_Enrichment_Score.txt"),
      append = T,
      sep = "\n")



  # graph mut in 1kb region
  p<-ggplot(data=mut.table.z,aes(x=pos,y=norm.mut))+
    geom_line()+# ylim(c(0,10))
    ggtitle(paste0("mut.table\n","n=",nrow(mut.table)))+
    xlab("CBS + flanking region (bp)")+
    ylab("somatic mutation nunber")+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line=element_line(colour="black"))

  return(p)

}
