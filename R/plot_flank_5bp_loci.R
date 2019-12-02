#' plot mutation direction in the region 0f 5bp
#'
#' @param motif.ovl oberlaped motif granges
#' @param somatic_COLO829 snp Granges
#' @return graph
#' @import ggplot2
#' @import IRanges GenomicRanges rtracklayer
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @importFrom stats aggregate
#' @import Biostrings
#' @importFrom Biostrings reverseComplement
#' @importClassesFrom Biostrings DNAString
#' @export
#' @examples
#'
plot_flank_5bp<-function(motif.ovl,somatic_COLO829){
  motifs=GRanges(seqnames=motif.ovl$seqnames,
                 IRanges(motif.ovl$start,motif.ovl$end),
                 dir=motif.ovl$dir,
                 pval=motif.ovl$pval,
                 motif=motif.ovl$motif)
  motifs=motifs+5
  ## overlap
  z=findOverlaps(somatic_COLO829,motifs)
  t=as.data.frame(somatic_COLO829[queryHits(z)]) # 43
  t$motif.seq=as.character(seqnames(motifs[subjectHits(z)]))
  t$motif.start=start(motifs[subjectHits(z)])
  t$motif.end=end(motifs[subjectHits(z)])
  t$dir=motifs[subjectHits(z)]$dir
  t$pval=motifs[subjectHits(z)]$pval
  t$motif=motifs[subjectHits(z)]$motif


  t=t[order(t$id,t$pval),]
  t=t[!duplicated(t$id),]
  t.pos=t[which(t$dir=="+"),]
  t.neg=t[which(t$dir=="-"),]

  ###### pos
  ref.seq.pos=GRanges(seqnames=t.pos$seqnames,
                      IRanges(t.pos$motif.start,t.pos$motif.end),
                      dir=t.pos$dir,
                      motif=t.pos$motif,
                      mut=t.pos$start,
                      alt=t.pos$Alt,   ## FOXME
                      ref=t.pos$Ref)	 ## FOXME
  dnaSeq.pos=as.data.frame(ref.seq.pos)

  dnaSeq.mut.pos=dnaSeq.pos
  dnaSeq.mut.pos$alt=as.character(dnaSeq.mut.pos$alt)
  dnaSeq.mut.pos$num=dnaSeq.mut.pos$mut-dnaSeq.mut.pos$start+1


  ##### neg
  ref.seq.neg=GRanges(seqnames=t.neg$seqnames,
                      IRanges(t.neg$motif.start,t.neg$motif.end),
                      dir=t.neg$dir,
                      motif=t.neg$motif,
                      mut=t.neg$start,
                      alt=t.neg$Alt,
                      ref=t.neg$Ref)
  dnaSeq.neg=as.data.frame(ref.seq.neg)

  dnaSeq.mut.neg=dnaSeq.neg
  dnaSeq.mut.neg$alt=as.character(dnaSeq.mut.neg$alt)
  dnaSeq.mut.neg$num=dnaSeq.mut.neg$end-dnaSeq.mut.neg$mut+1


  ### plot
  df=rbind(dnaSeq.mut.pos[,c("ref","alt","num")],dnaSeq.mut.neg[,c("ref","alt","num")])
  df$count=1
  df=aggregate(count~ref+alt+num,df,sum)

  require(Biostrings)
  for (i in 1:nrow(df)){
    if (!df$ref[i] %in% c("C","T")){
      df$alt[i]=as.character(reverseComplement(Biostrings::DNAString(df$alt[i])))
      df$ref[i]=as.character(reverseComplement(Biostrings::DNAString(df$ref[i])))
    }
  }

  df$mut=paste(df$ref,df$alt,sep="")
  df$mut=factor(df$mut,levels=c("TG","TC","TA","CT","CG","CA"))
  df$num=factor(df$num,levels=c(1:29))

  df=aggregate(count~ref+alt+num+mut,df,sum) # 20

  # patchwork

  P_5<-ggplot(df,aes(x=num,y=count,fill=mut))+
    geom_bar(stat="identity")+
    theme_classic()+
    scale_fill_manual(values = c("#33CC00","#33CCFF","#FF9933","#CC33FF","#FFFF33","#FF0000"))+
    scale_x_discrete(breaks=seq(1,29))+
    #scale_y_discrete(breaks=seq(1,10))+
    theme(legend.position="top")
  return(P_5)
}
