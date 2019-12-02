#' filter fimo_CTCF.
#'
#' @param dnase granges objects
#' @param fimo_GRanges ctcf_motif granges objects
#' @return filtered fimo GRanges objects
#' @importFrom  GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @export
#' @examples
#'
filter_fimo<-function(dnase,fimo_GRanges){
  z<-findOverlaps(dnase,fimo_GRanges)
  dnase.ovl=dnase[queryHits(z)]
  dnase.ovl=as.data.frame(dnase.ovl)
  dnase.ovl=unique(dnase.ovl)
  motif.ovl=fimo_GRanges[subjectHits(z)]
  motif.ovl=as.data.frame(motif.ovl)
  motif.ovl=unique(motif.ovl)
  return(motif.ovl)
}



#' extension 1000bp about center motif.
#'
#' @param motif.ovl
#'
#' @return site 1000 bp region
#' @importFrom  GenomicRanges GRanges
#' @export
#'
#' @examples
#'
filter_fimo_to_1000bp_gr<-function(motif.ovl){
  motif_center_site<-GRanges(seqnames=motif.ovl$seqnames,
          IRanges(start=(motif.ovl$start+motif.ovl$end)/2,end=(motif.ovl$start+motif.ovl$end)/2),
          pval=motif.ovl$pval,
          dir=motif.ovl$dir,
          motif=motif.ovl$motif)
  site=motif_center_site+1000
  return(site)
}


#' output position in 1000bp flank region.
#'
#' @param site flank-region granges objects
#' @param somatic_gastric granges objects
#'
#' @return mutation_pos_table
#' @importFrom  GenomicRanges findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom S4Vectors subjectHits
#' @export
#'
#' @examples
#'
mut_pos_tb<-function(site,somatic_gastric){
  mut=findOverlaps(site,somatic_gastric)
  mut.table=somatic_gastric[subjectHits(mut)]
  mut.table$region.start=start(site[queryHits(mut)])
  mut.table$region.end=end(site[queryHits(mut)])
  mut.table$region.chr=seqnames(site[queryHits(mut)])
  mut.table$pval=site[queryHits(mut)]$pval
  mut.table$dir=site[queryHits(mut)]$dir
  mut.table$motif=site[queryHits(mut)]$motif
  mut.table$id=c(1:length(mut.table))
  ###########
  mut.table=as.data.frame(mut.table)
  mut.table=mut.table[order(mut.table$id,mut.table$pval),]
  mut.table=mut.table[!duplicated(mut.table$id),]
  mut.table$seqnames=as.character(mut.table$seqnames)
  mut.table$region.chr=as.character(mut.table$region.chr)
  mut.table$motif.start=(mut.table$region.start+mut.table$region.end)/2
  return(mut.table)
}




