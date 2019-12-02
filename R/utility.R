
#' Title
#'
#' @param fimo_file FIMO_TSV.file
#' @param as output type
#' @param ... other
#'
#' @return data_strucure such Grange or bed
#' @import IRanges GenomicRanges rtracklayer tidyverse
#' @importFrom GenomicRanges GRanges
#' @importFrom readr read_delim
#' @export
#'
#' @examples
Fimo2GRanges<-function(fimo_file,as="GRanges",...){
  require(rtracklayer)
  require(tidyverse)
  as <- match.arg(as, c("GRanges", "data.frame"))  # 模糊匹配，as="G"
  ctcf.motif=read_delim(fimo_file,delim ="\t")
  if (as=="GRanges"){
    ctcf.Fimo.gr<-GRanges(seqnames=ctcf.motif$sequence_name,
                          IRanges(start=ctcf.motif$start,end=ctcf.motif$stop),
                          pval=ctcf.motif$`p-value`,
                          qval=ctcf.motif$`q-value`,
                          dir=ctcf.motif$strand,
                          motif=ctcf.motif$matched_sequence)
    return(ctcf.Fimo.gr)
  }
  return(ctcf.motif)   # data.frame
}





#' bed2GRanges : read bed to GRanges
#'
#' @param header
#' @param BEDFile bed file
#'
#' @return bed information, in GRanges object
#' @import IRanges GenomicRanges tidyverse
#' @importFrom readr read_delim
#' @importFrom GenomicRanges GRanges
#' @export
#' @examples
#' BEDFile="./test_data/E094_DNase.all.peaks.v2.bed"
#' b2G=bed2GRanges(BEDFile,header = FALSE)
#' head(b2G)


bed2GRanges <-function(BEDFile,header=FALSE){
  require(tidyverse)
  tmp.DF<-read_delim(BEDFile,delim ="\t",col_names = header)
  tmp.GR<-GRanges(seqnames=tmp.DF[[1]],
                       IRanges(start=tmp.DF[[2]],end=tmp.DF[[3]]))
  cn <- colnames(tmp.DF)
  if (length(cn) > 3) {
    for (i in 4:length(cn)) {
      mcols(tmp.GR)[[cn[i]]] <- tmp.DF[[cn[i]]]
    }
  }
  return(tmp.GR)
}



