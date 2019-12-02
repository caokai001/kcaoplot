############ UV

ctcf.motif=Fimo2GRanges("mynew.txt","GRanges")
NHEK.peak=bed2GRanges("./ctcf_peak/NHEK_merge_ctcf.bed",header=FALSE)
# snptogr
somatic_COLO829=bed2GRanges("hg19_uv_snp_pos.bed",header=FALSE)
names(mcols(somatic_COLO829)) <-c("Ref","Alt")
somatic_COLO829$id=c(1:length(somatic_COLO829))
# filter_fimo
motif.ovl=filter_fimo(NHEK.peak,ctcf.motif)
# 拓宽到1000bp
site=filter_fimo_to_1000bp_gr(motif.ovl)
# snp 落到区间位置
mut.table<-mut_pos_tb(site,somatic_COLO829)
# 画图1000bp
plot_mutation_loci(mut.table)
# 画图flank 5bp
plot_flank_5bp(motif.ovl,somatic_COLO829)

############ GS
ctcf.motif=Fimo2GRanges("mynew.txt","GRanges")
dnase.peak<-bed2GRanges("E094-DNase.all.peaks.v2.bed")
dnase.fdr<-bed2GRanges("E094-DNase.fdr0.01.peaks.v2.bed")
dnase.mac<-bed2GRanges("E094-DNase.macs2.narrowPeak")
dnase=c(dnase.peak,dnase.fdr,dnase.mac)
dnase=reduce(dnase)
somatic_gastric<-bed2GRanges("41467_2018_3828_MOESM6_ESM.txt",header=TRUE)
names(mcols(somatic_gastric))[1:2]<-c("Ref","Alt")
somatic_gastric$id=c(1:length(somatic_gastric))
motif.ovl<-filter_fimo(dnase,ctcf.motif)
site=filter_fimo_to_1000bp_gr(motif.ovl)
mut.table<-mut_pos_tb(site,somatic_gastric)
plot_mutation_loci(mut.table)
plot_flank_5bp(motif.ovl,somatic_gastric)


plot_flank_5bp(motif.ovl,somatic_gastric)