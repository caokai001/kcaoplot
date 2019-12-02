ctcf.motif=Fimo2GRanges("../TESTDATA/mynew.txt","GRanges")
# 开放区域合并
dnase.peak<-bed2GRanges("../TESTDATA/E094-DNase.all.peaks.v2.bed")
dnase.fdr<-bed2GRanges("../TESTDATA/E094-DNase.fdr0.01.peaks.v2.bed")
dnase.mac<-bed2GRanges("../TESTDATA/E094-DNase.macs2.narrowPeak")
dnase=c(dnase.peak,dnase.fdr,dnase.mac)
dnase=GenomicRanges::reduce(dnase)
# 读入snp信息
somatic_gastric<-bed2GRanges("../TESTDATA/41467_2018_3828_MOESM6_ESM.txt",header=TRUE)
names(mcols(somatic_gastric))[1:2]<-c("Ref","Alt")
somatic_gastric$id=c(1:length(somatic_gastric))
# 过滤FIMO
motif.ovl<-filter_fimo(dnase,ctcf.motif)
# 拓展到1kb
site=filter_fimo_to_1000bp_gr(motif.ovl)
# 与突变位置取交集
mut.table<-mut_pos_tb(site,somatic_gastric)
# 画1kp 区间突变分布图
#plot_mutation_loci(mut.table)
plot_mutation_loci(mut.table,center_expend_each=5,tf_name="NA",mutationScore_outdir="./")
