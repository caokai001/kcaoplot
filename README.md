# kcaoplot
R package for intersting🐭

>### install
```
library(devtools)
install_github("caokai001/kcaoplot")
library(kcaoplot)
```
---
# Part1 : try R package
> 尝试使用R包，写小函数🐻

### example 1

test_ggplot(6)

![2019年11月23日](https://upload-images.jianshu.io/upload_images/9589088-c0ea6db681989e36.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
- test_plot(4)
- hello()

>### contact: 
>[简书](https://www.jianshu.com/p/b5bc50d4c1d0)

---




# Part2: mutation analysis

> 目的：

​	从癌症组织中snp.vcf, 可以猜想snp 分布也许不是随机的，于是可以尝试研究是否到某个motif附近富集。

### example 2





> input file:

- ​          mynew.txt            :    Fimo.tsv file;use FIMO to scan hg19.fa
- NHEK_merge_ctcf.bed :    CTCF peak 文件或者Dnase peak 文件
- hg19_uv_snp_pos.bed ： snp 文件



```R
# UV
ctcf.motif=Fimo2GRanges("mynew.txt","GRanges")
NHEK.peak=bed2GRanges("./ctcf_peak/NHEK_merge_ctcf.bed",header=FALSE)
# snptogr
somatic_COLO829=bed2GRanges("hg19_uv_snp_pos.bed",header=FALSE)
names(mcols(somatic_COLO829)) <-c("Ref","Alt")
somatic_COLO829$id=c(1:length(somatic_COLO829))
# filter_fimo
motif.ovl=filter_fimo(NHEK.peak,ctcf.motif)
# 拓宽到1000bp
filter_fimo_to_1000bp_gr(motif.ovl)
# snp 落到区间位置
mut.table<-mut_pos_tb(site,somatic_COLO829)
# 画图1000bp
plot_mutation_loci(mut.table)
# 画图flank 5bp
plot_flank_5bp(motif.ovl,somatic_COLO829)

```

![](https://upload-images.jianshu.io/upload_images/9589088-d4adfa1ca0722af2.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

### example 3

```R
# 读入ctcf.motif
ctcf.motif=Fimo2GRanges("mynew.txt","GRanges")
# 开放区域合并
dnase.peak<-bed2GRanges("E094-DNase.all.peaks.v2.bed")
dnase.fdr<-bed2GRanges("E094-DNase.fdr0.01.peaks.v2.bed")
dnase.mac<-bed2GRanges("E094-DNase.macs2.narrowPeak")
dnase=c(dnase.peak,dnase.fdr,dnase.mac)
<<<<<<< HEAD
=======
library(GenomicRanges)
>>>>>>> d47606c6f9c4fcbfee5d37b391c8c5fe3d5fe782
dnase=reduce(dnase)
# 读入snp信息
somatic_gastric<-bed2GRanges("41467_2018_3828_MOESM6_ESM.txt",header=TRUE)
names(mcols(somatic_gastric))[1:2]<-c("Ref","Alt")
somatic_gastric$id=c(1:length(somatic_gastric))
# 过滤FIMO
motif.ovl<-filter_fimo(dnase,ctcf.motif)
# 拓展到1kb
site=filter_fimo_to_1000bp_gr(motif.ovl)
# 与突变位置取交集
mut.table<-mut_pos_tb(site,somatic_gastric)
# 画1kp 区间突变分布图
plot_mutation_loci(mut.table)
# 画flank_5bp 区间分布图
plot_flank_5bp(motif.ovl,somatic_gastric)
```

<img src="https://upload-images.jianshu.io/upload_images/9589088-45f2283c607b2157.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240" style="zoom: 67%;" />


