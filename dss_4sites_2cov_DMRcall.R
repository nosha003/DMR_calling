dss_4sites_2cov_DMRcall <- function(input_dir1, input_dir2, output_dir, file1, file2, genotype1, genotype2, dataset) {
setwd(input_dir1)
df1 <- read.delim(file1, header=T, sep="\t", stringsAsFactors = F)
df1$start <- df1$start - 1
setwd(input_dir2)
df2 <- read.delim(file2, header=T, sep="\t", stringsAsFactors = F)
df2$start <- df2$start - 1 

library(dplyr)
df1.1 <- mutate_all(df1, function(x) as.numeric(as.character(x)))
df2.1 <- mutate_all(df2, function(x) as.numeric(as.character(x)))

df1_dss_cg <- subset(df1.1, select=c(1,2,6,5))
colnames(df1_dss_cg) <- c("chr", "pos", "N", "X")
df2_dss_cg <- subset(df2.1, select=c(1,2,6,5))
colnames(df2_dss_cg) <- c("chr", "pos", "N", "X")
df1_dss_chg <- subset(df1.1, select=c(1,2,11,10))
colnames(df1_dss_chg) <- c("chr", "pos", "N", "X")
df2_dss_chg <- subset(df2.1, select=c(1,2,11,10))
colnames(df2_dss_chg) <- c("chr", "pos", "N", "X")
df1_dss_chh <- subset(df1.1, select=c(1,2,16,15))
colnames(df1_dss_chh) <- c("chr", "pos", "N", "X")
df2_dss_chh <- subset(df2.1, select=c(1,2,16,15))
colnames(df2_dss_chh) <- c("chr", "pos", "N", "X")

df1_dss_cg2 <- subset(df1_dss_cg, df1_dss_cg$N >= df1_dss_cg$X)
df2_dss_cg2 <- subset(df2_dss_cg, df2_dss_cg$N >= df2_dss_cg$X)
df1_dss_chg2 <- subset(df1_dss_chg, df1_dss_chg$N >= df1_dss_chg$X)
df2_dss_chg2 <- subset(df2_dss_chg, df2_dss_chg$N >= df2_dss_chg$X)
df1_dss_chh2 <- subset(df1_dss_chh, df1_dss_chh$N >= df1_dss_chh$X)
df2_dss_chh2 <- subset(df2_dss_chh, df2_dss_chh$N >= df2_dss_chh$X)

library(DSS)
require(bsseq)
library(dplyr)
library(tidygenomics)

setwd(output_dir)

bp_BSobj <- makeBSseqData(list(df1_dss_cg2, df2_dss_cg2), c("sample_1", "sample_2"))
bp_dmlTest <- DMLtest(bp_BSobj, group1=c("sample_1"), group2=c("sample_2"), smoothing=FALSE, equal.disp=TRUE)
write.table(bp_dmlTest, "BM_5genos_DSS_dmlTest_allchr_CG.txt", sep="\t", quote=F, row.names=F)
bp_dml_cg <- callDML(bp_dmlTest, delta=0.1, p.threshold=0.001)
write.table(bp_dml_cg, paste(genotype1, genotype2, "_", dataset, "_DSS_dml_nosmoothing_allchr_CG.txt", sep=""), sep="\t", quote=F, row.names=F)

bp_BSobj <- makeBSseqData(list(df1_dss_chg2, df2_dss_chg2), c("sample_1", "sample_2"))
bp_dmlTest <- DMLtest(bp_BSobj, group1=c("sample_1"), group2=c("sample_2"), smoothing=FALSE, equal.disp=TRUE)
write.table(bp_dmlTest, "BM_5genos_DSS_dmlTest_allchr_CHG.txt", sep="\t", quote=F, row.names=F)
bp_dml_chg <- callDML(bp_dmlTest, delta=0.1, p.threshold=0.001)
write.table(bp_dml_chg, paste(genotype1, genotype2, "_", dataset, "_DSS_dml_nosmoothing_allchr_CHG.txt", sep=""), sep="\t", quote=F, row.names=F)

bp_BSobj <- makeBSseqData(list(df1_dss_chh2, df2_dss_chh2), c("sample_1", "sample_2"))
bp_dmlTest <- DMLtest(bp_BSobj, group1=c("sample_1"), group2=c("sample_2"), smoothing=FALSE, equal.disp=TRUE)
write.table(bp_dmlTest, "BM_5genos_DSS_dmlTest_allchr_CHH.txt", sep="\t", quote=F, row.names=F)
bp_dml_chh <- callDML(bp_dmlTest, delta=0.1, p.threshold=0.001)
write.table(bp_dml_chh, paste(genotype1, genotype2, "_", dataset, "_DSS_dml_nosmoothing_allchr_CHH.txt", sep=""), sep="\t", quote=F, row.names=F)


# post-hoc filtering
library(dplyr)
library(ggplot2)

data <- full_join(df1.1, df2.1, by=c("chr", "start", "end"))

setwd(output_dir)
dss_cg <- bp_dml_cg
dss_chg <- bp_dml_chg
dss_chh <- bp_dml_chh

dss_cg <- dss_cg[,1:2]
colnames(dss_cg) <- c("chr", "start")
dss_cg$end <- dss_cg$start + 100
dss_data_cg <- left_join(dss_cg, data[,c(1:8,19:23)], by=c("chr", "start"))
dss_data_cg$diff <- abs(dss_data_cg$CG_ratio.x - dss_data_cg$CG_ratio.y)
dss_sites_cg <- subset(dss_data_cg, dss_data_cg$CG_sites.x >= 4 & dss_data_cg$CG_sites.y >= 4 & dss_data_cg$CG_cov.x >= 2 & dss_data_cg$CG_cov.y >= 2)
write.table(dss_sites_cg, paste(genotype1, genotype2, "_", dataset, "_DSS_dml_nosmoothing_allchr_CG_4sites_2cov.txt", sep=""), sep="\t", quote=F, row.names=F)

dss_chg <- dss_chg[,1:2]
colnames(dss_chg) <- c("chr", "start")
dss_chg$end <- dss_chg$start + 100
dss_data_chg <- left_join(dss_chg, data[,c(1:3,9:13,24:28)], by=c("chr", "start"))
dss_data_chg$diff <- abs(dss_data_chg$CHG_ratio.x - dss_data_chg$CHG_ratio.y)
dss_sites_chg <- subset(dss_data_chg, dss_data_chg$CHG_sites.x >= 4 & dss_data_chg$CHG_sites.y >= 4 & dss_data_chg$CHG_cov.x >= 2 & dss_data_chg$CHG_cov.y >= 2)
write.table(dss_sites_chg, paste(genotype1, genotype2, "_", dataset, "_DSS_dml_nosmoothing_allchr_CHG_4sites_2cov.txt", sep=""), sep="\t", quote=F, row.names=F)

dss_chh <- dss_chh[,1:2]
colnames(dss_chh) <- c("chr", "start")
dss_chh$end <- dss_chh$start + 100
dss_data_chh <- left_join(dss_chh, data[,c(1:3,14:18,29:33)], by=c("chr", "start"))
dss_data_chh$diff <- abs(dss_data_chh$CHH_ratio.x - dss_data_chh$CHH_ratio.y)
dss_sites_chh <- subset(dss_data_chh, dss_data_chh$CHH_sites.x >= 4 & dss_data_chh$CHH_sites.y >= 4 & dss_data_chh$CHH_cov.x >= 2 & dss_data_chh$CHH_cov.y >= 2)
write.table(dss_sites_chh, paste(genotype1, genotype2, "_", dataset, "_DSS_dml_nosmoothing_allchr_CHH_4sites_2cov.txt", sep=""), sep="\t", quote=F, row.names=F)
}


#dss_4sites_2cov_DMRcall <- function(input_dir1, input_dir2, file1, file2, genotype1, genotype2, dataset)

# GENOTYPE
#dss_4sites_2cov_DMRcall(
#  input_dir1 = "/scratch.global/nosha003/dmr/tiles_output",
#  input_dir2 = "/scratch.global/nosha003/dmr/tiles_output",
#  output_dir = "/scratch.global/nosha003/dmr/DMRcalls",
#  file1 = "B73_V2_100bp_cov.txt",
#  file2 = "PH207_V2_100bp_cov.txt",
#  genotype1 = "B73",
#  genotype2 = "PH207",
#  dataset = "genotype"
#)
