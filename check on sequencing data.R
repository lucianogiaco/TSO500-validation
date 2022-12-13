###SEQUENCING METRICS###

dati_paper_NGS<-read.csv2(file.choose(),sep=";",h=T,dec=",",na.strings="") #importa tab2_tab2.csv
head(dati_paper_NGS)
dim(dati_paper_NGS) #112 48
colnames(dati_paper_NGS)
summary(dati_paper_NGS)
table(dati_paper_NGS$run) #numero di run: 5

dati_paper_NGS[dati_paper_NGS$run2=="210715_A01423_0008_AH35CWDRXY",]
dati_paper_NGS[dati_paper_NGS$run2=="210729_A01423_0009_AH33WGDRXY",]
dati_paper_NGS[dati_paper_NGS$run2=="211022_A01423_0010_AHGYFYDRXY",]
dati_paper_NGS[dati_paper_NGS$run2=="211111_A01423_0011_AHH2Y2DRXY",]
dati_paper_NGS[dati_paper_NGS$run2=="211122_A01423_0012_AH2YWCDRXY",]

#IMPUTAZIONE DEL VALORE DI RNA
summary(dati_paper_NGS$rna_median_cv__gene_500x_perc)
dati_paper_NGS$rna_median_cv__gene_500x_perc[dati_paper_NGS$rna_median_cv__gene_500x_perc==117.08]<-100
summary(dati_paper_NGS$rna_median_cv__gene_500x_perc)

getwd()
#write.table(dati_paper_NGS,"dati_paper_NGS.csv",sep=";",dec=",",row.names=F,col.names=T,na="")



#########
###DNA###
#########

datiDNAsample<-data.frame(dati_paper_NGS$contamination_score,
dati_paper_NGS$contamination_p_value,
dati_paper_NGS$total_pf_reads,
dati_paper_NGS$median_target_coverage,
dati_paper_NGS$pct_chimeric_reads,
dati_paper_NGS$pct_exon_100x,
dati_paper_NGS$pct_read_enrichment,
dati_paper_NGS$pct_usable_umi_reads,
dati_paper_NGS$mean_target_coverage,
dati_paper_NGS$pct_aligned_reads,
dati_paper_NGS$pct_contamination_est,
dati_paper_NGS$pct_pf_uq_reads,
dati_paper_NGS$pct_target_04x_mean,
dati_paper_NGS$pct_target_100x,
dati_paper_NGS$pct_target_250x,
dati_paper_NGS$median_insert_size,
dati_paper_NGS$median_exon_coverage_count,
dati_paper_NGS$pct_exon_50x_perc,dati_paper_NGS$run)

dim(datiDNAsample) #112 19 
summary(datiDNAsample)

apply(datiDNAsample,2,function(x) sd(x,na.rm=T))



##########Kruskal-Wallis and Mann-Whitney###########

###Soglie per linee guida: metrics Output dal file
###trusight-oncology-500-local-app-v2.2-user-guide-1000000137777-01.pdf

options(scipen=999999999)

##############################
#dati_paper_NGS$pct_exon_100x#
##############################

summary(dati_paper_NGS$pct_exon_100x); #table(dati_paper_NGS$pct_exon_100x>...,useNA="ifany") 
sd(dati_paper_NGS$pct_exon_100x,na.rm=T)
tapply(dati_paper_NGS$pct_exon_100x,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$pct_exon_100x,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$pct_exon_100x,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$pct_exon_100x,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$pct_exon_100x~dati_paper_NGS$run, main="Exon 100x by run",ylab="Exon 100x (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$pct_exon_100x~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$pct_exon_100x, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$pct_exon_100x, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$pct_exon_100x, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$pct_exon_100x,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))

library(agricolae)
#Median.test(dati_paper_NGS$pct_exon_100x,factor(dati_paper_NGS$run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



library(coin)
wilcox_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_exon_100x[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_exon_100x)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr


#####################################
#dati_paper_NGS$pct_usable_umi_reads#
#####################################

summary(dati_paper_NGS$pct_usable_umi_reads); #table(dati_paper_NGS$pct_usable_umi_reads>...,useNA="ifany") 
sd(dati_paper_NGS$pct_usable_umi_reads,na.rm=T)
tapply(dati_paper_NGS$pct_usable_umi_reads,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$pct_usable_umi_reads,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$pct_usable_umi_reads,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$pct_usable_umi_reads,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$pct_usable_umi_reads~dati_paper_NGS$run, main="Reads that contain usable UMI information by run",ylab="Reads that contain usable UMI information (%)", xlab="RUN",axes=F,ylim=c(99,100))
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$pct_usable_umi_reads~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$pct_usable_umi_reads, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$pct_usable_umi_reads, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$pct_usable_umi_reads, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$pct_usable_umi_reads,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))


#Median.test(dati_paper_NGS$pct_usable_umi_reads,factor(dati_paper_NGS$run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



library(coin)
wilcox_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_usable_umi_reads[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_usable_umi_reads)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



################################
#dati_paper_NGS$pct_pf_uq_reads#
################################

summary(dati_paper_NGS$pct_pf_uq_reads); table(dati_paper_NGS$pct_pf_uq_reads>=55,useNA="ifany") 
sd(dati_paper_NGS$pct_pf_uq_reads,na.rm=T)
tapply(dati_paper_NGS$pct_pf_uq_reads,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$pct_pf_uq_reads,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$pct_pf_uq_reads,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$pct_pf_uq_reads,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$pct_pf_uq_reads~dati_paper_NGS$run, main="Unique reads passing filter by run",
ylab="Unique reads passing filter(%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$pct_pf_uq_reads~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$pct_pf_uq_reads, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$pct_pf_uq_reads, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$pct_pf_uq_reads, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$pct_pf_uq_reads,dati_paper_NGS$run,function(x) table(x>=55,useNA="ifany"))

#Median.test(dati_paper_NGS$pct_pf_uq_reads,factor(dati_paper_NGS$run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



library(coin)
wilcox_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_pf_uq_reads[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_pf_uq_reads)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



################################
#dati_paper_NGS$pct_target_100x#
################################

summary(dati_paper_NGS$pct_target_100x); #table(dati_paper_NGS$pct_target_100x>...,useNA="ifany") 
sd(dati_paper_NGS$pct_target_100x,na.rm=T)
tapply(dati_paper_NGS$pct_target_100x,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$pct_target_100x,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$pct_target_100x,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$pct_target_100x,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$pct_target_100x~dati_paper_NGS$run, main="Target 100x by run",
ylab="Target 100x (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$pct_target_100x~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$pct_target_100x, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$pct_target_100x, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$pct_target_100x, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$pct_target_100x,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))


#Median.test(dati_paper_NGS$pct_target_100x,factor(dati_paper_NGS$run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



library(coin)
wilcox_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_target_100x[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_100x)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



################################
#dati_paper_NGS$pct_target_250x#
################################

summary(dati_paper_NGS$pct_target_250x); #table(dati_paper_NGS$pct_target_250x>...,useNA="ifany") 
sd(dati_paper_NGS$pct_target_250x,na.rm=T)
tapply(dati_paper_NGS$pct_target_250x,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$pct_target_250x,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$pct_target_250x,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$pct_target_250x,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$pct_target_250x~dati_paper_NGS$run, main="Target 250x by run",
ylab="Target 250x (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$pct_target_250x~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$pct_target_250x, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$pct_target_250x, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$pct_target_250x, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$pct_target_250x,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))


#Median.test(dati_paper_NGS$pct_target_250x,factor(dati_paper_NGS$run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



library(coin)
wilcox_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_target_250x[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$pct_target_250x)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr




#######################################
#dati_paper_NGS$rna_pct_chimeric_reads#
#######################################

summary(dati_paper_NGS$rna_pct_chimeric_reads); #table(dati_paper_NGS$rna_pct_chimeric_reads>...,useNA="ifany") 
sd(dati_paper_NGS$rna_pct_chimeric_reads,na.rm=T)
tapply(dati_paper_NGS$rna_pct_chimeric_reads,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$rna_pct_chimeric_reads,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$rna_pct_chimeric_reads,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$rna_pct_chimeric_reads,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$rna_pct_chimeric_reads~dati_paper_NGS$run, main="Chimeric reads by run",
ylab="Chimeric reads (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper_NGS$rna_pct_chimeric_reads~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$rna_pct_chimeric_reads, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$rna_pct_chimeric_reads, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$rna_pct_chimeric_reads, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$rna_pct_chimeric_reads,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))


#Median.test(dati_paper_NGS$rna_pct_chimeric_reads,factor(dati_paper_NGS$run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]~
factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]),
data = dati_paper_NGS, distribution = "exact")

###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



library(coin)
wilcox_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==2)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==3)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==1 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==3)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==2 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==4)]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==3 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$rna_pct_chimeric_reads[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]~factor(dati_paper_NGS$run[is.na(dati_paper_NGS$rna_pct_chimeric_reads)==F & (dati_paper_NGS$run==4 | dati_paper_NGS$run==5)]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



