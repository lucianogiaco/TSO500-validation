####################################################
#########PAPER VALIDAZIONE SU DESCRITTIVE###########
####################################################


#install.packages("vioplot")
#install.packages("coin")

#library(vioplot)
#library(ggplot2)
library(coin)
library(fmsb)

citation(package = "coin", lib.loc = NULL)
citation(package = "fmsb", lib.loc = NULL)

#giustificazione per l'utilizzo del median test (Brown-Mood median test)
#https://v8doc.sas.com/sashtml/stat/chap47/sect17.htm
?wilcox_test (link su exact del pacchetto coin)


###WET METRICS###

dati_paper<-read.csv2(file.choose(),sep=";",h=T,dec=",",na.strings="") #importa FILE VALIDAZIONE_25-3-22_rev_DG_overall_senza 4 run RNA_250322.csv
head(dati_paper)
dim(dati_paper) # 130  19
colnames(dati_paper)
summary(dati_paper)
table(dati_paper$RUN[dati_paper$Source=="DNA"],useNA = "ifany")
table(dati_paper$RUN[dati_paper$Source=="RNA"],useNA = "ifany")

###FARE L'ANALISI IN MANIERA STRATIFICATA: DNA SAMPLES E RNA SAMPLES###

table(dati_paper$Source,useNA = "ifany")
sum(table(dati_paper$Source,useNA = "ifany"))

length(unique(dati_paper$ID[dati_paper$Source=="DNA"])) #71
length(unique(dati_paper$ID[dati_paper$Source=="RNA"])) #52


###STATISTICHE DESCRITTIVE PER LE WET METRICS###

#DNA
par(mfrow=c(2,2))
#DNA..ng.ul....3.5..12uL.
summary(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"]); table(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"]>3.5,useNA="ifany") 
sd(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"],na.rm=T)
tapply(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],summary)
tapply(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],length)
tapply(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) length(x[is.na(x)==T]))
tapply(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) table(x>3.5,useNA="ifany"))

boxplot(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"], main="DNA quantification by run",ylab="DNA quantification (ng/ul)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"])
pairwise.wilcox.test(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$DNA..ng.ul....3.5..12uL.[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "BH")

library(coin)
fwt<-function(x,fac,tipo,type,i,j){
zbau<-wilcox_test(x[(fac==i | fac==j) & tipo==type]~factor(fac[(fac==i | fac==j) & tipo==type]),
distribution = "exact", conf.int = TRUE)
return(pvalue(zbau))
}

fwt(x=dati_paper$DNA..ng.ul....3.5..12uL.,fac=dati_paper$RUN,tipo=dati_paper$Source,type="DNA",i=1,j=2)
fwt(x=dati_paper$DNA..ng.ul....3.5..12uL.,fac=dati_paper$RUN,tipo=dati_paper$Source,type="DNA",i=1,j=3,tipo=dati_paper$Source)
fwt(x=dati_paper$DNA..ng.ul....3.5..12uL.,fac=dati_paper$RUN,tipo=dati_paper$Source,type="DNA",i=1,j=4,tipo=dati_paper$Source)
fwt(x=dati_paper$DNA..ng.ul....3.5..12uL.,fac=dati_paper$RUN,tipo=dati_paper$Source,type="DNA",i=1,j=5,tipo=dati_paper$Source)


#dati_paper$A.260.280
summary(dati_paper$A.260.280[dati_paper$Source=="DNA"]); table(dati_paper$A.260.280[dati_paper$Source=="DNA"]>2,useNA="ifany")
sd(dati_paper$A.260.280[dati_paper$Source=="DNA"],na.rm=T) 
tapply(dati_paper$A.260.280[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],summary)
tapply(dati_paper$A.260.280[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$A.260.280[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],length)
tapply(dati_paper$A.260.280[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$A.260.280[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"], main="DNA Protein Nanodrop Purity Qualification by run",ylab="Protein Nanodrop Purity Qualification (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=29)"))
axis(2)
kruskal.test(dati_paper$A.260.280[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"])
pairwise.wilcox.test(dati_paper$A.260.280[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$A.260.280[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$A.260.280[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "BH")
tapply(dati_paper$A.260.280[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) table(x>2,useNA="ifany"))


#dati_paper$A.260.230
summary(dati_paper$A.260.230[dati_paper$Source=="DNA"]); table(dati_paper$A.260.230[dati_paper$Source=="DNA"]>2,useNA="ifany")
sd(dati_paper$A.260.230[dati_paper$Source=="DNA"],na.rm=T)  
tapply(dati_paper$A.260.230[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],summary)
tapply(dati_paper$A.260.230[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$A.260.230[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],length)
tapply(dati_paper$A.260.230[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$A.260.230[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"], main="Carbohydrate Nanodrop Purity Qualification by run",ylab="Carbohydrate Nanodrop Purity Qualification (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper$A.260.230[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"])
pairwise.wilcox.test(dati_paper$A.260.230[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$A.260.230[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$A.260.230[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "BH")
tapply(dati_paper$A.260.230[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) table(x>2,useNA="ifany"))



#dati_paper$Delta.Ct...5.
summary(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"]); table(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"]<5,useNA="ifany")
sd(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"],na.rm=T)   
tapply(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],summary)
tapply(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],length)
tapply(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"], main="DNA Qualification by RT-PCR by run",ylab="DNA Qualification by RT-PCR (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"])
pairwise.wilcox.test(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "BH")
tapply(dati_paper$Delta.Ct...5.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) table(x<5,useNA="ifany"))


par(mfrow=c(1,3))
library(dplyr)
#dati_paper$QC.post.Fragmentation.150.300bp
summary(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"]); table(between(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"],150,300),useNA="ifany")
sd(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"],na.rm=T) 
tapply(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],summary)
tapply(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],length)
tapply(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"], main="DNA Fragmentation by sonication by run",ylab="DNA Fragmentation by sonication (bp)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"])
pairwise.wilcox.test(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "BH")
tapply(dati_paper$QC.post.Fragmentation.150.300bp[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) table(between(x,150,300)))


#IBRIDIZZAZZIONE DNA
#dati_paper$Pre.Hyb.ng.ul....30.
summary(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"]); table(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"]>20,useNA="ifany")
sd(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"],na.rm=T) 
tapply(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],summary)
tapply(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],length)
tapply(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"], main="DNA pre-hybridization measurement by run",ylab="DNA pre-hybridization measurement (ng/ul)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"])
pairwise.wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "BH")
tapply(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) table(x>20,useNA="ifany"))


#dati_paper$Post.ng.ul....3.
summary(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"]); table(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"]>3,useNA="ifany")
sd(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"],na.rm=T) 
tapply(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],summary)
tapply(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],length)
tapply(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"], main="DNA post-hybridization measurement by run",ylab="DNA post-hybridization measurement (ng/ul)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"]~dati_paper$RUN[dati_paper$Source=="DNA"])
pairwise.wilcox.test(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"], dati_paper$RUN[dati_paper$Source=="DNA"], p.adj= "BH")
tapply(dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"],dati_paper$RUN[dati_paper$Source=="DNA"],function(x) table(x>3,useNA="ifany"))


wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA"],dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA"],paired=T)
wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA" & dati_paper$RUN==1],dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA" & dati_paper$RUN==1],paired=T)
wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA" & dati_paper$RUN==2],dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA" & dati_paper$RUN==2],paired=T)
wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA" & dati_paper$RUN==3],dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA" & dati_paper$RUN==3],paired=T)
wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA" & dati_paper$RUN==4],dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA" & dati_paper$RUN==4],paired=T)
wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="DNA" & dati_paper$RUN==5],dati_paper$Post.ng.ul....3.[dati_paper$Source=="DNA" & dati_paper$RUN==5],paired=T)




library(dplyr)
#dati_paper$QC.post.Hybridization..250bp
summary(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3]); table(between(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3],250,300),useNA="ifany")
sd(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3],na.rm=T) 
tapply(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3],dati_paper$RUN[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3],summary)
tapply(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3],dati_paper$RUN[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3],function(x) sd(x,na.rm=T))
tapply(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3],dati_paper$RUN[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3],length)
tapply(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3],dati_paper$RUN[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3]~dati_paper$RUN[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3], main="QC post Hybridization by run",ylab="QC post Hybridization (bp)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=0)","2 (n=0)","3 (n=0)","4 (n=0)","5 (n=2)"))
axis(2)
kruskal.test(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3]~dati_paper$RUN[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3])
pairwise.wilcox.test(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3], dati_paper$RUN[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3], p.adj= "none")
pairwise.wilcox.test(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3], dati_paper$RUN[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3], dati_paper$RUN[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3], p.adj= "BH")
tapply(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3],dati_paper$RUN[dati_paper$Source=="DNA" & dati_paper$Post.ng.ul....3.<=3],function(x) table(between(x,250,300)))




#RNA
par(mfrow=c(2,3))
#dati_paper$RNA..ng.ul...10ng.ul..8.5ul.
summary(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"]); table(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"]>10.5,useNA="ifany") 
sd(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"],na.rm=T)
tapply(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],summary)
tapply(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],length)
tapply(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"]~dati_paper$RUN[dati_paper$Source=="RNA"], main="RNA quantification by run",ylab="RNA Qualification (ng/ul)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"]~dati_paper$RUN[dati_paper$Source=="RNA"])
pairwise.wilcox.test(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "BH")
tapply(dati_paper$RNA..ng.ul...10ng.ul..8.5ul.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) table(x>10.5,useNA="ifany"))


#dati_paper$RNA.A260.280
summary(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"]); table(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"]>2,useNA="ifany")
sd(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"],na.rm=T) 
tapply(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],summary)
tapply(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],length)
tapply(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"]~dati_paper$RUN[dati_paper$Source=="RNA"], main="RNA Protein Purity Qualification by run",ylab="Protein Nanodrop Purity Qualification (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"]~dati_paper$RUN[dati_paper$Source=="RNA"])
pairwise.wilcox.test(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "BH")
tapply(dati_paper$RNA.A260.280[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) table(x>2,useNA="ifany"))


#dati_paper$RNA.A260.230
summary(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"]); table(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"]>2,useNA="ifany")
sd(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"],na.rm=T) 
tapply(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],summary)
tapply(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],length)
tapply(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"]~dati_paper$RUN[dati_paper$Source=="RNA"], main="RNA Carbohydrate Purity Qualification by run",ylab="Carbohydrate Purity Qualification (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"]~dati_paper$RUN[dati_paper$Source=="RNA"])
pairwise.wilcox.test(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "BH")
tapply(dati_paper$RNA.A260.230[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) table(x>2,useNA="ifany"))


#dati_paper$X.DV200...20..
summary(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"]); table(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"]>20,useNA="ifany")
sd(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"],na.rm=T) 
tapply(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],summary)
tapply(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],length)
tapply(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"]~dati_paper$RUN[dati_paper$Source=="RNA"], main="RNA Peak signal measurement by run",ylab="RNA Peak signal measurement (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"]~dati_paper$RUN[dati_paper$Source=="RNA"])
pairwise.wilcox.test(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "BH")
tapply(dati_paper$X.DV200...20..[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) table(x>20,useNA="ifany"))


#IBRIDIZZAZZIONE RNA
#dati_paper$Pre.Hyb.ng.ul....30.
summary(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"]); table(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"]>20,useNA="ifany")
sd(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"],na.rm=T) 
tapply(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],summary)
tapply(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],length)
tapply(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"]~dati_paper$RUN[dati_paper$Source=="RNA"], main="RNA pre-hybridization measurement by run",ylab="RNA pre-hybridization measurement (ng/ul)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"]~dati_paper$RUN[dati_paper$Source=="RNA"])
pairwise.wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "BH")
tapply(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) table(x>20,useNA="ifany"))

library(coin)
wilcox_test(dati_paper$Pre.Hyb.ng.ul....30.~dati_paper$Source,distribution = "exact", conf.int = TRUE)


#dati_paper$Post.ng.ul....3.
summary(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"]); table(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"]>3,useNA="ifany")
sd(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"],na.rm=T) 
tapply(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],summary)
tapply(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) sd(x,na.rm=T))
tapply(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],length)
tapply(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"]~dati_paper$RUN[dati_paper$Source=="RNA"], main="RNA post-hybridization measurement by run",ylab="RNA post-hybridization measurement (ng/ul)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"]~dati_paper$RUN[dati_paper$Source=="RNA"])
pairwise.wilcox.test(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "none")
pairwise.wilcox.test(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"], dati_paper$RUN[dati_paper$Source=="RNA"], p.adj= "BH")
tapply(dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"],dati_paper$RUN[dati_paper$Source=="RNA"],function(x) table(x>3,useNA="ifany"))

library(coin)
wilcox_test(dati_paper$Post.ng.ul....3.~dati_paper$Source,distribution = "exact", conf.int = TRUE)

wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA"],dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA"],paired=T)
wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA" & dati_paper$RUN==1],dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA" & dati_paper$RUN==1],paired=T)
wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA" & dati_paper$RUN==2],dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA" & dati_paper$RUN==2],paired=T)
wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA" & dati_paper$RUN==3],dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA" & dati_paper$RUN==3],paired=T)
wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA" & dati_paper$RUN==4],dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA" & dati_paper$RUN==4],paired=T)
wilcox.test(dati_paper$Pre.Hyb.ng.ul....30.[dati_paper$Source=="RNA" & dati_paper$RUN==5],dati_paper$Post.ng.ul....3.[dati_paper$Source=="RNA" & dati_paper$RUN==5],paired=T)


library(dplyr)
#dati_paper$QC.post.Hybridization..250bp
summary(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3]); table(between(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3],250,300),useNA="ifany")
sd(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3],na.rm=T) 
tapply(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3],dati_paper$RUN[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3],summary)
tapply(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3],dati_paper$RUN[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3],function(x) sd(x,na.rm=T))
tapply(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3],dati_paper$RUN[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3],length)
tapply(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3],dati_paper$RUN[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3],function(x) length(x[is.na(x)==T]))
boxplot(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3]~dati_paper$RUN[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3], main="QC post Hybridization by run",ylab="QC post Hybridization (bp)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=2)","2 (n=0)","3 (n=4)","4 (n=1)","5 (n=0)"))
axis(2)
kruskal.test(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3]~dati_paper$RUN[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3])
pairwise.wilcox.test(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3], dati_paper$RUN[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3], p.adj= "none")
pairwise.wilcox.test(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3], dati_paper$RUN[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3], p.adj= "fdr")
pairwise.wilcox.test(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3], dati_paper$RUN[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3], p.adj= "BH")
tapply(dati_paper$QC.post.Hybridization..250bp[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3],dati_paper$RUN[dati_paper$Source=="RNA" & dati_paper$Post.ng.ul....3.<=3],function(x) table(between(x,250,300),useNA="ifany"))






###SEQUENCING METRICS###

dati_paper_NGS<-read.csv2(file.choose(),sep="\t",h=T,dec=".",na.strings="NA") #importa tab2.tab
head(dati_paper_NGS)
dim(dati_paper_NGS) #112 47
colnames(dati_paper_NGS)
summary(dati_paper_NGS)
table(dati_paper_NGS$run) #numero di run: 5

dati_paper_NGS[dati_paper_NGS$run=="210715_A01423_0008_AH35CWDRXY",]
dati_paper_NGS[dati_paper_NGS$run=="210729_A01423_0009_AH33WGDRXY",]
dati_paper_NGS[dati_paper_NGS$run=="211022_A01423_0010_AHGYFYDRXY",]
dati_paper_NGS[dati_paper_NGS$run=="211111_A01423_0011_AHH2Y2DRXY",]
dati_paper_NGS[dati_paper_NGS$run=="211122_A01423_0012_AH2YWCDRXY",]

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

dim(datiDNAsample) #112 18 
summary(datiDNAsample)

apply(datiDNAsample,2,function(x) sd(x,na.rm=T))



##########Kruskal-Wallis and Mann-Whitney###########

###Soglie per linee guida: metrics Output dal file
###trusight-oncology-500-local-app-v2.2-user-guide-1000000137777-01.pdf

options(scipen=999999999)


par(mfrow=c(2,3))
#dati_paper_NGS$total_pf_reads
summary(dati_paper_NGS$total_pf_reads/1000000); #table(dati_paper_NGS$total_pf_reads>...,useNA="ifany") 
sd(dati_paper_NGS$total_pf_reads/1000000,na.rm=T)
tapply(dati_paper_NGS$total_pf_reads/1000000,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$total_pf_reads/1000000,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$total_pf_reads/1000000,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$total_pf_reads/1000000,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$total_pf_reads/1000000~dati_paper_NGS$run, main="Total number of reads passing filter by run",
ylab="Total number of reads passing filter (1000000u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$total_pf_reads/1000000~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$total_pf_reads/1000000, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$total_pf_reads/1000000, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$total_pf_reads/1000000, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$total_pf_reads/1000000,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))



#dati_paper_NGS$median_target_coverage
summary(dati_paper_NGS$median_target_coverage); #table(dati_paper_NGS$median_target_coverage>...,useNA="ifany") 
sd(dati_paper_NGS$median_target_coverage,na.rm=T)
tapply(dati_paper_NGS$median_target_coverage,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$median_target_coverage,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$median_target_coverage,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$median_target_coverage,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$median_target_coverage~dati_paper_NGS$run, main="Median target coverage by run",ylab="Median target coverage (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$median_target_coverage~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$median_target_coverage, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$median_target_coverage, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$median_target_coverage, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$median_target_coverage,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))


#dati_paper_NGS$pct_chimeric_reads
summary(dati_paper_NGS$pct_chimeric_reads); #table(dati_paper_NGS$pct_chimeric_reads>...,useNA="ifany") 
sd(dati_paper_NGS$pct_chimeric_reads,na.rm=T)
tapply(dati_paper_NGS$pct_chimeric_reads,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$pct_chimeric_reads,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$pct_chimeric_reads,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$pct_chimeric_reads,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$pct_chimeric_reads~dati_paper_NGS$run, main="Chimeric reads by run",ylab="Chimeric reads (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$pct_chimeric_reads~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$pct_chimeric_reads, dati_paper$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$pct_chimeric_reads, dati_paper$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$pct_chimeric_reads, dati_paper$run, p.adj= "BH")
#tapply(dati_paper_NGS$pct_chimeric_reads,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))


#dati_paper_NGS$pct_exon_100x
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

#Median.test(dati_cov_SNV$Q2x500,factor(dati_cov_SNV$Run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==1 | dati_paper_NGS$run==2]~
factor(dati_paper_NGS$run[dati_paper_NGS$run==1 | dati_paper_NGS$run==2]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==1 | dati_paper_NGS$run==3]~
factor(dati_paper_NGS$run[dati_paper_NGS$run==1 | dati_paper_NGS$run==3]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==1 | dati_paper_NGS$run==4]~
factor(dati_paper_NGS$run[dati_paper_NGS$run==1 | dati_paper_NGS$run==4]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==1 | dati_paper_NGS$run==5]~
factor(dati_paper_NGS$run[dati_paper_NGS$run==1 | dati_paper_NGS$run==5]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==2 | dati_paper_NGS$run==3]~
factor(dati_paper_NGS$run[dati_paper_NGS$run==2 | dati_paper_NGS$run==3]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==2 | dati_paper_NGS$run==4]~
factor(dati_paper_NGS$run[dati_paper_NGS$run==2 | dati_paper_NGS$run==4]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==2 | dati_paper_NGS$run==5]~
factor(dati_paper_NGS[dati_paper_NGS$Run==2 | dati_paper_NGS$Run==5]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==3 | dati_paper_NGS$run==4]~
factor(dati_paper_NGS$Run[dati_paper_NGS$run==3 | dati_paper_NGS$run==4]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==3 | dati_paper_NGS$run==5]~
factor(dati_paper_NGS$run[dati_paper_NGS$run==3 | dati_paper_NGS$run==5]),
data = dati_paper_NGS, distribution = "exact")

median_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==4 | dati_paper_NGS$run==5]~
factor(dati_paper_NGS$run[dati_paper_NGS$run==4 | dati_paper_NGS$run==5]),
data = dati_paper_NGS, distribution = "exact")

###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



library(coin)
wilcox_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==1 | dati_paper_NGS$run==2]~factor(dati_paper_NGS$run[dati_paper_NGS$run==1 | dati_paper_NGS$run==2]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==1 | dati_paper_NGS$run==3]~factor(dati_paper_NGS$run[dati_paper_NGS$run==1 | dati_paper_NGS$run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==1 | dati_paper_NGS$run==4]~factor(dati_paper_NGS$run[dati_paper_NGS$run==1 | dati_paper_NGS$run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==1 | dati_paper_NGS$run==5]~factor(dati_paper_NGS$run[dati_paper_NGS$run==1 | dati_paper_NGS$run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==2 | dati_paper_NGS$run==3]~factor(dati_paper_NGS$run[dati_paper_NGS$run==2 | dati_paper_NGS$run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==2 | dati_paper_NGS$run==4]~factor(dati_paper_NGS$run[dati_paper_NGS$run==2 | dati_paper_NGS$run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==2 | dati_paper_NGS$run==5]~factor(dati_paper_NGS$run[dati_paper_NGS$run==2 | dati_paper_NGS$run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==3 | dati_paper_NGS$run==4]~factor(dati_paper_NGS$run[dati_paper_NGS$run==3 | dati_paper_NGS$run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==3 | dati_paper_NGS$run==5]~factor(dati_paper_NGS$run[dati_paper_NGS$run==3 | dati_paper_NGS$run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_paper_NGS$pct_exon_100x[dati_paper_NGS$run==4 | dati_paper_NGS$run==5]~factor(dati_paper_NGS$run[dati_paper_NGS$run==4 | dati_paper_NGS$run==5]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr





#dati_paper_NGS$pct_read_enrichment
summary(dati_paper_NGS$pct_read_enrichment); #table(dati_paper_NGS$pct_read_enrichment>...,useNA="ifany") 
sd(dati_paper_NGS$pct_read_enrichment,na.rm=T)
tapply(dati_paper_NGS$pct_read_enrichment,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$pct_read_enrichment,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$pct_read_enrichment,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$pct_read_enrichment,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$pct_read_enrichment~dati_paper_NGS$run, main="Read enrichment by run",ylab="Read enrichment (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$pct_read_enrichment~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$pct_read_enrichment, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$pct_read_enrichment, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$pct_read_enrichment, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$pct_read_enrichment,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))



#dati_paper_NGS$pct_usable_umi_reads
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


par(mfrow=c(2,3))
#dati_paper_NGS$mean_target_coverage
summary(dati_paper_NGS$mean_target_coverage); #table(dati_paper_NGS$mean_target_coverage>...,useNA="ifany") 
sd(dati_paper_NGS$mean_target_coverage,na.rm=T)
tapply(dati_paper_NGS$mean_target_coverage,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$mean_target_coverage,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$mean_target_coverage,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$mean_target_coverage,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$mean_target_coverage~dati_paper_NGS$run, main="Mean target coverage by run",ylab="Mean target coverage (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$mean_target_coverage~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$mean_target_coverage, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$mean_target_coverage, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$mean_target_coverage, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$mean_target_coverage,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))



#dati_paper_NGS$pct_aligned_reads
summary(dati_paper_NGS$pct_aligned_reads); #table(dati_paper_NGS$pct_aligned_reads>...,useNA="ifany") 
sd(dati_paper_NGS$pct_aligned_reads,na.rm=T)
tapply(dati_paper_NGS$pct_aligned_reads,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$pct_aligned_reads,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$pct_aligned_reads,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$pct_aligned_reads,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$pct_aligned_reads~dati_paper_NGS$run, main="Aligned reads by run",ylab="Aligned reads (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$pct_aligned_reads~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$pct_aligned_reads, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$pct_aligned_reads, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$pct_aligned_reads, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$pct_aligned_reads,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))



#dati_paper_NGS$pct_contamination_est
summary(dati_paper_NGS$pct_contamination_est); #table(dati_paper_NGS$pct_contamination_est>...,useNA="ifany") 
sd(dati_paper_NGS$pct_contamination_est,na.rm=T)
tapply(dati_paper_NGS$pct_contamination_est,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$pct_contamination_est,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$pct_contamination_est,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$pct_contamination_est,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$pct_contamination_est~dati_paper_NGS$run, main="Percent of foreign DNA estimated to be in the sample by run",
ylab="Percent of foreign DNA estimated to be in the sample (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$pct_contamination_est~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$pct_contamination_est, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$pct_contamination_est, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$pct_contamination_est, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$pct_contamination_est,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))



#dati_paper_NGS$pct_pf_uq_reads
table(dati_paper_NGS$PCT_PF_READS_guidelines)
tapply(dati_paper_NGS$PCT_PF_READS_guidelines,dati_paper_NGS$run, function(x) table(x,useNA="ifany"))

data.frame(dati_paper_NGS$PCT_PF_READS_perc,dati_paper_NGS$pct_pf_uq_reads)

#dati_paper_NGS$pct_pf_uq_reads
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
pairwise.wilcox.test(dati_paper_NGS$pct_pf_uq_reads, dati_paper$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$pct_pf_uq_reads, dati_paper$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$pct_pf_uq_reads, dati_paper$run, p.adj= "BH")
#tapply(dati_paper_NGS$pct_pf_uq_reads,dati_paper_NGS$run,function(x) table(x>=55,useNA="ifany"))



#dati_paper_NGS$pct_target_04x_mean
summary(dati_paper_NGS$pct_target_04x_mean); #table(dati_paper_NGS$pct_target_04x_mean>...,useNA="ifany") 
sd(dati_paper_NGS$pct_target_04x_mean,na.rm=T)
tapply(dati_paper_NGS$pct_target_04x_mean,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$pct_target_04x_mean,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$pct_target_04x_mean,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$pct_target_04x_mean,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$pct_target_04x_mean~dati_paper_NGS$run, main="Target 04x mean by run",
ylab="Target 04x mean (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$pct_target_04x_mean~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$pct_target_04x_mean, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$pct_target_04x_mean, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$pct_target_04x_mean, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$pct_target_04x_mean,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))



#dati_paper_NGS$pct_target_100x
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


par(mfrow=c(2,3))
#dati_paper_NGS$pct_target_250x
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


#dati_paper_NGS$contamination_score
summary(dati_paper_NGS$contamination_score); table(dati_paper_NGS$contamination_score<=3106,useNA="ifany") #53
summary(dati_paper_NGS$contamination_score[dati_paper_NGS$contamination_p_value<=0.049 & dati_paper_NGS$contamination_score>3106]) #15
length(na.omit(dati_paper_NGS$contamination_score[dati_paper_NGS$contamination_p_value<=0.049 & dati_paper_NGS$contamination_score>3106])) #15
#data.frame(dati_paper_NGS$contamination_score[dati_paper_NGS$contamination_p_value<=0.049 & dati_paper_NGS$contamination_score>3106],
#dati_paper_NGS$Contamination[dati_paper_NGS$contamination_p_value<=0.049 & dati_paper_NGS$contamination_score>3106])

sd(dati_paper_NGS$contamination_score,na.rm=T)
tapply(dati_paper_NGS$contamination_score,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$contamination_score,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$contamination_score,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$contamination_score,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$contamination_score~dati_paper_NGS$run, main="Contamination score by run",
ylab="Contamination score (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$contamination_score~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$contamination_score, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$contamination_score, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$contamination_score, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$contamination_score,dati_paper_NGS$run,function(x) table(x<=3106,useNA="ifany"))

#dati_paper_NGS$Contamination
tapply(dati_paper_NGS$Contamination,dati_paper_NGS$run, function(x) table(x,useNA="ifany"))
tapply(dati_paper_NGS$contamination_score,dati_paper_NGS$Contamination,summary)
tapply(dati_paper_NGS$contamination_score,dati_paper_NGS$Contamination,length)
table(dati_paper_NGS$contamination_score[dati_paper_NGS$Contamination==1]<=3106,useNA="ifany") #15
table(dati_paper_NGS$contamination_score[dati_paper_NGS$Contamination==1]>3106,useNA="ifany") #15

table(is.na(dati_paper_NGS$contamination_p_value[dati_paper_NGS$Contamination==1 & dati_paper_NGS$contamination_score>3106])==F)
summary(dati_paper_NGS$contamination_p_value[dati_paper_NGS$Contamination==1 & dati_paper_NGS$contamination_score>3106]) #15

tapply(dati_paper_NGS$contamination_p_value,dati_paper_NGS$Contamination,summary)
tapply(dati_paper_NGS$contamination_p_value,dati_paper_NGS$Contamination,length)


na.omit(data.frame(dati_paper_NGS$contamination_score[dati_paper_NGS$Contamination==0],
dati_paper_NGS$contamination_p_value[dati_paper_NGS$Contamination==0]))

na.omit(data.frame(dati_paper_NGS$contamination_score[dati_paper_NGS$Contamination==1],
dati_paper_NGS$contamination_p_value[dati_paper_NGS$Contamination==1]))



#dati_paper_NGS$contamination_p_value
summary(dati_paper_NGS$contamination_p_value) 
sd(dati_paper_NGS$contamination_p_value,na.rm=T)
tapply(dati_paper_NGS$contamination_p_value,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$contamination_p_value,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$contamination_p_value,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$contamination_p_value,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$contamination_p_value~dati_paper_NGS$run, main="Contamination P value by run",
ylab="Contamination P value (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$contamination_p_value~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$contamination_p_value, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$contamination_p_value, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$contamination_p_value, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$contamination_p_value,dati_paper_NGS$run,function(x) table(x<=0.049 & dati_paper_NGS$contamination_score>3106,useNA="ifany"))




#dati_paper_NGS$coverage_mad_guidelines
tapply(dati_paper_NGS$coverage_mad_guidelines,dati_paper_NGS$run, function(x) table(x,useNA="ifany"))

#dati_paper_NGS$coverage_mad_count
summary(dati_paper_NGS$coverage_mad_count); table(dati_paper_NGS$coverage_mad_count<=0.21,useNA="ifany") 
sd(dati_paper_NGS$coverage_mad_count,na.rm=T)
tapply(dati_paper_NGS$coverage_mad_count,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$coverage_mad_count,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$coverage_mad_count,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$coverage_mad_count,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$coverage_mad_count~dati_paper_NGS$run, main="Coverage MAD by run",
ylab="Coverage MAD (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$coverage_mad_count~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$coverage_mad_count, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$coverage_mad_count, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$coverage_mad_count, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$coverage_mad_count,dati_paper_NGS$run,function(x) table(x<=0.21,useNA="ifany"))



#dati_paper_NGS$median_bin_count_cnv_target_guidelines
tapply(dati_paper_NGS$median_bin_count_cnv_target_guidelines,dati_paper_NGS$run, function(x) table(x,useNA="ifany"))

#dati_paper_NGS$median_bin_count_cnv_target_count
summary(dati_paper_NGS$median_bin_count_cnv_target_count); table(dati_paper_NGS$median_bin_count_cnv_target_count>=1,useNA="ifany") 
sd(dati_paper_NGS$median_bin_count_cnv_target_count,na.rm=T)
tapply(dati_paper_NGS$median_bin_count_cnv_target_count,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$median_bin_count_cnv_target_count,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$median_bin_count_cnv_target_count,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$median_bin_count_cnv_target_count,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$median_bin_count_cnv_target_count~dati_paper_NGS$run, main="Median bin CNV target by run",
ylab="Median bin CNV target (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$median_bin_count_cnv_target_count~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$median_bin_count_cnv_target_count, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$median_bin_count_cnv_target_count, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$median_bin_count_cnv_target_count, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$median_bin_count_cnv_target_count,dati_paper_NGS$run,function(x) table(x>=1,useNA="ifany"))



#dati_paper_NGS$median_insert_size_guidelines
tapply(dati_paper_NGS$median_insert_size_guidelines,dati_paper_NGS$run, function(x) table(x,useNA="ifany"))

#dati_paper_NGS$median_insert_size
summary(dati_paper_NGS$median_insert_size); table(dati_paper_NGS$median_insert_size>=70,useNA="ifany") 
sd(dati_paper_NGS$median_insert_size,na.rm=T)
tapply(dati_paper_NGS$median_insert_size,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$median_insert_size,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$median_insert_size,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$median_insert_size,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$median_insert_size~dati_paper_NGS$run, main="DNA median insert size by run",
ylab="DNA median_insert_size (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$median_insert_size~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$median_insert_size, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$median_insert_size, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$median_insert_size, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$median_insert_size,dati_paper_NGS$run,function(x) table(x>=70,useNA="ifany"))


#dati_paper_NGS$median_exon_coverage_guidelines
tapply(dati_paper_NGS$median_exon_coverage_guidelines,dati_paper_NGS$run, function(x) table(x,useNA="ifany"))


par(mfrow=c(2,3))
#dati_paper_NGS$median_exon_coverage_count
summary(dati_paper_NGS$median_exon_coverage_count); table(dati_paper_NGS$median_exon_coverage_count>=150,useNA="ifany") 
sd(dati_paper_NGS$median_exon_coverage_count,na.rm=T)
tapply(dati_paper_NGS$median_exon_coverage_count,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$median_exon_coverage_count,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$median_exon_coverage_count,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$median_exon_coverage_count,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$median_exon_coverage_count~dati_paper_NGS$run, main="Median exon coverage by run",
ylab="Median exon coverage (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$median_exon_coverage_count~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$median_exon_coverage_count, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$median_exon_coverage_count, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$median_exon_coverage_count, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$median_exon_coverage_count,dati_paper_NGS$run,function(x) table(x>=150,useNA="ifany"))


#dati_paper_NGS$pct_exon_50x_guidelines
tapply(dati_paper_NGS$pct_exon_50x_guidelines,dati_paper_NGS$run, function(x) table(x,useNA="ifany"))

par(mfrow=c(1,2))
#dati_paper_NGS$pct_exon_50x_perc
summary(dati_paper_NGS$pct_exon_50x_perc); table(dati_paper_NGS$pct_exon_50x_perc>=90,useNA="ifany") 
sd(dati_paper_NGS$pct_exon_50x_perc,na.rm=T)
tapply(dati_paper_NGS$pct_exon_50x_perc,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$pct_exon_50x_perc,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$pct_exon_50x_perc,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$pct_exon_50x_perc,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$pct_exon_50x_perc~dati_paper_NGS$run, main="Percentage of exon 50x by run",
ylab="Exon 50x (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$pct_exon_50x_perc~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$pct_exon_50x_perc, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$pct_exon_50x_perc, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$pct_exon_50x_perc, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$pct_exon_50x_perc,dati_paper_NGS$run,function(x) table(x>=90,useNA="ifany"))


#dati_paper_NGS$unstable_msi_sites_guidelines
tapply(dati_paper_NGS$unstable_msi_sites_guidelines,dati_paper_NGS$run, function(x) table(x,useNA="ifany"))

#dati_paper_NGS$unstable_msi_sites
summary(dati_paper_NGS$unstable_msi_sites); table(dati_paper_NGS$unstable_msi_sites>=40,useNA="ifany") 
sd(dati_paper_NGS$unstable_msi_sites,na.rm=T)
tapply(dati_paper_NGS$unstable_msi_sites,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$unstable_msi_sites,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$unstable_msi_sites,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$unstable_msi_sites,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$unstable_msi_sites~dati_paper_NGS$run, main="Unstable (MSI) sites by run",
ylab="Unstable (MSI) sites (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
kruskal.test(dati_paper_NGS$unstable_msi_sites~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$unstable_msi_sites, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$unstable_msi_sites, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$unstable_msi_sites, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$unstable_msi_sites,dati_paper_NGS$run,function(x) table(x>=40,useNA="ifany"))



#########
###RNA###
#########


datiRNAsample<-data.frame(dati_paper_NGS$rna_total_pf_reads,
dati_paper_NGS$rna_scaled_median_gene_coverage,
dati_paper_NGS$rna_pct_chimeric_reads,
dati_paper_NGS$rna_pct_on_target_reads,
dati_paper_NGS$rna_median_cv__gene_500x_perc)

dim(datiRNAsample) #112 5
summary(datiRNAsample)


apply(datiRNAsample,2,function(x) sd(x,na.rm=T))



##########Kruskal-Wallis and Mann-Whitney###########

###Soglie per linee guida: metrics Output dal file
###trusight-oncology-500-local-app-v2.2-user-guide-1000000137777-01.pdf

options(scipen=999999999)


par(mfrow=c(2,2))
#dati_paper_NGS$rna_total_pf_reads
summary(dati_paper_NGS$rna_total_pf_reads/1000000); #table(dati_paper_NGS$rna_total_pf_reads>...,useNA="ifany") 
sd(dati_paper_NGS$rna_total_pf_reads/1000000,na.rm=T)
tapply(dati_paper_NGS$rna_total_pf_reads/1000000,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$rna_total_pf_reads/1000000,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$rna_total_pf_reads/1000000,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$rna_total_pf_reads/1000000,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$rna_total_pf_reads/1000000~dati_paper_NGS$run, main="Total number of reads passing filter by run",
ylab="Total number of reads passing filter (1000000u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper_NGS$rna_total_pf_reads/1000000~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$rna_total_pf_reads/1000000, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$rna_total_pf_reads/1000000, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$rna_total_pf_reads/1000000, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$rna_total_pf_reads/1000000,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))


#dati_paper_NGS$rna_scaled_median_gene_coverage
summary(dati_paper_NGS$rna_scaled_median_gene_coverage); #table(dati_paper_NGS$rna_scaled_median_gene_coverage>...,useNA="ifany") 
sd(dati_paper_NGS$rna_scaled_median_gene_coverage,na.rm=T)
tapply(dati_paper_NGS$rna_scaled_median_gene_coverage,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$rna_scaled_median_gene_coverage,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$rna_scaled_median_gene_coverage,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$rna_scaled_median_gene_coverage,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$rna_scaled_median_gene_coverage~dati_paper_NGS$run, main="Scaled median gene coverage by run",
ylab="Scaled median gene coverage (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper_NGS$rna_scaled_median_gene_coverage~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$rna_scaled_median_gene_coverage, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$rna_scaled_median_gene_coverage, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$rna_scaled_median_gene_coverage, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$rna_scaled_median_gene_coverage,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))


#dati_paper_NGS$rna_total_on_target_reads_guidelines
tapply(dati_paper_NGS$rna_total_on_target_reads_guidelines,dati_paper_NGS$run, function(x) table(x,useNA="ifany"))

#dati_paper_NGS$rna_total_on_target_reads_count
summary(dati_paper_NGS$rna_total_on_target_reads_count/1000000); table(dati_paper_NGS$rna_total_on_target_reads_count>=9000000,useNA="ifany") 
sd(dati_paper_NGS$rna_total_on_target_reads_count/1000000,na.rm=T)
tapply(dati_paper_NGS$rna_total_on_target_reads_count/1000000,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$rna_total_on_target_reads_count/1000000,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$rna_total_on_target_reads_count/1000000,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$rna_total_on_target_reads_count/1000000,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$rna_total_on_target_reads_count/1000000~dati_paper_NGS$run, 
main="The total number of reads that map to the target regions by run",
ylab="Total number of reads (1000000u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper_NGS$rna_total_on_target_reads_count/1000000~dati_paper_NGS$run)
#kruskal.test(dati_paper_NGS$rna_total_on_target_reads_count~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$rna_total_on_target_reads_count/1000000, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$rna_total_on_target_reads_count/1000000, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$rna_total_on_target_reads_count/1000000, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$rna_total_on_target_reads_count/1000000,dati_paper_NGS$run,function(x) table(x>=9,useNA="ifany"))



#dati_paper_NGS$rna_median_insert_size_guidelines
tapply(dati_paper_NGS$rna_median_insert_size_guidelines,dati_paper_NGS$run, function(x) table(x,useNA="ifany"))

#dati_paper_NGS$rna_median_insert_size_count
summary(dati_paper_NGS$rna_median_insert_size_count); table(dati_paper_NGS$rna_median_insert_size_count>=80,useNA="ifany") 
sd(dati_paper_NGS$rna_median_insert_size_count,na.rm=T)
tapply(dati_paper_NGS$rna_median_insert_size_count,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$rna_median_insert_size_count,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$rna_median_insert_size_count,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$rna_median_insert_size_count,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$rna_median_insert_size_count~dati_paper_NGS$run, main="RNA Median insert size by run",
ylab="RNA Median insert size (u)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper_NGS$rna_median_insert_size_count~dati_paper_NGS$run)
#kruskal.test(dati_paper_NGS$rna_median_insert_size_count~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$rna_median_insert_size_count, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$rna_median_insert_size_count, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$rna_median_insert_size_count, dati_paper_NGSr$run, p.adj= "BH")
#tapply(dati_paper_NGS$rna_median_insert_size_count,dati_paper_NGS$run,function(x) table(x>=80,useNA="ifany"))


par(mfrow=c(1,3))
#dati_paper_NGS$rna_pct_chimeric_reads
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



#dati_paper_NGS$rna_pct_on_target_reads
summary(dati_paper_NGS$rna_pct_on_target_reads); #table(dati_paper_NGS$rna_pct_on_target_reads>...,useNA="ifany") 
sd(dati_paper_NGS$rna_pct_on_target_reads,na.rm=T)
tapply(dati_paper_NGS$rna_pct_on_target_reads,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$rna_pct_on_target_reads,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$rna_pct_on_target_reads,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$rna_pct_on_target_reads,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$rna_pct_on_target_reads~dati_paper_NGS$run, main="Target reads by run",
ylab="Target reads (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper_NGS$rna_pct_on_target_reads~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$rna_pct_on_target_reads, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$rna_pct_on_target_reads, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$rna_pct_on_target_reads, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$rna_pct_on_target_reads,dati_paper_NGS$run,function(x) table(x>...,useNA="ifany"))



#dati_paper_NGS$rna_median_cv__gene_500x_guidelines
tapply(dati_paper_NGS$rna_median_cv__gene_500x_guidelines,dati_paper_NGS$run, function(x) table(x,useNA="ifany"))

#dati_paper_NGS$rna_median_cv__gene_500x_perc
summary(dati_paper_NGS$rna_median_cv__gene_500x_perc); table(dati_paper_NGS$rna_median_cv__gene_500x_perc<=93,useNA="ifany") 
sd(dati_paper_NGS$rna_median_cv__gene_500x_perc,na.rm=T)
tapply(dati_paper_NGS$rna_median_cv__gene_500x_perc,dati_paper_NGS$run,summary)
tapply(dati_paper_NGS$rna_median_cv__gene_500x_perc,dati_paper_NGS$run,function(x) sd(x,na.rm=T))
tapply(dati_paper_NGS$rna_median_cv__gene_500x_perc,dati_paper_NGS$run,length)
tapply(dati_paper_NGS$rna_median_cv__gene_500x_perc,dati_paper_NGS$run,function(x) length(x[is.na(x)==T]))
boxplot(dati_paper_NGS$rna_median_cv__gene_500x_perc~dati_paper_NGS$run, main="Median gene coverage 500x by run",
ylab="Median gene coverage 500x (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
kruskal.test(dati_paper_NGS$rna_median_cv__gene_500x_perc~dati_paper_NGS$run)
pairwise.wilcox.test(dati_paper_NGS$rna_median_cv__gene_500x_perc, dati_paper_NGS$run, p.adj= "none")
pairwise.wilcox.test(dati_paper_NGS$rna_median_cv__gene_500x_perc, dati_paper_NGS$run, p.adj= "fdr")
pairwise.wilcox.test(dati_paper_NGS$rna_median_cv__gene_500x_perc, dati_paper_NGS$run, p.adj= "BH")
#tapply(dati_paper_NGS$rna_median_cv__gene_500x_perc,dati_paper_NGS$run,function(x) table(x<=93,useNA="ifany"))



###Overall run metrics###

PF_READS_perc_overall <-c(86.0,81.4,89.0,86.2,86.7)
mean(PF_READS_perc_overall)
sd(PF_READS_perc_overall)
median(PF_READS_perc_overall)
min(PF_READS_perc_overall); max(PF_READS_perc_overall)


PCT_Q30_R1_overall <-c(93.9,92.9,94.4,93.7,93.9)
mean(PCT_Q30_R1_overall)
sd(PCT_Q30_R1_overall)
median(PCT_Q30_R1_overall)
min(PCT_Q30_R1_overall); max(PCT_Q30_R1_overall)


PCT_Q30_R2_overall <-c(92.0,91.3,94.8,93.5,93.8)
mean(PCT_Q30_R2_overall)
sd(PCT_Q30_R2_overall)
median(PCT_Q30_R2_overall)
min(PCT_Q30_R2_overall); max(PCT_Q30_R2_overall)


###RADAR PLOT###

#install.packages("fmsb")
library(fmsb)

citation(package = "fmsb", lib.loc = NULL)

data <- as.data.frame(matrix(c(0.55,0.81,0.91,0.6), ncol=4))
colnames(data) <- c("Dimension 1 (from ...)" , "Dimension 2 (from ...)" , "Dimension 3 (from ...)" , "Dimension (from ...)")
 
# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
data <- rbind(rep(1,4) , rep(0,4) , data)
data<-rbind(data,rbeta(4,1,2))


# Check your data, it has to look like this!
# head(data)

# The default radar chart 
par(mar=c(2, 0, 2, 0))
radarchart(data, vlabels=c("Dimension 1\n(from ...)", 
"Dimension 2\n(from ...)", 
"Dimension 3\n(from ...)", 
"Dimension 4\n(from ...)"),vlcex=0.65,axistype=4,axislabcol="black",calcex=0.7)
title("Radar chart of the performance \n(run level)",cex.main=0.85)


# The default radar chart 
par(mar=c(2, 0, 2, 0))
radarchart(data, vlabels=c("Dimension 1\n(from ...)", 
"Dimension 2\n(from ...)", 
"Dimension 3\n(from ...)", 
"Dimension 4\n(from -.)"),vlcex=0.65,axistype=4,axislabcol="black",calcex=0.7)
title("Radar chart of the performance \n(region level)",cex.main=0.85)


# The default radar chart 
par(mar=c(2, 0, 2, 0))
radarchart(data, vlabels=c("Dimension 1\n(from ...)", 
"Dimension 2\n(from ...)", 
"Dimension 3\n(from ...)", 
"Dimension 4\n(from -.)"),vlcex=0.65,axistype=4,axislabcol="black",calcex=0.7)
title("Radar chart of the performance \n(gene level)",cex.main=0.85)


# The default radar chart 
par(mar=c(2, 0, 2, 0))
radarchart(data, vlabels=c("Dimension 1\n(from ...)", 
"Dimension 2\n(from ...)", 
"Dimension 3\n(from ...)", 
"Dimension 4\n(from -.)"),vlcex=0.65,axistype=4,axislabcol="black",calcex=0.7)
title("Radar chart of the performance \n(variant level)",cex.main=0.85)



## The default radar chart 
#par(mar=c(2, 0, 2, 0))
#radarchart(data, vlabels=c("Dimension 1\n(from ...)", 
#"Dimension 2\n(from ...)", 
#"Dimension 3\n(from ...)", 
#"Dimension 4\n(from ...)"),vlcex=0.65,axistype=4,axislabcol="black",calcex=0.7)
#title("Radar chart of the performance \n(subject level)",cex.main=0.85)


summary(dati_paper_NGS$pct_chimeric_reads)
summary(dati_paper_NGS$pct_exon_100x)
summary(dati_paper_NGS$pct_read_enrichment)
summary(dati_paper_NGS$pct_aligned_reads)
summary(dati_paper_NGS$pct_target_100x)
summary(dati_paper_NGS$pct_target_250x)


length(na.omit(dati_paper_NGS$pct_chimeric_reads))
length(na.omit(dati_paper_NGS$pct_exon_100x))
length(na.omit(dati_paper_NGS$pct_read_enrichment))
length(na.omit(dati_paper_NGS$pct_aligned_reads))
length(na.omit(dati_paper_NGS$pct_target_100x))
length(na.omit(dati_paper_NGS$pct_target_250x))


inf_DNA<-rep(0,6); inf_DNA
sup_DNA<-rep(100,6); sup_DNA


med_DNA<-c(median(dati_paper_NGS$pct_chimeric_reads,na.rm=T),median(dati_paper_NGS$pct_exon_100x,na.rm=T),median(dati_paper_NGS$pct_read_enrichment,na.rm=T),
           median(dati_paper_NGS$pct_aligned_reads,na.rm=T),median(dati_paper_NGS$pct_target_100x,na.rm=T),median(dati_paper_NGS$pct_target_250x,na.rm=T)); med_DNA

min_DNA<-c(min(dati_paper_NGS$pct_chimeric_reads,na.rm=T),min(dati_paper_NGS$pct_exon_100x,na.rm=T),min(dati_paper_NGS$pct_read_enrichment,na.rm=T),
           min(dati_paper_NGS$pct_aligned_reads,na.rm=T),min(dati_paper_NGS$pct_target_100x,na.rm=T),min(dati_paper_NGS$pct_target_250x,na.rm=T)); min_DNA

max_DNA<-c(max(dati_paper_NGS$pct_chimeric_reads,na.rm=T),max(dati_paper_NGS$pct_exon_100x,na.rm=T),max(dati_paper_NGS$pct_read_enrichment,na.rm=T),
           max(dati_paper_NGS$pct_aligned_reads,na.rm=T),max(dati_paper_NGS$pct_target_100x,na.rm=T),max(dati_paper_NGS$pct_target_250x,na.rm=T)); max_DNA


data_DNA <- as.data.frame(matrix(sup_DNA, ncol=6))
colnames(data_DNA) <- c("Chimeric reads (%)" , "Exon 100x (%)" , "Read enrichment (%)" , 
"Aligned reads (%)", "Target 100x (%)", "Target 250x (%)")
 
data_DNA

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
data_DNA <- rbind(data_DNA,inf_DNA,min_DNA, med_DNA, max_DNA); data_DNA


# The default radar chart 
#par(mfrow=c(1,2), mar=c(2, 0, 2, 0))
radarchart(data_DNA, vlabels=c("Chimeric reads (%)", 
"Exon 100x (%)", 
"Read enrichment (%)", 
"Aligned reads (%)",
"Target 100x (%)",
"Target 250x (%)"),vlcex=1.3,axistype=4,axislabcol="black",calcex=0.8,caxislabels=seq(0,100,25),plwd=4)
title("Radar chart of the DNA performance \n(sample level, n=71)",cex.main=1)
legend("bottomleft",lty=c(1,2,3),col=c(1,2,3),lwd=c(2,2,2),c("min","median","max"),cex=1)




#############################################################################################################
#####Radar plot per le features in giallo di RNA (con anche la feature RNA median cv gene 500x perc)#########
#############################################################################################################

summary(dati_paper_NGS$rna_pct_chimeric_reads)
summary(dati_paper_NGS$rna_median_cv__gene_500x_perc)
summary(dati_paper_NGS$rna_pct_on_target_reads)

length(na.omit(dati_paper_NGS$rna_pct_chimeric_reads))
length(na.omit(dati_paper_NGS$rna_median_cv__gene_500x_perc))
length(na.omit(dati_paper_NGS$rna_pct_on_target_reads))


inf_RNA<-rep(0,3); inf_RNA
sup_RNA<-rep(100,3); sup_RNA

med_RNA<-c(median(dati_paper_NGS$rna_pct_chimeric_reads,na.rm=T),median(dati_paper_NGS$rna_median_cv__gene_500x_perc,na.rm=T),median(dati_paper_NGS$rna_pct_on_target_reads,na.rm=T))
med_RNA


min_RNA<-c(min(dati_paper_NGS$rna_pct_chimeric_reads,na.rm=T),min(dati_paper_NGS$rna_median_cv__gene_500x_perc,na.rm=T),min(dati_paper_NGS$rna_pct_on_target_reads,na.rm=T))
min_RNA


max_RNA<-c(max(dati_paper_NGS$rna_pct_chimeric_reads,na.rm=T),max(dati_paper_NGS$rna_median_cv__gene_500x_perc,na.rm=T),max(dati_paper_NGS$rna_pct_on_target_reads,na.rm=T))
max_RNA


data_RNA <- as.data.frame(matrix(sup_RNA, ncol=3))
colnames(data_RNA) <- c("Chimeric reads (%)", "Gene median coverage 500x (%)", "Target reads (%)")

data_RNA

# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each topic to show on the plot!
data_RNA <- rbind(data_RNA,inf_RNA,min_RNA, med_RNA, max_RNA); data_RNA



# The default radar chart 
#par(mar=c(2, 0, 2, 0))
radarchart(data_RNA, vlabels=c("Chimeric reads (%)", 
"Gene median coverage 500x (%)", "Target reads (%)"),vlcex=1.3,axistype=4,axislabcol="black",calcex=0.8,caxislabels=seq(0,100,25),plwd=4)
title("Radar chart of the RNA performance \n(sample level, n=59)",cex.main=1)
legend("bottomleft",lty=c(1,2,3),col=c(1,2,3),lwd=c(2,2,2),c("min","median","max"),cex=1)



#######################
###COVERAGE ANALYSIS###
#######################

#####
#SNV#
#####

dati_cov_SNV<-read.csv2(file.choose(),h=T,sep=";",dec=",",na.strings="") #importo TSO500_SNV_coverage_byPatient.csv
dim(dati_cov_SNV) # 71 28
head(dati_cov_SNV)
colnames(dati_cov_SNV)
summary(dati_cov_SNV)

par(mfrow=c(2,2))
#Q2x50
summary(dati_cov_SNV$Q2x50); table(dati_cov_SNV$Q2x50>=75,useNA="ifany") 
sd(dati_cov_SNV$Q2x50,na.rm=T)
tapply(dati_cov_SNV$Q2x50,dati_cov_SNV$Run,summary)
tapply(dati_cov_SNV$Q2x50,dati_cov_SNV$Run,function(x) sd(x,na.rm=T))
tapply(dati_cov_SNV$Q2x50,dati_cov_SNV$Run,length)
tapply(dati_cov_SNV$Q2x50,dati_cov_SNV$Run,function(x) length(x[is.na(x)==T]))
boxplot(dati_cov_SNV$Q2x50~dati_cov_SNV$Run, main="Median SNV coverage percent at 50X depth by run",
ylab="Median SNV coverage at 50X depth by run (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
abline(h=75,lwd=2,lty=2)
kruskal.test(dati_cov_SNV$Q2x50~dati_cov_SNV$Run)
pairwise.wilcox.test(dati_cov_SNV$Q2x50, dati_cov_SNV$Run, p.adj= "none")
pairwise.wilcox.test(dati_cov_SNV$Q2x50, dati_cov_SNV$Run, p.adj= "fdr")
pairwise.wilcox.test(dati_cov_SNV$Q2x50, dati_cov_SNV$Run, p.adj= "BH")
#tapply(dati_cov_SNV$Q2x50,dati_cov_SNA$Run,function(x) table(x>=75,useNA="ifany"))


#Q2x100
summary(dati_cov_SNV$Q2x100); table(dati_cov_SNV$Q2x100>=75,useNA="ifany") 
sd(dati_cov_SNV$Q2x100,na.rm=T)
tapply(dati_cov_SNV$Q2x100,dati_cov_SNV$Run,summary)
tapply(dati_cov_SNV$Q2x100,dati_cov_SNV$Run,function(x) sd(x,na.rm=T))
tapply(dati_cov_SNV$Q2x100,dati_cov_SNV$Run,length)
tapply(dati_cov_SNV$Q2x100,dati_cov_SNV$Run,function(x) length(x[is.na(x)==T]))
boxplot(dati_cov_SNV$Q2x100~dati_cov_SNV$Run, main="Median SNV coverage percent at 100X depth by run",
ylab="Median SNV coverage at 100X depth by run (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
abline(h=75,lwd=2,lty=2)
kruskal.test(dati_cov_SNV$Q2x100~dati_cov_SNV$Run)
library(coin)
kruskal_test(dati_cov_SNV$Q2x100~factor(dati_cov_SNV$Run))
pairwise.wilcox.test(dati_cov_SNV$Q2x100, dati_cov_SNV$Run, p.adj= "none")
pairwise.wilcox.test(dati_cov_SNV$Q2x100, dati_cov_SNV$Run, p.adj= "fdr")
pairwise.wilcox.test(dati_cov_SNV$Q2x100, dati_cov_SNV$Run, p.adj= "BH")
#tapply(dati_cov_SNV$Q2x100,dati_cov_SNV$Run,function(x) table(x>=75,useNA="ifany"))


#install.packages("agricolae")
library(agricolae)

Median.test(dati_cov_SNV$Q2x100,dati_cov_SNV$Run,alpha=0.05, group = FALSE)

library(coin)
median_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==2]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==2]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==3]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==3]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==4]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==4]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==5]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==5]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==3]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==3]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==4]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==4]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==5]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==5]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==4]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==4]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==5]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==5]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==4 | dati_cov_SNV$Run==5]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==4 | dati_cov_SNV$Run==5]),
data = dati_cov_SNV, distribution = "exact")

###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr


library(coin)
wilcox_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==2]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==2]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==3]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==4]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==5]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==3]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==4]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==5]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==4]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==5]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_SNV$Q2x100[dati_cov_SNV$Run==4 | dati_cov_SNV$Run==5]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==4 | dati_cov_SNV$Run==5]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL TEST DI MANN-WHITNEY
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



#Q2x250
summary(dati_cov_SNV$Q2x250); table(dati_cov_SNV$Q2x250>=75,useNA="ifany") 
sd(dati_cov_SNV$Q2x250,na.rm=T)
tapply(dati_cov_SNV$Q2x250,dati_cov_SNV$Run,summary)
tapply(dati_cov_SNV$Q2x250,dati_cov_SNV$Run,function(x) sd(x,na.rm=T))
tapply(dati_cov_SNV$Q2x250,dati_cov_SNV$Run,length)
tapply(dati_cov_SNV$Q2x250,dati_cov_SNV$Run,function(x) length(x[is.na(x)==T]))
boxplot(dati_cov_SNV$Q2x250~dati_cov_SNV$Run, main="Median SNV coverage percent at 250X depth by run",
ylab="Median SNV coverage at 250X depth by run (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
abline(h=75,lwd=2,lty=2)
kruskal.test(dati_cov_SNV$Q2x250~dati_cov_SNV$Run)
pairwise.wilcox.test(dati_cov_SNV$Q2x250, dati_cov_SNV$Run, p.adj= "none")
pairwise.wilcox.test(dati_cov_SNV$Q2x250, dati_cov_SNV$Run, p.adj= "fdr")
pairwise.wilcox.test(dati_cov_SNV$Q2x250, dati_cov_SNV$Run, p.adj= "BH")
#tapply(dati_cov_SNV$Q2x250,dati_cov_SNV$Run,function(x) table(x>=75,useNA="ifany"))

#install.packages("agricolae")
library(agricolae)

#Median.test(dati_cov_SNV$Q2x250,factor(dati_cov_SNV$Run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==2]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==2]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==3]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==3]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==4]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==4]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==5]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==5]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==3]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==3]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==4]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==4]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==5]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==5]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==4]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==4]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==5]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==5]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==4 | dati_cov_SNV$Run==5]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==4 | dati_cov_SNV$Run==5]),
data = dati_cov_SNV, distribution = "exact")

###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



library(coin)
wilcox_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==2]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==2]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==3]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==4]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==5]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==3]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==4]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==5]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==4]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==5]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_SNV$Q2x250[dati_cov_SNV$Run==4 | dati_cov_SNV$Run==5]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==4 | dati_cov_SNV$Run==5]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



#Q2x500
summary(dati_cov_SNV$Q2x500); table(dati_cov_SNV$Q2x500>=75,useNA="ifany") 
sd(dati_cov_SNV$Q2x500,na.rm=T)
tapply(dati_cov_SNV$Q2x500,dati_cov_SNV$Run,summary)
tapply(dati_cov_SNV$Q2x500,dati_cov_SNV$Run,function(x) sd(x,na.rm=T))
tapply(dati_cov_SNV$Q2x500,dati_cov_SNV$Run,length)
tapply(dati_cov_SNV$Q2x500,dati_cov_SNV$Run,function(x) length(x[is.na(x)==T]))
boxplot(dati_cov_SNV$Q2x500~dati_cov_SNV$Run, main="Median SNV coverage percent at 500X depth by run",
ylab="Median SNV coverage at 500X depth by run (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
abline(h=75,lwd=2,lty=2)
kruskal.test(dati_cov_SNV$Q2x500~dati_cov_SNV$Run)
pairwise.wilcox.test(dati_cov_SNV$Q2x500, dati_cov_SNV$Run, p.adj= "none")
pairwise.wilcox.test(dati_cov_SNV$Q2x500, dati_cov_SNV$Run, p.adj= "fdr")
pairwise.wilcox.test(dati_cov_SNV$Q2x500, dati_cov_SNV$Run, p.adj= "BH")
#tapply(dati_cov_SNV$Q2x500,dati_cov_SNV$Run,function(x) table(x>=75,useNA="ifany"))


#Median.test(dati_cov_SNV$Q2x500,factor(dati_cov_SNV$Run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==2]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==2]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==3]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==3]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==4]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==4]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==5]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==5]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==3]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==3]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==4]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==4]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==5]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==5]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==4]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==4]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==5]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==5]),
data = dati_cov_SNV, distribution = "exact")

median_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==4 | dati_cov_SNV$Run==5]~
factor(dati_cov_SNV$Run[dati_cov_SNV$Run==4 | dati_cov_SNV$Run==5]),
data = dati_cov_SNV, distribution = "exact")

###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



library(coin)
wilcox_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==2]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==2]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==3]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==4]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==5]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==1 | dati_cov_SNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==3]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==4]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==5]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==2 | dati_cov_SNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==4]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==5]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==3 | dati_cov_SNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_SNV$Q2x500[dati_cov_SNV$Run==4 | dati_cov_SNV$Run==5]~factor(dati_cov_SNV$Run[dati_cov_SNV$Run==4 | dati_cov_SNV$Run==5]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



###DEPTH###
depth_SNV<-c(rep(50,71),rep(100,71),rep(250,71),rep(500,71))
depth_level_SNV<-c(dati_cov_SNV$Q2x50,dati_cov_SNV$Q2x100,dati_cov_SNV$Q2x250,dati_cov_SNV$Q2x500)
tapply(depth_level_SNV,depth_SNV,summary)
tapply(depth_level_SNV,depth_SNV,function(x) sd(x,na.rm=T))

kruskal.test(depth_level_SNV~depth_SNV)

boxplot(depth_level_SNV~depth_SNV, main="Median SNV coverage percent by depth",
ylab="Median SNV coverage by depth (%)", xlab="Depth",axes=F)
axis(1, at=c(1:4), labels=c("50X (n=71)","100X (n=71)","250X (n=71)","500X (n=71)"))
axis(2)
abline(h=75,lwd=2,lty=2)
pairwise.wilcox.test(depth_level_SNV, depth_SNV, p.adj= "none")
pairwise.wilcox.test(depth_level_SNV, depth_SNV, p.adj= "fdr")
pairwise.wilcox.test(depth_level_SNV, depth_SNV, p.adj= "BH")



Median.test(depth_level_SNV,depth_SNV,alpha=0.05, group = FALSE)

library(coin)
median_test(depth_level_SNV[depth_SNV==50 | depth_SNV==100]~
factor(depth_SNV[depth_SNV==50 | depth_SNV==100]), distribution = "exact")

median_test(depth_level_SNV[depth_SNV==50 | depth_SNV==250]~
factor(depth_SNV[depth_SNV==50 | depth_SNV==250]), distribution = "exact")

median_test(depth_level_SNV[depth_SNV==50 | depth_SNV==500]~
factor(depth_SNV[depth_SNV==50 | depth_SNV==500]), distribution = "exact")

median_test(depth_level_SNV[depth_SNV==100 | depth_SNV==250]~
factor(depth_SNV[depth_SNV==100 | depth_SNV==250]), distribution = "exact")

median_test(depth_level_SNV[depth_SNV==100 | depth_SNV==500]~
factor(depth_SNV[depth_SNV==100 | depth_SNV==500]), distribution = "exact")

median_test(depth_level_SNV[depth_SNV==250 | depth_SNV==500]~
factor(depth_SNV[depth_SNV==250 | depth_SNV==500]), distribution = "exact")

###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 6)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr


library(coin)
wilcox_test(depth_level_SNV[depth_SNV==50 | depth_SNV==100]~factor(depth_SNV[depth_SNV==50 | depth_SNV==100]),distribution = "exact", conf.int = TRUE)
wilcox_test(depth_level_SNV[depth_SNV==50 | depth_SNV==250]~factor(depth_SNV[depth_SNV==50 | depth_SNV==250]),distribution = "exact", conf.int = TRUE)
wilcox_test(depth_level_SNV[depth_SNV==50 | depth_SNV==500]~factor(depth_SNV[depth_SNV==50 | depth_SNV==500]),distribution = "exact", conf.int = TRUE)

wilcox_test(depth_level_SNV[depth_SNV==100 | depth_SNV==250]~factor(depth_SNV[depth_SNV==100 | depth_SNV==250]),distribution = "exact", conf.int = TRUE)
wilcox_test(depth_level_SNV[depth_SNV==100 | depth_SNV==500]~factor(depth_SNV[depth_SNV==100 | depth_SNV==500]),distribution = "exact", conf.int = TRUE)

wilcox_test(depth_level_SNV[depth_SNV==250 | depth_SNV==500]~factor(depth_SNV[depth_SNV==250 | depth_SNV==500]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL TEST DI MANN-WHITNEY
p.adjust(p, method = "fdr", n = 6)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr




#####
#CNV#
#####

dati_cov_CNV<-read.csv2(file.choose(),h=T,sep=";",dec=",",na.strings="") #importo TSO500_CNV_coverage_byPatient.csv
dim(dati_cov_CNV) #71 28
head(dati_cov_CNV)
colnames(dati_cov_CNV)
summary(dati_cov_CNV)


par(mfrow=c(2,2))
#Q2x50
summary(dati_cov_CNV$Q2x50); table(dati_cov_CNV$Q2x50>=75,useNA="ifany") 
sd(dati_cov_CNV$Q2x50,na.rm=T)
tapply(dati_cov_CNV$Q2x50,dati_cov_CNV$Run,summary)
tapply(dati_cov_CNV$Q2x50,dati_cov_CNV$Run,function(x) sd(x,na.rm=T))
tapply(dati_cov_CNV$Q2x50,dati_cov_CNV$Run,length)
tapply(dati_cov_CNV$Q2x50,dati_cov_CNV$Run,function(x) length(x[is.na(x)==T]))
boxplot(dati_cov_CNV$Q2x50~dati_cov_CNV$Run, main="Median CNV coverage percent at 50X depth by run",
ylab="Median CNV coverage at 50X depth by run (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
abline(h=75,lwd=2,lty=2)
kruskal.test(dati_cov_CNV$Q2x50~dati_cov_CNV$Run)
pairwise.wilcox.test(dati_cov_CNV$Q2x50, dati_cov_CNV$Run, p.adj= "none")
pairwise.wilcox.test(dati_cov_CNV$Q2x50, dati_cov_CNV$Run, p.adj= "fdr")
pairwise.wilcox.test(dati_cov_CNV$Q2x50, dati_cov_CNV$Run, p.adj= "BH")
#tapply(dati_cov_CNV$Q2x50,dati_cov_RNA$Run,function(x) table(x>=75,useNA="ifany"))


#Median.test(dati_cov_CNV$Q2x50,factor(dati_cov_CNV$Run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")


###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



library(coin)
wilcox_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_SNV$Run==1 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_CNV$Q2x50[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



#Q2x100
summary(dati_cov_CNV$Q2x100); table(dati_cov_CNV$Q2x100>=75,useNA="ifany") 
sd(dati_cov_CNV$Q2x100,na.rm=T)
tapply(dati_cov_CNV$Q2x100,dati_cov_CNV$Run,summary)
tapply(dati_cov_CNV$Q2x100,dati_cov_CNV$Run,function(x) sd(x,na.rm=T))
tapply(dati_cov_CNV$Q2x100,dati_cov_CNV$Run,length)
tapply(dati_cov_CNV$Q2x100,dati_cov_CNV$Run,function(x) length(x[is.na(x)==T]))
boxplot(dati_cov_CNV$Q2x100~dati_cov_CNV$Run, main="Median CNV coverage percent at 100X depth by run",
ylab="Median CNV coverage at 100X depth by run (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
abline(h=75,lwd=2,lty=2)
kruskal.test(dati_cov_CNV$Q2x100~dati_cov_CNV$Run)
pairwise.wilcox.test(dati_cov_CNV$Q2x100, dati_cov_CNV$Run, p.adj= "none")
pairwise.wilcox.test(dati_cov_CNV$Q2x100, dati_cov_CNV$Run, p.adj= "fdr")
pairwise.wilcox.test(dati_cov_CNV$Q2x100, dati_cov_CNV$Run, p.adj= "BH")
#tapply(dati_cov_CNV$Q2x100,dati_cov_CNV$Run,function(x) table(x>=75,useNA="ifany"))


#Median.test(dati_cov_CNV$Q2x100,factor(dati_cov_CNV$Run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")


###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr


library(coin)
wilcox_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_SNV$Run==1 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_CNV$Q2x100[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)


###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr




#Q2x250
summary(dati_cov_CNV$Q2x250); table(dati_cov_CNV$Q2x250>=75,useNA="ifany") 
sd(dati_cov_CNV$Q2x250,na.rm=T)
tapply(dati_cov_CNV$Q2x250,dati_cov_CNV$Run,summary)
tapply(dati_cov_CNV$Q2x250,dati_cov_CNV$Run,function(x) sd(x,na.rm=T))
tapply(dati_cov_CNV$Q2x250,dati_cov_CNV$Run,length)
tapply(dati_cov_CNV$Q2x250,dati_cov_CNV$Run,function(x) length(x[is.na(x)==T]))
boxplot(dati_cov_CNV$Q2x250~dati_cov_CNV$Run, main="Median CNV coverage percent at 250X depth by run",
ylab="Median CNV coverage at 250X depth by run (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
abline(h=75,lwd=2,lty=2)
kruskal.test(dati_cov_CNV$Q2x250~dati_cov_CNV$Run)
pairwise.wilcox.test(dati_cov_CNV$Q2x250, dati_cov_CNV$Run, p.adj= "none")
pairwise.wilcox.test(dati_cov_CNV$Q2x250, dati_cov_CNV$Run, p.adj= "fdr")
pairwise.wilcox.test(dati_cov_CNV$Q2x250, dati_cov_CNV$Run, p.adj= "BH")
#tapply(dati_cov_CNV$Q2x250,dati_cov_CNV$Run,function(x) table(x>=75,useNA="ifany"))

#Median.test(dati_cov_CNV$Q2x250,factor(dati_cov_CNV$Run),alpha=0.05, group = FALSE)


library(coin)
median_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")


###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr


library(coin)
wilcox_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_SNV$Run==1 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_CNV$Q2x250[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)



###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr




#Q2x500
summary(dati_cov_CNV$Q2x500); table(dati_cov_CNV$Q2x500>=75,useNA="ifany") 
sd(dati_cov_CNV$Q2x500,na.rm=T)
tapply(dati_cov_CNV$Q2x500,dati_cov_CNV$Run,summary)
tapply(dati_cov_CNV$Q2x500,dati_cov_CNV$Run,function(x) sd(x,na.rm=T))
tapply(dati_cov_CNV$Q2x500,dati_cov_CNV$Run,length)
tapply(dati_cov_CNV$Q2x500,dati_cov_CNV$Run,function(x) length(x[is.na(x)==T]))
boxplot(dati_cov_CNV$Q2x500~dati_cov_CNV$Run, main="Median CNV coverage percent at 500X depth by run",
ylab="Median CNV coverage at 500X depth by run (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=8)","2 (n=8)","3 (n=8)","4 (n=16)","5 (n=31)"))
axis(2)
abline(h=75,lwd=2,lty=2)
kruskal.test(dati_cov_CNV$Q2x500~dati_cov_CNV$Run)
pairwise.wilcox.test(dati_cov_CNV$Q2x500, dati_cov_CNV$Run, p.adj= "none")
pairwise.wilcox.test(dati_cov_CNV$Q2x500, dati_cov_CNV$Run, p.adj= "fdr")
pairwise.wilcox.test(dati_cov_CNV$Q2x500, dati_cov_CNV$Run, p.adj= "BH")
#tapply(dati_cov_CNV$Q2x500,dati_cov_CNV$Run,function(x) table(x>=75,useNA="ifany"))


#Median.test(dati_cov_CNV$Q2x500,factor(dati_cov_CNV$Run),alpha=0.05, group = FALSE)


library(coin)
median_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")

median_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]~
factor(dati_cov_CNV$Run[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]),
data = dati_cov_CNV, distribution = "exact")


###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



library(coin)
wilcox_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==2]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==1 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_SNV$Run==1 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==2 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==3 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_CNV$Q2x500[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]~factor(dati_cov_CNV$Run[dati_cov_CNV$Run==4 | dati_cov_CNV$Run==5]),distribution = "exact", conf.int = TRUE)


###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr





###DEPTH###
depth_CNV<-c(rep(50,71),rep(100,71),rep(250,71),rep(500,71))
depth_level_CNV<-c(dati_cov_CNV$Q2x50,dati_cov_CNV$Q2x100,dati_cov_CNV$Q2x250,dati_cov_CNV$Q2x500)
tapply(depth_level_CNV,depth_CNV,summary)
tapply(depth_level_CNV,depth_CNV,function(x) sd(x,na.rm=T))

kruskal.test(depth_level_CNV~depth_CNV)

boxplot(depth_level_CNV~depth_CNV, main="Median CNV coverage percent by depth",
ylab="Median CNV coverage by depth (%)", xlab="Depth",axes=F)
axis(1, at=c(1:4), labels=c("50X (n=71)","100X (n=71)","250X (n=71)","500X (n=71)"))
axis(2)
abline(h=75,lwd=2,lty=2)
pairwise.wilcox.test(depth_level_CNV, depth_CNV, p.adj= "none")
pairwise.wilcox.test(depth_level_CNV, depth_CNV, p.adj= "fdr")
pairwise.wilcox.test(depth_level_CNV, depth_CNV, p.adj= "BH")

library(agricolae)
Median.test(depth_level_CNV,depth_CNV,alpha=0.05, group = FALSE)

library(coin)
median_test(depth_level_CNV[depth_CNV==50 | depth_CNV==100]~
factor(depth_CNV[depth_CNV==50 | depth_CNV==100]), distribution = "exact")

median_test(depth_level_CNV[depth_CNV==50 | depth_CNV==250]~
factor(depth_CNV[depth_CNV==50 | depth_CNV==250]), distribution = "exact")

median_test(depth_level_CNV[depth_CNV==50 | depth_CNV==500]~
factor(depth_CNV[depth_CNV==50 | depth_CNV==500]), distribution = "exact")

median_test(depth_level_CNV[depth_CNV==100 | depth_CNV==250]~
factor(depth_CNV[depth_CNV==100 | depth_CNV==250]), distribution = "exact")

median_test(depth_level_CNV[depth_CNV==100 | depth_CNV==500]~
factor(depth_CNV[depth_CNV==100 | depth_CNV==500]), distribution = "exact")

median_test(depth_level_CNV[depth_CNV==250 | depth_CNV==500]~
factor(depth_CNV[depth_CNV==250 | depth_CNV==500]), distribution = "exact")


###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 6)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



library(coin)
wilcox_test(depth_level_CNV[depth_CNV==50 | depth_CNV==100]~factor(depth_CNV[depth_CNV==50 | depth_CNV==100]),distribution = "exact", conf.int = TRUE)
wilcox_test(depth_level_CNV[depth_CNV==50 | depth_CNV==250]~factor(depth_CNV[depth_CNV==50 | depth_CNV==250]),distribution = "exact", conf.int = TRUE)
wilcox_test(depth_level_CNV[depth_CNV==50 | depth_CNV==500]~factor(depth_CNV[depth_CNV==50 | depth_CNV==500]),distribution = "exact", conf.int = TRUE)

wilcox_test(depth_level_CNV[depth_CNV==100 | depth_CNV==250]~factor(depth_CNV[depth_CNV==100 | depth_CNV==250]),distribution = "exact", conf.int = TRUE)
wilcox_test(depth_level_CNV[depth_CNV==100 | depth_CNV==500]~factor(depth_CNV[depth_CNV==100 | depth_CNV==500]),distribution = "exact", conf.int = TRUE)

wilcox_test(depth_level_CNV[depth_CNV==250 | depth_CNV==500]~factor(depth_CNV[depth_CNV==250 | depth_CNV==500]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 6)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr


#####
#RNA#
#####

dati_cov_RNA<-read.csv2(file.choose(),h=T,sep=";",dec=",",na.strings="") #importo TSO500_RNA_sampleCoverage_Validation_V1.csv
dim(dati_cov_RNA) #59 40
head(dati_cov_RNA)
colnames(dati_cov_RNA)
summary(dati_cov_RNA)


par(mfrow=c(1,3))
#Q2x5
summary(dati_cov_RNA$Q2x5); table(dati_cov_RNA$Q2x5>=75,useNA="ifany") 
sd(dati_cov_RNA$Q2x5,na.rm=T)
tapply(dati_cov_RNA$Q2x5,dati_cov_RNA$Run,summary)
tapply(dati_cov_RNA$Q2x5,dati_cov_RNA$Run,function(x) sd(x,na.rm=T))
tapply(dati_cov_RNA$Q2x5,dati_cov_RNA$Run,length)
tapply(dati_cov_RNA$Q2x5,dati_cov_RNA$Run,function(x) length(x[is.na(x)==T]))
boxplot(dati_cov_RNA$Q2x5~dati_cov_RNA$Run, main="Median RNA coverage percent at 5X depth by run",
ylab="Median RNA coverage at 5X depth by run (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
abline(h=75,lwd=2,lty=2)
kruskal.test(dati_cov_RNA$Q2x5~dati_cov_RNA$Run)
pairwise.wilcox.test(dati_cov_RNA$Q2x5, dati_cov_RNA$Run, p.adj= "none")
pairwise.wilcox.test(dati_cov_RNA$Q2x5, dati_cov_RNA$Run, p.adj= "fdr")
pairwise.wilcox.test(dati_cov_RNA$Q2x5, dati_cov_RNA$Run, p.adj= "BH")
#tapply(dati_cov_RNA$Q2x5,dati_cov_RNA$Run,function(x) table(x>=75,useNA="ifany"))


#Median.test(dati_cov_RNA$Q2x5,factor(dati_cov_RNA$Run),alpha=0.05, group = FALSE)


library(coin)
median_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==2]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==2]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==3]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==3]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==4]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==4]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==5]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==5]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==3]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==3]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==4]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==4]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==5]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==5]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==4]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==4]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==5]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==5]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==4 | dati_cov_RNA$Run==5]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==4 | dati_cov_RNA$Run==5]),
data = dati_cov_RNA, distribution = "exact")


###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr




library(coin)
wilcox_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==2]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==2]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==3]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==4]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==5]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==3]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==4]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==5]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==4]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==5]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_RNA$Q2x5[dati_cov_RNA$Run==4 | dati_cov_RNA$Run==5]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==4 | dati_cov_RNA$Run==5]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr





#Q2x10
summary(dati_cov_RNA$Q2x10); table(dati_cov_RNA$Q2x10>=75,useNA="ifany") 
sd(dati_cov_RNA$Q2x10,na.rm=T)
tapply(dati_cov_RNA$Q2x10,dati_cov_RNA$Run,summary)
tapply(dati_cov_RNA$Q2x10,dati_cov_RNA$Run,function(x) sd(x,na.rm=T))
tapply(dati_cov_RNA$Q2x10,dati_cov_RNA$Run,length)
tapply(dati_cov_RNA$Q2x10,dati_cov_RNA$Run,function(x) length(x[is.na(x)==T]))
boxplot(dati_cov_RNA$Q2x10~dati_cov_RNA$Run, main="Median RNA coverage percent at 10X depth by run",
ylab="Median RNA coverage at 10X depth by run (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
abline(h=75,lwd=2,lty=2)
kruskal.test(dati_cov_RNA$Q2x10~dati_cov_RNA$Run)
pairwise.wilcox.test(dati_cov_RNA$Q2x10, dati_cov_RNA$Run, p.adj= "none")
pairwise.wilcox.test(dati_cov_RNA$Q2x10, dati_cov_RNA$Run, p.adj= "fdr")
pairwise.wilcox.test(dati_cov_RNA$Q2x10, dati_cov_RNA$Run, p.adj= "BH")
#tapply(dati_cov_RNA$Q2x10,dati_cov_RNA$Run,function(x) table(x>=75,useNA="ifany"))


#Median.test(dati_cov_RNA$Q2x10,factor(dati_cov_RNA$Run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==2]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==2]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==3]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==3]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==4]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==4]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==5]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==5]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==3]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==3]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==4]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==4]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==5]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==5]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==4]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==4]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==5]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==5]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==4 | dati_cov_RNA$Run==5]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==4 | dati_cov_RNA$Run==5]),
data = dati_cov_RNA, distribution = "exact")


###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr




library(coin)
wilcox_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==2]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==2]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==3]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==4]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==5]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==3]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==4]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==5]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==4]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==5]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_RNA$Q2x10[dati_cov_RNA$Run==4 | dati_cov_RNA$Run==5]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==4 | dati_cov_RNA$Run==5]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr




#Q2x50
summary(dati_cov_RNA$Q2x50); table(dati_cov_RNA$Q2x50>=75,useNA="ifany") 
sd(dati_cov_RNA$Q2x50,na.rm=T)
tapply(dati_cov_RNA$Q2x50,dati_cov_RNA$Run,summary)
tapply(dati_cov_RNA$Q2x50,dati_cov_RNA$Run,function(x) sd(x,na.rm=T))
tapply(dati_cov_RNA$Q2x50,dati_cov_RNA$Run,length)
tapply(dati_cov_RNA$Q2x50,dati_cov_RNA$Run,function(x) length(x[is.na(x)==T]))
boxplot(dati_cov_RNA$Q2x50~dati_cov_RNA$Run, main="Median RNA coverage percent at 50X depth by run",
ylab="Median RNA coverage at 50X depth by run (%)", xlab="RUN",axes=F)
axis(1, at=c(1:5), labels=c("1 (n=6)","2 (n=6)","3 (n=8)","4 (n=15)","5 (n=24)"))
axis(2)
abline(h=75,lwd=2,lty=2)
kruskal.test(dati_cov_RNA$Q2x50~dati_cov_RNA$Run)
pairwise.wilcox.test(dati_cov_RNA$Q2x50, dati_cov_RNA$Run, p.adj= "none")
pairwise.wilcox.test(dati_cov_RNA$Q2x50, dati_cov_RNA$Run, p.adj= "fdr")
pairwise.wilcox.test(dati_cov_RNA$Q2x50, dati_cov_RNA$Run, p.adj= "BH")
#tapply(dati_cov_RNA$Q2x50,dati_cov_RNA$Run,function(x) table(x>=75,useNA="ifany"))


#Median.test(dati_cov_RNA$Q2x50,factor(dati_cov_RNA$Run),alpha=0.05, group = FALSE)

library(coin)
median_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==2]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==2]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==3]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==3]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==4]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==4]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==5]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==5]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==3]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==3]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==4]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==4]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==5]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==5]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==4]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==4]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==5]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==5]),
data = dati_cov_RNA, distribution = "exact")

median_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==4 | dati_cov_RNA$Run==5]~
factor(dati_cov_RNA$Run[dati_cov_RNA$Run==4 | dati_cov_RNA$Run==5]),
data = dati_cov_RNA, distribution = "exact")


###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr




library(coin)
wilcox_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==2]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==2]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==3]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==4]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==5]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==1 | dati_cov_RNA$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==3]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==3]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==4]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==5]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==2 | dati_cov_RNA$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==4]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==4]),distribution = "exact", conf.int = TRUE)
wilcox_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==5]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==3 | dati_cov_RNA$Run==5]),distribution = "exact", conf.int = TRUE)

wilcox_test(dati_cov_RNA$Q2x50[dati_cov_RNA$Run==4 | dati_cov_RNA$Run==5]~factor(dati_cov_RNA$Run[dati_cov_RNA$Run==4 | dati_cov_RNA$Run==5]),distribution = "exact", conf.int = TRUE)

###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 10)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr



###DEPTH###
depth_RNA<-c(rep(5,59),rep(10,59),rep(50,59))
depth_level_RNA<-c(dati_cov_RNA$Q2x5,dati_cov_RNA$Q2x10,dati_cov_RNA$Q2x50)
tapply(depth_level_RNA,depth_RNA,summary)
tapply(depth_level_RNA,depth_RNA,function(x) sd(x,na.rm=T))

kruskal.test(depth_level_RNA~depth_RNA)

boxplot(depth_level_RNA~depth_RNA, main="Median RNA coverage percent by depth",
ylab="Median RNA coverage by depth (%)", xlab="Depth",axes=F)
axis(1, at=c(1:3), labels=c("5X (n=59)","10X (n=59)","50X (n=59)"))
axis(2)
abline(h=75,lwd=2,lty=2)
pairwise.wilcox.test(depth_level_RNA, depth_RNA, p.adj= "none")
pairwise.wilcox.test(depth_level_RNA, depth_RNA, p.adj= "fdr")
pairwise.wilcox.test(depth_level_RNA, depth_RNA, p.adj= "BH")


Median.test(depth_level_RNA,depth_RNA,alpha=0.05, group = FALSE)

library(coin)
median_test(depth_level_RNA[depth_RNA==5 | depth_RNA==10]~
factor(depth_RNA[depth_RNA==5 | depth_RNA==10]), distribution = "exact")

median_test(depth_level_RNA[depth_RNA==5 | depth_RNA==50]~
factor(depth_RNA[depth_RNA==5 | depth_RNA==50]), distribution = "exact")

median_test(depth_level_RNA[depth_RNA==10 | depth_RNA==50]~
factor(depth_RNA[depth_RNA==10 | depth_RNA==50]), distribution = "exact")

###AGGIUSTAMENTO FDR DEL MEDIAN TEST
p.adjust(p, method = "fdr", n = 6)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr

library(coin)
wilcox_test(depth_level_RNA[depth_RNA==5 | depth_RNA==10]~factor(depth_RNA[depth_RNA==5 | depth_RNA==10]),distribution = "exact", conf.int = TRUE)
wilcox_test(depth_level_RNA[depth_RNA==5 | depth_RNA==50]~factor(depth_RNA[depth_RNA==5 | depth_RNA==50]),distribution = "exact", conf.int = TRUE)
wilcox_test(depth_level_RNA[depth_RNA==10 | depth_RNA==50]~factor(depth_RNA[depth_RNA==10 | depth_RNA==50]),distribution = "exact", conf.int = TRUE)


###AGGIUSTAMENTO FDR DEL MANN-WHITNEY TEST
p.adjust(p, method = "fdr", n = 6)
library(fdrtool)
fdr = fdrtool(pval$raw.p.val, statistic="pvalue")
fdr$qval # estimated Fdr values
fdr$lfdr # estimated local fdr




###ANALISI DI CORRELAZIONE###

dati_corr_wet_cov<-read.csv2(file.choose(),h=T,sep=";",dec=",",na.strings="") #importo WET & COVERAGE_250322.csv
dim(dati_corr_wet_cov) #130  56
head(dati_corr_wet_cov)
colnames(dati_corr_wet_cov)
summary(dati_corr_wet_cov)


library(Hmisc)
library(RVAideMemoire)
library(spearmanCI)
library(DescTools)

#rcorr(as.matrix(dati_corr_wet_cov[dati_corr_wet_cov$Source=="DNA",c(6:10,28:31,40:43)]),type="kendall")
#rcorr(as.matrix(dati[,10:18]),type="pearson")

#########
###DNA###
#########
library(DescTools)


#########
###SNV###
#########

#############
###OVERALL###
#############

###100X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)


###250X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)


###500X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)


############
###1a RUN###
############

###100X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)


###250X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)


###500X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
conf.level=0.95)


cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
conf.level=0.95)


cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
conf.level=0.95)


############
###2a RUN###
############

###100X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)


###250X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)


###500X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
conf.level=0.95)



############
###3a RUN###
############

###100X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)


###250X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)


###500X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
conf.level=0.95)


cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)



############
###4a RUN###
############

###100X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)


###250X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)


###500X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)



############
###5a RUN###
############

###100X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)


###250X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)


###500X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)


###Scatter plot###

colnames(dati_corr_wet_cov)


SNV<-data.frame(
run=factor(dati_corr_wet_cov$RUN[dati_corr_wet_cov$Source=="DNA"]),
dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"])


#pairs(SNV[,c(1,2,9)], col = SNV$run, oma=c(2,2,3.5,15),pch=19, main="Scatter plots", line.main=1.5)
#par(xpd = TRUE)
#legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)))


###SNV 50X###

par(mfrow=c(2,2))
plot(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Quantification (ng/ul)",ylab="SNV Coverage 50x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Quantification vs SNV Coverage 50x")

plot(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Protein Nanodrop (u)",ylab="SNV Coverage 50x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Protein Nanodrop vs SNV Coverage 50x")

plot(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Carbohydrate Nanodrop (u)",ylab="SNV Coverage 50x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Carbohydrate Nanodrop vs SNV Coverage 50x")

plot(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Qualification by RT-PCR (u)",ylab="SNV Coverage 50x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Qualification by RT-PCR vs SNV Coverage 50x")


par(mfrow=c(2,2))
plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Fragmentation by Sonication (bp)",ylab="SNV Coverage 50x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Fragmentation by Sonication vs SNV Coverage 50x")

plot(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA pre-hybridization (ng/ul)",ylab="SNV Coverage 50x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA pre-hybridization vs SNV Coverage 50x")

plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA post-hybridization (ng/ul)",ylab="SNV Coverage 50x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA post-hybridization vs SNV Coverage 50x")


###SNV 100X###

par(mfrow=c(2,2))
plot(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Quantification (ng/ul)",ylab="SNV Coverage 100x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Quantification vs SNV Coverage 100x")

plot(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Protein Nanodrop (u)",ylab="SNV Coverage 100x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Protein Nanodrop vs SNV Coverage 100x")

plot(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Carbohydrate Nanodrop (u)",ylab="SNV Coverage 100x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Carbohydrate Nanodrop vs SNV Coverage 100x")


plot(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Qualification by RT-PCR (u)",ylab="SNV Coverage 100x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Qualification by RT-PCR vs SNV Coverage 100x")


par(mfrow=c(1,3))
plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Fragmentation by Sonication (bp)",ylab="SNV Coverage 100x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Fragmentation by Sonication vs SNV Coverage 100x")

plot(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA pre-hybridization (ng/ul)",ylab="SNV Coverage 100x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA pre-hybridization vs SNV Coverage 100x")

plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA post-hybridization (ng/ul)",ylab="SNV Coverage 100x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA post-hybridization vs SNV Coverage 100x")




###SNV 250X###

par(mfrow=c(2,2))
plot(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Quantification (ng/ul)",ylab="SNV Coverage 250x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Quantification vs SNV Coverage 250x")

plot(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Protein Nanodrop (u)",ylab="SNV Coverage 250x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Protein Nanodrop vs SNV Coverage 250x")

plot(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Carbohydrate Nanodrop (u)",ylab="SNV Coverage 250x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Carbohydrate Nanodrop vs SNV Coverage 250x")


plot(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Qualification by RT-PCR (u)",ylab="SNV Coverage 250x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Qualification by RT-PCR vs SNV Coverage 250x")

par(mfrow=c(1,3))
plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Fragmentation by Sonication (bp)",ylab="SNV Coverage 250x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Fragmentation by Sonication vs SNV Coverage 250x")

plot(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA pre-hybridization (ng/ul)",ylab="SNV Coverage 250x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA pre-hybridization vs SNV Coverage 250x")

plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA post-hybridization (ng/ul)",ylab="SNV Coverage 250x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA post-hybridization vs SNV Coverage 250x")



###SNV 500X###

par(mfrow=c(2,2))
plot(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Quantification (ng/ul)",ylab="SNV Coverage 500x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Quantification vs SNV Coverage 500x")

plot(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Protein Nanodrop (u)",ylab="SNV Coverage 500x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Protein Nanodrop vs SNV Coverage 500x")

plot(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Carbohydrate Nanodrop (u)",ylab="SNV Coverage 500x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Carbohydrate Nanodrop vs SNV Coverage 500x")


plot(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Qualification by RT-PCR (u)",ylab="SNV Coverage 500x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Qualification by RT-PCR vs SNV Coverage 500x")

par(mfrow=c(1,3))
plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA Fragmentation by Sonication (bp)",ylab="SNV Coverage 500x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Fragmentation by Sonication vs SNV Coverage 500x")

plot(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA pre-hybridization (ng/ul)",ylab="SNV Coverage 500x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA pre-hybridization vs SNV Coverage 500x")

plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(SNV$run)[SNV$run],xlab="DNA post-hybridization (ng/ul)",ylab="SNV Coverage 500x (%)")
legend("bottomright", fill = unique(SNV$run), legend = c(levels(SNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA post-hybridization vs SNV Coverage 500x")





#########
###CNV###
#########

#############
###OVERALL###
#############

###100X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)


###250X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)


###500X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"], conf.level=0.95)


############
###1a RUN###
############

###100X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)


###250X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)


###500X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)



############
###2a RUN###
############

###100X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)


###250X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)


###500X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==2],
conf.level=0.95)



############
###3a RUN###
############

###100X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)


###250X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)


###500X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
conf.level=0.95)


cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)



############
###4a RUN###
############

###100X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)


###250X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)


###500X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)



############
###5a RUN###
############

###100X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)


###250X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)


###500X###

cor.test(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)



##################
###Scatter plot###
##################

CNV<-data.frame(
run=factor(dati_corr_wet_cov$RUN[dati_corr_wet_cov$Source=="DNA"]),
dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"])



par(mfrow=c(2,2))
plot(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Quantification (ng/ul)",ylab="CNV Coverage 50x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Quantification vs CNV Coverage 50x")

plot(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Protein Nanodrop (u)",ylab="CNV Coverage 50x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Protein Nanodrop vs CNV Coverage 50x")

plot(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Carbohydrate Nanodrop (u)",ylab="CNV Coverage 50x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Carbohydrate Nanodrop vs CNV Coverage 50x")

plot(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Qualification by RT-PCR (u)",ylab="CNV Coverage 50x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Qualification by RT-PCR vs CNV Coverage 50x")


par(mfrow=c(2,2))
plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Fragmentation by Sonication (bp)",ylab="CNV Coverage 50x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Fragmentation by Sonication vs CNV Coverage 50x")

plot(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA pre-hybridization (ng/ul)",ylab="CNV Coverage 50x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA pre-hybridization vs CNV Coverage 50x")

plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA post-hybridization (ng/ul)",ylab="CNV Coverage 50x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA post-hybridization vs CNV Coverage 50x")


###CNV 100X###

par(mfrow=c(2,2))
plot(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Quantification (ng/ul)",ylab="CNV Coverage 100x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Quantification vs CNV Coverage 100x")

plot(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Protein Nanodrop (u)",ylab="CNV Coverage 100x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Protein Nanodrop vs CNV Coverage 100x")

plot(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Carbohydrate Nanodrop (u)",ylab="CNV Coverage 100x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Carbohydrate Nanodrop vs CNV Coverage 100x")


plot(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Qualification by RT-PCR (u)",ylab="CNV Coverage 100x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Qualification by RT-PCR vs CNV Coverage 100x")


par(mfrow=c(1,3))
plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Fragmentation by Sonication (bp)",ylab="CNV Coverage 100x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Fragmentation by Sonication vs CNV Coverage 100x")

plot(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA pre-hybridization (ng/ul)",ylab="CNV Coverage 100x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA pre-hybridization vs CNV Coverage 100x")

plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA post-hybridization (ng/ul)",ylab="CNV Coverage 100x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA post-hybridization vs CNV Coverage 100x")




###SNV 250X###

par(mfrow=c(2,2))
plot(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Quantification (ng/ul)",ylab="CNV Coverage 250x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Quantification vs CNV Coverage 250x")

plot(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Protein Nanodrop (u)",ylab="CNV Coverage 250x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Protein Nanodrop vs CNV Coverage 250x")

plot(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Carbohydrate Nanodrop (u)",ylab="CNV Coverage 250x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Carbohydrate Nanodrop vs CNV Coverage 250x")


plot(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Qualification by RT-PCR (u)",ylab="CNV Coverage 250x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Qualification by RT-PCR vs CNV Coverage 250x")

par(mfrow=c(1,3))
plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Fragmentation by Sonication (bp)",ylab="CNV Coverage 250x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Fragmentation by Sonication vs CNV Coverage 250x")

plot(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA pre-hybridization (ng/ul)",ylab="CNV Coverage 250x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA pre-hybridization vs CNV Coverage 250x")

plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA post-hybridization (ng/ul)",ylab="CNV Coverage 250x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA post-hybridization vs CNV Coverage 250x")



###SNV 500X###

par(mfrow=c(2,2))
plot(dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Quantification (ng/ul)",ylab="CNV Coverage 500x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Quantification vs CNV Coverage 500x")

plot(dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Protein Nanodrop (u)",ylab="CNV Coverage 500x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Protein Nanodrop vs CNV Coverage 500x")

plot(dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Carbohydrate Nanodrop (u)",ylab="CNV Coverage 500x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Carbohydrate Nanodrop vs CNV Coverage 500x")


plot(dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Qualification by RT-PCR (u)",ylab="CNV Coverage 500x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Qualification by RT-PCR vs CNV Coverage 500x")

par(mfrow=c(1,3))
plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA Fragmentation by Sonication (bp)",ylab="CNV Coverage 500x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA Fragmentation by Sonication vs CNV Coverage 500x")

plot(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA pre-hybridization (ng/ul)",ylab="CNV Coverage 500x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA pre-hybridization vs CNV Coverage 500x")

plot(dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"], 
dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
pch = 19, col = levels(CNV$run)[CNV$run],xlab="DNA post-hybridization (ng/ul)",ylab="CNV Coverage 500x (%)")
legend("bottomright", fill = unique(CNV$run), legend = c(levels(CNV$run)),title="Run",cex=0.8)
title("Scatter plot DNA post-hybridization vs CNV Coverage 5




#########
###RNA###
#########

#############
###OVERALL###
#############

colnames(dati_corr_wet_cov)

###5X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)




###10X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"],
conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)


###50X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA"],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"], conf.level=0.95)



############
###1a RUN###
############


###5X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"  & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"  & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
conf.level=0.95)


cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)




###10X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)


###50X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==1],
conf.level=0.95)




############
###2a RUN###
############

###5X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"  & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"  & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)




###10X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
conf.level=0.95)


cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)


###50X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==2], conf.level=0.95)



############
###3a RUN###
############

###5X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"  & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"  & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"  & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
conf.level=0.95)


cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(
dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
conf.level=0.95)




###10X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)


###50X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==3], conf.level=0.95)



############
###4a RUN###
############

###5X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"  & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"  & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
conf.level=0.95)



###10X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
conf.level=0.95)


cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)


###50X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)
#source("Kendall_TauB_T_95%CI.R")
KendallTauB_T(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==4], conf.level=0.95)



############
###5a RUN###
############

###5X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"  & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"  & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)




###10X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)


###50X###

cor.test(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)

cor.test(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], method = "kendall",exact = FALSE)
KendallTauB(dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5],
dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA" & dati_corr_wet_cov$RUN==5], conf.level=0.95)




######################
###Correlation plot###
######################

colnames(dati_corr_wet_cov)


#######
##DNA##
#######

DNA_corr<-data.frame(DNA=dati_corr_wet_cov$DNA..ng.ul....3.5..12uL.[dati_corr_wet_cov$Source=="DNA"],
A260280=dati_corr_wet_cov$A.260.280[dati_corr_wet_cov$Source=="DNA"],
A260230=dati_corr_wet_cov$A.260.230[dati_corr_wet_cov$Source=="DNA"],
DeltaCt=dati_corr_wet_cov$Delta.Ct...5.[dati_corr_wet_cov$Source=="DNA"],
QC_post_frag=dati_corr_wet_cov$QC.post.Fragmentation.150.300bp[dati_corr_wet_cov$Source=="DNA"],
Pre_Hyb=dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="DNA"],
Post_Hyb=dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="DNA"],
#SNVQ2x50=dati_corr_wet_cov$SNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
SNVQ2x100=dati_corr_wet_cov$SNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
SNVQ2x250=dati_corr_wet_cov$SNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
SNVQ2x500=dati_corr_wet_cov$SNVQ2x500[dati_corr_wet_cov$Source=="DNA"],
#CNVQ2x50=dati_corr_wet_cov$CNVQ2x50[dati_corr_wet_cov$Source=="DNA"],
CNVQ2x100=dati_corr_wet_cov$CNVQ2x100[dati_corr_wet_cov$Source=="DNA"],
CNVQ2x250=dati_corr_wet_cov$CNVQ2x250[dati_corr_wet_cov$Source=="DNA"],
CNVQ2x500=dati_corr_wet_cov$CNVQ2x500[dati_corr_wet_cov$Source=="DNA"]
)

library(reshape2)
rdati_DNA<-cor(DNA_corr,method = "kendall"); rdati_DNA
melted_cormat2_DNA <- melt(rdati_DNA); melted_cormat2_DNA

#colnames(melted_cormat2_DNA)<-c("Var2","Var1","value")
head(melted_cormat2_DNA)
melted_cormat2_DNA[,1]<-factor(melted_cormat2_DNA[,1],levels=colnames(rdati_DNA))
head(melted_cormat2_DNA)

library(ggplot2)

#ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + geom_tile()

ggplot(data = melted_cormat2_DNA, aes(x=Var2, y=Var1, fill=value)) + 
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Kendall's Tau") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() + ggtitle("Heatmap of DNA features: Kendall's Tau correlations") +
theme(plot.title = element_text(hjust = 0.5))


#######
##RNA##
#######

colnames(dati_corr_wet_cov)

RNA_corr<-data.frame(RNA=dati_corr_wet_cov$RNA..ng.ul...10ng.ul..8.5ul.[dati_corr_wet_cov$Source=="RNA"],
A260280=dati_corr_wet_cov$RNA.A260.280[dati_corr_wet_cov$Source=="RNA"],
A260230=dati_corr_wet_cov$RNA.A260.230[dati_corr_wet_cov$Source=="RNA"],
DV=dati_corr_wet_cov$X.DV200...20..[dati_corr_wet_cov$Source=="RNA"],
Pre_Hyb=dati_corr_wet_cov$Pre.Hyb.ng.ul....30.[dati_corr_wet_cov$Source=="RNA"],
Post_Hyb=dati_corr_wet_cov$Post.ng.ul....3.[dati_corr_wet_cov$Source=="RNA"],
RNAQ2x5=dati_corr_wet_cov$RNAQ2x5[dati_corr_wet_cov$Source=="RNA"],
RNAQ2x10=dati_corr_wet_cov$RNAQ2x10[dati_corr_wet_cov$Source=="RNA"],
RNAQ2x50=dati_corr_wet_cov$RNAQ2x50[dati_corr_wet_cov$Source=="RNA"]
)

library(reshape2)
rdati_RNA<-cor(RNA_corr,method = "kendall"); rdati_RNA
melted_cormat2_RNA <- melt(rdati_RNA); melted_cormat2_RNA

#colnames(melted_cormat2_RNA)<-c("Var2","Var1","value")
head(melted_cormat2_RNA)
melted_cormat2_RNA[,1]<-factor(melted_cormat2_RNA[,1],levels=colnames(rdati_RNA))
head(melted_cormat2_RNA)

library(ggplot2)

#ggplot(data = melted_cormat2, aes(x=Var2, y=Var1, fill=value)) + geom_tile()

ggplot(data = melted_cormat2_RNA, aes(x=Var2, y=Var1, fill=value)) + 
 geom_tile(color = "white")+
 scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
   midpoint = 0, limit = c(-1,1), space = "Lab", 
   name="Kendall's Tau") +
  theme_minimal()+ 
 theme(axis.text.x = element_text(angle = 45, vjust = 1, 
    size = 12, hjust = 1))+
 coord_fixed() + ggtitle("Heatmap of RNA features: Kendall's Tau correlations") +
theme(plot.title = element_text(hjust = 0.5))

