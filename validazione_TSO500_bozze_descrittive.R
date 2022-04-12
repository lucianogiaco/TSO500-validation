
######################################
#########VALIDAZIONE TSO500###########
######################################


###VALUTARE LE PERFORMANCE DI QUALITA' DEL PANNELLO###


###FARE STATISTICA INCROCIANDO I DATI RELATIVI A:
#-METRICHE DEL SEQUENZIAMENTO
#-COVERAGE
#-VARIABILI CHE CI DEVONO FORNIRE ANGELO ET AL. DALLA GENOMICA RELATIVE:
#   * RUN FATTE A MANO O MEDIANTE HAMILTON (QUESTA POTREBBE ESSERE UNA VARIABILE DI GRUPPO O STRATIFICAZIONE)
#   * TEMPO DI NON REFRIGERAZIONE

#QC: QUALITY CONTROL

#CREARE VARIABILI DICOTOMICHE RISPETTO ALLE LINEE GUIDA:
#1= SODDISFA LE LINEE GUIDA; 0=NON SODDISFA LE LINEE GUIDA (2 TIPI, UNO PER GLI SCORE E UNO PER IL P-VALUE)


#LA VALIDAZIONE (VALUTAZIONE DELLE PERFORMANCE) ANDRA' FATTA A LIVELLO DI:
#-RUN
#-CAMPIONE (SOGGETTO)
#-REGIONE
#-GENE
#-VARIANTI
#...eventualmente incrociando i livelli


###Possibili output:
#-3 radar chart, 1 per livello, in cui si presentano le singole caratteristiche di performance
#-Produrre una black list delle regioni difficili:
#    * Esempio: la precentuale di drop



################
###RADAR PLOT###
################

#DECIDERE LE DIMENSIONI DA PLOTTARE, E SE SONO OSSERVATE O SINTETIZZATE IN QUALCHE "VARIABILE NASCOSTA" (ES. VARIABILE COMPOSITA O LATENTE)

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




###PYRAMID CHART

# Load required R packages
library(dplyr)
library(highcharter) 
# Set highcharter options
options(highcharter.theme = hc_theme_smpl(tooltip = list(valueDecimals = 2)))

df <- data.frame(
        x = c(0, 1, 2, 3, 4, 5),
        y = c(100, 75, 50, 25, 12.5, 6.25),
        name = as.factor(c("run", "subjects", "regions", "genes", "exons", "variants"))
) %>%
  arrange(-y)
df


hc <- df %>%
  hchart(
    "pyramid", hcaes(x = name, y = y),
    name = "Genomic validation levels"
    )



################
####ANALISI#####
################

#install.packages("vioplot")
#install.packages("coin")

library(vioplot)
library(ggplot2)
library(coin)


dati<-read.csv2(file.choose(),sep="\t",h=T,dec=".",na.strings="NA") #importa tab.tab
head(dati)
dim(dati) #63 47
colnames(dati)
summary(dati)


dati[dati$Sample=="37AP",]




###DNA CONTAMINATION SCORES###
colnames(dati)
summary(data.frame(dati$contamination_score,dati$contamination_p_value,dati$Contamination))

par(mfrow=c(1,2))
boxplot(dati$contamination_score)
ggplot(dati, aes(x=factor(rep("",nrow(dati))) , y=contamination_score)) + 
  geom_violin(trim=FALSE, fill="gray")+
  labs(title="Plot of contamination score",x=" ", y = "Contamination score")+
  geom_boxplot(width=0.1)+
  theme_classic()


boxplot(dati$contamination_p_value)

ggplot(ToothGrowth, aes(x=dose, y=len)) + 
  geom_violin(trim=FALSE, fill="gray")+
  labs(title="Plot of contamination p_value",x=" ", y = "Contamination P-value")+
  geom_boxplot(width=0.1)+
  theme_classic()

subset(dati, Contamination==0)

table(dati$Contamination)
table(is.na(dati$Contamination))

tapply(dati$contamination_score,dati$Contamination,summary)
tapply(dati$contamination_score,dati$Contamination,length)

tapply(dati$contamination_p_value,dati$Contamination,summary)
tapply(dati$contamination_p_value,dati$Contamination,length)




par(mfrow=c(1,2))
boxplot(dati$contamination_score, ylab="Contamination score",main="Boxplot of contamination score (all samples)")
vioplot(dati$contamination_score, ylab="Contamination score",names = "",main="Violin plot of contamination score (all samples)")


par(mfrow=c(1,2))
boxplot(dati$contamination_score[dati$Contamination==1], ylab="Contamination score",main="Boxplot of contamination score (No contamination group)")
vioplot(dati$contamination_score[dati$Contamination==1], ylab="Contamination score",names = "",main="Violin plot of contamination score (No contamination group)")




#0 NO CONTAMINAZIONE
#1 CONTAMINAZIONE MA P-VALUE < SOGLIA (LINEEGUIDA) #NO CONTAMINAZIONE MA PROBABILMENTE PRESENZA DI COPY-NUMBER
#1 CONTAMINAZIONE MA P-VALUE > SOGLIA (LINEEGUIDA) #SI CONTAMINAZIONE


###########################
###COPY NUMBER VARIATION###
###########################

#N.b. APPROFONDIRE SOPRATTUTTO LE STATISTICHE SUI COPY NUMBER VARIATION#

summary(dati$contamination_score)
summary(dati$contamination_p_value)
summary(dati$pct_aligned_reads)
summary(dati$pct_read_enrichment)
summary(dati$pct_chimeric_reads)

dati$contamination_score[dati$Contamination==0]
dati$contamination_p_value[dati$Contamination==0]
dati$pct_aligned_reads[dati$Contamination==0]
dati$pct_read_enrichment[dati$Contamination==0]
dati$pct_chimeric_reads[dati$Contamination==0]

length(na.omit(dati$contamination_score[dati$Contamination==1]))

summary(dati$contamination_score[dati$Contamination==1])
summary(dati$contamination_p_value[dati$Contamination==1])
summary(dati$pct_aligned_reads[dati$Contamination==1])
summary(dati$pct_read_enrichment[dati$Contamination==1])
summary(dati$pct_chimeric_reads[dati$Contamination==1])


potent_CNV<-subset(dati, Contamination==1 &  contamination_score >3106); potent_CNV
dim(potent_CNV) #9 47
summary(potent_CNV)
summary(potent_CNV$contamination_score)
summary(potent_CNV$contamination_p_value)
summary(potent_CNV$pct_aligned_reads)
summary(potent_CNV$pct_read_enrichment)
summary(potent_CNV$pct_chimeric_reads)

#write.table(potent_CNV,"potent_CNV.csv",sep=";",col.names=T,row.names=F,dec=",",na="")


###CREO LA VARIABILE BINARIA COPY NUMBER VARIATION (1=CNV, 0=NO CNV)
potCNV<-ifelse((dati$Contamination==1 &  dati$contamination_score >3106 & dati$Contamination!=0),1,0)
potCNV[dati$Contamination==0]<-NA
table(potCNV)
summary(potCNV)

tapply(dati$contamination_score,potCNV,summary)
tapply(dati$contamination_score,potCNV,length)
boxplot(dati$contamination_score~potCNV,ylab="Contamination score", xlab="Potential CNV", names=c("No","Yes"),main="Contamination score") #DA INCLUDERE IN PRESENTAZIONE PPT
text(1.5,65,"P-value=0.881")
wilcox.test(dati$contamination_score~potCNV)

tapply(dati$contamination_p_value,potCNV,summary)
tapply(dati$contamination_p_value,potCNV,length)
boxplot(dati$contamination_p_value~potCNV,ylab="contamination_p_value", xlab="Potential CNV", names=c("No","Yes"),main="contamination_p_value") #DA INCLUDERE IN PRESENTAZIONE PPT
text(1.5,65,"P-value=0.881")
wilcox.test(dati$contamination_p_value~potCNV)



par(mfrow=c(1,2))
tapply(dati$pct_aligned_reads,potCNV,summary)
tapply(dati$pct_aligned_reads,potCNV,length)
boxplot(dati$pct_aligned_reads~potCNV,ylab="Aligned reads (%)", xlab="Potential CNV", names=c("No","Yes"),main="Aligned read by potential CNVs") #DA INCLUDERE IN PRESENTAZIONE PPT
text(1.5,65,"P-value=0.881")
wilcox.test(dati$pct_aligned_reads~potCNV)

tapply(dati$pct_read_enrichment,potCNV,summary)
tapply(dati$pct_read_enrichment,potCNV,length)
boxplot(dati$pct_read_enrichment~potCNV)
wilcox.test(dati$pct_read_enrichment~potCNV)

tapply(dati$pct_chimeric_reads,potCNV,summary)
tapply(dati$pct_chimeric_reads,potCNV,length)
boxplot(dati$pct_chimeric_reads~potCNV,ylab="Chimeric reads (%)", xlab="Potential CNV", names=c("No","Yes"),main="Chimeric read by potential CNVs") #DA INCLUDERE IN PRESENTAZIONE PPT
text(1.5,4.5,"P-value=0.026")
wilcox.test(dati$pct_chimeric_reads~potCNV)




#I POTENZIALI COPY NUMBER VARIATION VANNO INCROCIATI (VALUTARE ASSOCIAZIONE) CON I CAMPIONI IN CUI E' STATA RILEVATA UN COPY NUMBER

summary(potent_CNV$coverage_mad_count)
summary(potent_CNV$median_bin_count_cnv_target_count)
summary(potent_CNV$pct_pf_uq_reads) #pct_pf_uq_reads =qualità delle reads
summary(potent_CNV$pct_aligned_reads)
summary(potent_CNV$pct_aligned_reads)


table(dati$coverage_mad_guidelines)
prop.test(table(dati$coverage_mad_guidelines))

summary(dati$coverage_mad_count)
summary(dati$median_bin_count_cnv_target_count)


table(dati$coverage_mad_guidelines)
table(is.na(dati$coverage_mad_guidelines))




#ESTRARRE I 6 CAMPIONI
dati[dati$coverage_mad_guidelines==0,]
COVERAGE_MAD0<-subset(dati, coverage_mad_guidelines==0); COVERAGE_MAD0
getwd()
write.table(COVERAGE_MAD0,"COVERAGE_MAD0.csv",sep=";",row.names=F,col.names=T)

tapply(dati$coverage_mad_count,dati$coverage_mad_guidelines,summary)
tapply(dati$coverage_mad_count,dati$coverage_mad_guidelines,length)
wilcox.test(dati$coverage_mad_count~dati$coverage_mad_guidelines)
wilcox_test(coverage_mad_count~factor(coverage_mad_guidelines), data = dati,distribution = "exact", conf.int = TRUE)




tapply(dati$median_bin_count_cnv_target_count,dati$coverage_mad_guidelines,summary)
tapply(dati$median_bin_count_cnv_target_count,dati$coverage_mad_guidelines,length)
wilcox.test(dati$median_bin_count_cnv_target_count~dati$coverage_mad_guidelines)
wilcox_test(median_bin_count_cnv_target_count~factor(coverage_mad_guidelines), data = dati,distribution = "exact", conf.int = TRUE)


par(mfrow=c(1,2))
boxplot(dati$coverage_mad_count, ylab="Coverage MAD",main="Boxplot of coverage MAD (all samples)")
vioplot(dati$coverage_mad_count, ylab="coverage MAD",names = "",main="Violin plot of coverage MAD (all samples)")

par(mfrow=c(1,2))
boxplot(dati$coverage_mad_count[dati$coverage_mad_guidelines==1], ylab="Coverage MAD",main="Boxplot of coverage MAD (No group)")
vioplot(dati$coverage_mad_count[dati$coverage_mad_guidelines==1], ylab="Coverage MAD",names = "",main="Violin plot of coverage MAD (No group)")

par(mfrow=c(1,2))
boxplot(dati$coverage_mad_count[dati$coverage_mad_guidelines==0], ylab="Coverage MAD",main="Boxplot of coverage MAD (Yes group)")
vioplot(dati$coverage_mad_count[dati$coverage_mad_guidelines==0], ylab="Coverage MAD",names = "",main="Violin plot of coverage MAD (Yes group)")


par(mfrow=c(1,2))
boxplot(dati$median_bin_count_cnv_target_count, ylab="Median bin CNV target",main="Boxplot of Median bin CNV target (all samples)")
vioplot(dati$median_bin_count_cnv_target_count, ylab="Median bin CNV target",names = "",main="Violin plot of Median bin CNV target (all samples)")

par(mfrow=c(1,2))
boxplot(dati$median_bin_count_cnv_target_count[dati$coverage_mad_guidelines==1], ylab="Median bin CNV target",main="Boxplot of Median bin CNV target (No group)")
vioplot(dati$median_bin_count_cnv_target_count[dati$coverage_mad_guidelines==1], ylab="Median bin CNV target",names = "",main="Violin plot of Median bin CNV target (No group)")

par(mfrow=c(1,2))
boxplot(dati$median_bin_count_cnv_target_count[dati$coverage_mad_guidelines==0], ylab="Median bin CNV target",main="Boxplot of Median bin CNV target (Yes group)")
vioplot(dati$median_bin_count_cnv_target_count[dati$coverage_mad_guidelines==0], ylab="Median bin CNV target",names = "",main="Violin plot of Median bin CNV target (Yes group)")





###DNA METRICS###
colnames(dati)


#TMB
colnames(dati)

summary(dati$median_insert_size)
table(dati$median_insert_size_guidelines)
table(is.na(dati$median_insert_size_guidelines))
#prop.test(table(dati$median_insert_size_guidelines))
#tapply(dati$median_insert_size,dati$median_insert_size_guidelines,summary)
#tapply(dati$median_insert_size,dati$median_insert_size_guidelines,length)
#wilcox.test(dati$median_insert_size~median_insert_size_guidelines)
#wilcox_test(median_insert_size~factor(median_insert_size_guidelines), data = dati,distribution = "exact", conf.int = TRUE)

summary(dati$median_exon_coverage_count)
table(dati$median_exon_coverage_guidelines)
table(is.na(dati$median_exon_coverage_guidelines))
#prop.test(table(dati$median_exon_coverage_guidelines))
tapply(dati$median_exon_coverage_count,dati$median_exon_coverage_guidelines,summary)
tapply(dati$median_exon_coverage_count,dati$median_exon_coverage_guidelines,length)
#wilcox.test(dati$median_exon_coverage_count~dati$median_exon_coverage_guidelines)
#wilcox_test(median_exon_coverage_count~factor(median_exon_coverage_guidelines), data = dati,distribution = "exact", conf.int = TRUE)

summary(dati$pct_exon_50x_perc)
table(dati$pct_exon_50x_guidelines)
table(is.na(dati$coverage_mad_guidelines))
#prop.test(table(dati$pct_exon_50x_guidelines))
tapply(dati$pct_exon_50x_perc,dati$pct_exon_50x_guidelines,summary)
tapply(dati$pct_exon_50x_perc,dati$pct_exon_50x_guidelines,length)
#wilcox.test(dati$pct_exon_50x_perc~dati$pct_exon_50x_guidelines)
#wilcox_test(pct_exon_50x_perc~factor(pct_exon_50x_guidelines), data = dati,distribution = "exact", conf.int = TRUE)




#MSI
colnames(dati)

summary(dati$unstable_msi_sites)
table(dati$unstable_msi_sites_guidelines)
table(dati$unstable_msi_sites_guidelines,dati$unstable_msi_sites>=40)
table(is.na(dati$unstable_msi_sites_guidelines))
prop.test(table(dati$unstable_msi_sites_guidelines))
tapply(dati$unstable_msi_sites,dati$unstable_msi_sites_guidelines,summary)
tapply(dati$unstable_msi_sites,dati$unstable_msi_sites_guidelines,length)
#wilcox.test(dati$unstable_msi_sites~dati$unstable_msi_sites_guidelines)
#wilcox_test(unstable_msi_sites~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)


#summary(dati$mean_family_size)
#tapply(dati$mean_family_size,dati$unstable_msi_sites_guidelines,summary)

summary(mean_target_coverage)
tapply(dati$mean_target_coverage,dati$unstable_msi_sites_guidelines,summary)

summary(median_target_coverage)
tapply(dati$median_target_coverage,dati$unstable_msi_sites_guidelines,summary)

summary(dati$pct_chimeric_reads)
tapply(dati$pct_chimeric_reads,dati$unstable_msi_sites_guidelines,summary)

summary(dati$pct_exon_100x)
tapply(dati$pct_exon_100x,dati$unstable_msi_sites_guidelines,summary)

summary(dati$pct_read_enrichment)
summary(pct_usable_umi_reads)
summary(pct_aligned_reads)

summary(dati$pct_pf_uq_reads)
summary(dati$pct_target_04x_mean)
summary(dati$pct_target_100x)
summary(dati$pct_target_250x)
summary(dati$pct_contamination_est)




###RNA METRICS###
colnames(dati)

rna_median_cv__gene_500x_perc
rna_total_on_target_reads_count
rna_median_insert_size_count
rna_median_cv__gene_500x_guidelines
rna_total_on_target_reads_guidelines
rna_median_insert_size_guidelines
total_pf_reads
rna_pct_chimeric_reads
rna_pct_on_target_reads
rna_scaled_median_gene_coverage
rna_total_pf_reads



###RUN###
colnames(dati)

PCT_PF_READS_perc
PCT_PF_READS_guidelines
PCT_Q30_R1_perc
PCT_Q30_R1_guidelines
PCT_Q30_R2_perc
PCT_Q30_R2_guidelines
run


##########################
###DNA Expanded Metrics### DA METTERE IN PPT
##########################
colnames(dati)

summary(dati$total_pf_reads)
summary(dati$median_target_coverage)
summary(dati$pct_chimeric_reads)
summary(dati$pct_exon_100x)
summary(dati$pct_read_enrichment)
summary(dati$pct_usable_umi_reads)
summary(dati$mean_target_coverage)
summary(dati$pct_aligned_reads)
summary(dati$pct_contamination_est)
summary(dati$pct_pf_uq_reads)
summary(dati$pct_target_04x_mean)
summary(dati$pct_target_100x)
summary(dati$pct_target_250x)



##########################
###RNA Expanded Metrics### DA METTERE IN PPT
##########################
colnames(dati)

summary(dati$rna_total_pf_reads)
summary(dati$rna_scaled_median_gene_coverage)
summary(dati$rna_pct_chimeric_reads)
summary(dati$rna_pct_on_target_reads)




		



#####################################
########ASSOCIATION ANALYSIS#########
#####################################

colnames(dati)

####################################################
###DNA Expanded Metrics with MSI sites_guidelines###
####################################################

#tapply(dati$mean_family_size,dati$unstable_msi_sites_guidelines,summary)
#wilcox.test(dati$mean_family_size~dati$unstable_msi_sites_guidelines)
#wilcox_test(mean_family_size~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)
#boxplot(dati$mean_family_size~dati$unstable_msi_sites_guidelines, main="Mean Family Size",ylab="Mean Family Size", xlab="Unstable MSI sites guidelines")

par(mfrow=c(1,2))
vioplot(dati$mean_family_size[dati$unstable_msi_sites_guidelines==0], ylab="mean family size",names = "unstable msi_sites guidelines==0",main="Violin plot of mean family size")
vioplot(dati$mean_family_size[dati$unstable_msi_sites_guidelines==1], ylab="mean family size",names = "unstable msi_sites guidelines==1",main="Violin plot of mean family size")


tapply(dati$total_pf_reads,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$total_pf_reads~dati$unstable_msi_sites_guidelines)
wilcox_test(total_pf_reads~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)
boxplot(dati$total_pf_reads~dati$unstable_msi_sites_guidelines, main="total_pf_reads",ylab="total_pf_reads", xlab="Unstable MSI sites guidelines")

par(mfrow=c(1,2))
vioplot(dati$total_pf_reads[dati$unstable_msi_sites_guidelines==0], ylab="total_pf_reads",names = "unstable msi_sites guidelines==0",main="Violin plot of total_pf_reads")
vioplot(dati$total_pf_reads[dati$unstable_msi_sites_guidelines==1], ylab="total_pf_reads",names = "unstable msi_sites guidelines==1",main="Violin plot of total_pf_reads")



tapply(dati$median_target_coverage,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$median_target_coverage~dati$unstable_msi_sites_guidelines)
wilcox_test(median_target_coverage~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)
boxplot(dati$median_target_coverage~dati$unstable_msi_sites_guidelines, main="median_target_coverage",ylab="median_target_coverage", xlab="Unstable MSI sites guidelines")

par(mfrow=c(1,2))
vioplot(dati$median_target_coverage[dati$unstable_msi_sites_guidelines==0], ylab="median_target_coverage",names = "unstable msi_sites guidelines==0",main="Violin plot of median_target_coverage")
vioplot(dati$median_target_coverage[dati$unstable_msi_sites_guidelines==1], ylab="median_target_coverage",names = "unstable msi_sites guidelines==1",main="Violin plot of median_target_coverage")


###DA INSERIRE IN WORD (#DA INCLUDERE IN PRESENTAZIONE PPT)
summary(dati$pct_chimeric_reads)
tapply(dati$pct_chimeric_reads,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$pct_chimeric_reads~dati$unstable_msi_sites_guidelines)
wilcox_test(pct_chimeric_reads~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)
boxplot(dati$pct_chimeric_reads~dati$unstable_msi_sites_guidelines, main="pct_chimeric_reads",ylab="pct_chimeric_reads", xlab="Unstable MSI sites guidelines")

par(mfrow=c(1,2))
vioplot(dati$pct_chimeric_reads[dati$unstable_msi_sites_guidelines==0], ylab="pct_chimeric_reads",names = "unstable msi_sites guidelines==0",main="Violin plot of pct_chimeric_reads")
vioplot(dati$pct_chimeric_reads[dati$unstable_msi_sites_guidelines==1], ylab="pct_chimeric_reads",names = "unstable msi_sites guidelines==1",main="Violin plot of pct_chimeric_reads")



tapply(dati$pct_exon_100x,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$pct_exon_100x~dati$unstable_msi_sites_guidelines)
wilcox_test(pct_exon_100x~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)
boxplot(dati$pct_exon_100x~dati$unstable_msi_sites_guidelines, main="pct_exon_100x",ylab="pct_chimeric_reads", xlab="Unstable MSI sites guidelines")

par(mfrow=c(1,2))
vioplot(dati$pct_exon_100x[dati$unstable_msi_sites_guidelines==0], ylab="pct_exon_100x",names = "unstable msi_sites guidelines==0",main="Violin plot of pct_exon_100x")
vioplot(dati$pct_exon_100x[dati$unstable_msi_sites_guidelines==1], ylab="pct_exon_100x",names = "unstable msi_sites guidelines==1",main="Violin plot of pct_exon_100x")



tapply(dati$pct_read_enrichment,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$pct_read_enrichment~dati$unstable_msi_sites_guidelines)
wilcox_test(pct_read_enrichment~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)

tapply(dati$pct_usable_umi_reads,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$pct_usable_umi_reads~dati$unstable_msi_sites_guidelines)
wilcox_test(pct_usable_umi_reads~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)


tapply(dati$mean_target_coverage,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$mean_target_coverage~dati$unstable_msi_sites_guidelines)
wilcox_test(mean_target_coverage~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)


###DA INSERIRE IN WORD (#DA INCLUDERE IN PRESENTAZIONE PPT)
summary(dati$pct_aligned_reads)
tapply(dati$pct_aligned_reads,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$pct_aligned_reads~dati$unstable_msi_sites_guidelines)
wilcox_test(pct_aligned_reads~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)

tapply(dati$pct_contamination_est,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$pct_contamination_est~dati$unstable_msi_sites_guidelines)
wilcox_test(pct_contamination_est~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)


###DA INSERIRE IN WORD (#DA INCLUDERE IN PRESENTAZIONE PPT)
summary(dati$pct_pf_uq_reads)
tapply(dati$pct_pf_uq_reads,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$pct_pf_uq_reads~dati$unstable_msi_sites_guidelines)
wilcox_test(pct_pf_uq_reads~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)

tapply(dati$pct_target_04x_mean,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$pct_target_04x_mean~dati$unstable_msi_sites_guidelines)
wilcox_test(pct_target_04x_mean~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)

tapply(dati$pct_target_100x,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$pct_target_100x~dati$unstable_msi_sites_guidelines)
wilcox_test(pct_target_100x~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)

tapply(dati$pct_target_250x,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$pct_target_250x~dati$unstable_msi_sites_guidelines)
wilcox_test(pct_target_250x~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)


####################################################
###RNA Expanded Metrics with MSI sites_guidelines###
####################################################

summary(dati$rna_total_pf_reads)
summary(dati$rna_scaled_median_gene_coverage)
summary(dati$rna_pct_chimeric_reads)
summary(dati$rna_pct_on_target_reads)


###DA INSERIRE IN WORD (#DA INCLUDERE IN PRESENTAZIONE PPT)
summary(dati$rna_total_pf_reads)
tapply(dati$rna_total_pf_reads,dati$unstable_msi_sites_guidelines,summary)
tapply(dati$rna_total_pf_reads,dati$unstable_msi_sites_guidelines,length)

wilcox.test(dati$rna_total_pf_reads~dati$unstable_msi_sites_guidelines)
wilcox_test(rna_total_pf_reads~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)
boxplot(dati$rna_total_pf_reads~dati$unstable_msi_sites_guidelines, main="rna_total_pf_reads",ylab="rna_total_pf_reads", xlab="Unstable MSI sites guidelines")

par(mfrow=c(1,2))
vioplot(dati$rna_total_pf_reads[dati$unstable_msi_sites_guidelines==0], ylab="rna_total_pf_reads",names = "unstable msi_sites guidelines==0",main="Violin plot of rna_total_pf_reads")
vioplot(dati$rna_total_pf_reads[dati$unstable_msi_sites_guidelines==1], ylab="rna_total_pf_reads",names = "unstable msi_sites guidelines==1",main="Violin plot of rna_total_pf_reads")


###DA INSERIRE IN WORD (#DA INCLUDERE IN PRESENTAZIONE PPT)
summary(dati$rna_scaled_median_gene_coverage)
tapply(dati$rna_scaled_median_gene_coverage,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$rna_scaled_median_gene_coverage~dati$unstable_msi_sites_guidelines)
wilcox_test(rna_scaled_median_gene_coverage~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)
boxplot(dati$rna_scaled_median_gene_coverage~dati$unstable_msi_sites_guidelines, main="rna_scaled_median_gene_coverage",ylab="rna_scaled_median_gene_coverage", xlab="Unstable MSI sites guidelines")

par(mfrow=c(1,2))
vioplot(dati$rna_scaled_median_gene_coverage[dati$unstable_msi_sites_guidelines==0], ylab="rna_scaled_median_gene_coverage",names = "unstable msi_sites guidelines==0",main="Violin plot of rna_scaled_median_gene_coverage")
vioplot(dati$rna_scaled_median_gene_coverage[dati$unstable_msi_sites_guidelines==1], ylab="rna_scaled_median_gene_coverage",names = "unstable msi_sites guidelines==1",main="Violin plot of rna_scaled_median_gene_coverage")


###DA INSERIRE IN WORD (#DA INCLUDERE IN PRESENTAZIONE PPT)
summary(dati$rna_pct_chimeric_reads)
tapply(dati$rna_pct_chimeric_reads,dati$unstable_msi_sites_guidelines,summary)
wilcox.test(dati$rna_pct_chimeric_reads~dati$unstable_msi_sites_guidelines)
wilcox_test(rna_pct_chimeric_reads~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)
boxplot(dati$rna_pct_chimeric_reads~dati$unstable_msi_sites_guidelines, main="rna_pct_chimeric_reads",ylab="rna_scaled_median_gene_coverage", xlab="Unstable MSI sites guidelines")

par(mfrow=c(1,2))
vioplot(dati$rna_pct_chimeric_reads[dati$unstable_msi_sites_guidelines==0], ylab="rna_pct_chimeric_reads",names = "unstable msi_sites guidelines==0",main="Violin plot of rna_pct_chimeric_reads")
vioplot(dati$rna_pct_chimeric_reads[dati$unstable_msi_sites_guidelines==1], ylab="rna_pct_chimeric_reads",names = "unstable msi_sites guidelines==1",main="Violin plot of rna_pct_chimeric_reads")





#############################################
#####CHECK SU NUOVI CAMPIONI DI LUCIANO######
#############################################


dati_new<-read.csv2(file.choose(),sep="\t",h=T,dec=".",na.strings="NA") #importa table check LUCIANO.tab
head(dati_new)
dim(dati_new) # 112  47
colnames(dati_new)
summary(dati_new)


summary(dati$total_pf_reads)
summary(dati$median_target_coverage)
summary(dati$pct_chimeric_reads)
summary(dati$pct_exon_100x)
summary(dati$pct_read_enrichment)
summary(dati$pct_usable_umi_reads)
summary(dati$mean_target_coverage)
summary(dati$pct_aligned_reads)
summary(dati$pct_contamination_est)
summary(dati$pct_pf_uq_reads)
summary(dati$pct_target_04x_mean)
summary(dati$pct_target_100x)
summary(dati$pct_target_250x)

###Check per Luciano###
RA_DNA<-dati_new[dati_new$Sample=="RA_DNA",]
#write.table(RA_DNA,"RA_DNA.csv",sep=";",col.names=T,row.names=F,dec=",",na="")



#####

dati_new<-read.csv2(file.choose(),sep="\t",h=T,dec=".",na.strings="NA") #importa table check table2.tab
head(dati_new)
dim(dati_new) # 112  47
colnames(dati_new)
summary(dati_new)

dati_new$Sample

sample3<-rbind(
dati_new[dati_new$Sample=="MD13",],
dati_new[dati_new$Sample=="89AP_DNA",],
dati_new[dati_new$Sample=="S0705_DNA_FFT",],
dati_new[dati_new$Sample=="101_BM_DNA",],
dati_new[dati_new$Sample=="104_BM_DNA",],
dati_new[dati_new$Sample=="109_BM_DNA",],
dati_new[dati_new$Sample=="58_AP_DNA",],
dati_new[dati_new$Sample=="75_AP_DNA",],
dati_new[dati_new$Sample=="80AP_DNA",],
dati_new[dati_new$Sample=="84AP_DNA",],
dati_new[dati_new$Sample=="87AP_DNA",],
dati_new[dati_new$Sample=="S0733_DNA_FFT",],
dati_new[dati_new$Sample=="S0705_DNA_FFT",],
dati_new[dati_new$Sample=="78AP_DNA",],
dati_new[dati_new$Sample=="89AP_DNA",]
)


dim(sample3) #14 47


getwd()
#write.table(sample3,"14 campioni.csv",sep=";",dec=",",col.names=T,row.names=F,na="")

#Tutti DNA
#Controllo su:
#-total_pf_reads
#-pct_exon_enrichment
#-pct_aligned_reads
#-PCT_PF_UQ_READS (%)
#-PCT_TARGET_0.4X_MEAN
#-PCT_TARGET_100X
#-PCT_TARGET_250X

colnames(sample3)
class(sample3)

par(mfrow=c(2,2))
summary(dati_new$total_pf_reads)
boxplot(dati_new$total_pf_reads,main="Total number of reads passing filter", ylab="units")
legend("bottomleft",legend=c("MD13","500x"),pch=c(20,20),col=c(3,2),cex=0.75,pt.cex=c(1.75,1.75))
points(sample3$total_pf_reads[1],col=3,pch=20,cex=2.3)
points(sample3$total_pf_reads[2],col=2,pch=20,cex=2.3)
points(sample3$total_pf_reads[3],col=2,pch=20,cex=2.3)
points(sample3$total_pf_reads[4],col=2,pch=20,cex=2.3)
points(sample3$total_pf_reads[5],col=2,pch=20,cex=2.3)
points(sample3$total_pf_reads[6],col=2,pch=20,cex=2.3)
points(sample3$total_pf_reads[7],col=2,pch=20,cex=2.3)
points(sample3$total_pf_reads[8],col=2,pch=20,cex=2.3)
points(sample3$total_pf_reads[9],col=2,pch=20,cex=2.3)
points(sample3$total_pf_reads[10],col=2,pch=20,cex=2.3)
points(sample3$total_pf_reads[11],col=2,pch=20,cex=2.3)
points(sample3$total_pf_reads[12],col=2,pch=20,cex=2.3)
points(sample3$total_pf_reads[13],col=2,pch=20,cex=2.3)
points(sample3$total_pf_reads[14],col=2,pch=20,cex=2.3)
mean(sample3$total_pf_reads); sd(sample3$total_pf_reads)
sd(sample3$total_pf_reads)/mean(sample3$total_pf_reads)

summary(dati_new$pct_read_enrichment)
boxplot(dati_new$pct_read_enrichment,main="Read enrichment", ylab="%")
legend("bottomleft",legend=c("MD13","500x"),pch=c(20,20),col=c(3,2),cex=0.75,pt.cex=c(1.75,1.75))
points(sample3$pct_read_enrichment[1],col=3,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[2],col=2,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[3],col=2,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[4],col=2,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[5],col=2,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[6],col=2,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[7],col=2,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[8],col=2,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[9],col=2,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[10],col=2,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[11],col=2,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[12],col=2,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[13],col=2,pch=20,cex=2.3)
points(sample3$pct_read_enrichment[14],col=2,pch=20,cex=2.3)
mean(sample3$pct_read_enrichment); sd(sample3$pct_read_enrichment)
sd(sample3$pct_read_enrichment)/mean(sample3$pct_read_enrichment)

summary(dati_new$pct_aligned_reads)
boxplot(dati_new$pct_aligned_reads,main="Aligned reads", ylab="%")
legend("bottomleft",legend=c("MD13","500x"),pch=c(20,20),col=c(3,2),cex=0.75,pt.cex=c(1.75,1.75))
points(sample3$pct_aligned_reads[1],col=3,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[2],col=2,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[3],col=2,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[4],col=2,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[5],col=2,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[6],col=2,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[7],col=2,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[8],col=2,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[9],col=2,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[10],col=2,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[11],col=2,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[12],col=2,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[13],col=2,pch=20,cex=2.3)
points(sample3$pct_aligned_reads[14],col=2,pch=20,cex=2.3)
mean(sample3$pct_aligned_reads); sd(sample3$pct_aligned_reads)
sd(sample3$pct_aligned_reads)/mean(sample3$pct_aligned_reads)

summary(dati_new$pct_pf_uq_reads)
boxplot(dati_new$pct_pf_uq_reads,main="Unique reads passing filter", ylab="%")
legend("bottomleft",legend=c("MD13","500x"),pch=c(20,20),col=c(3,2),cex=0.75,pt.cex=c(1.75,1.75))
#points(sample3$pct_pf_uq_reads[1],col=2,pch=20,cex=2.3)
#points(sample3$pct_pf_uq_reads[2],col=2,pch=20,cex=2.3)
#points(sample3$pct_pf_uq_reads[3],col=3,pch=20,cex=2.3)
mean(sample3$pct_pf_uq_reads); sd(sample3$pct_pf_uq_reads)
sd(sample3$pct_pf_uq_reads)/mean(sample3$pct_pf_uq_reads)


par(mfrow=c(2,2))
summary(dati_new$pct_target_04x_mean)
boxplot(dati_new$pct_target_04x_mean,main="Target 04x mean", ylab="%")
legend("bottomleft",legend=c("MD13","500x"),pch=c(20,20),col=c(3,2),cex=0.75,pt.cex=c(1.75,1.75))
points(sample3$pct_target_04x_mean[1],col=3,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[2],col=2,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[3],col=2,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[4],col=2,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[5],col=2,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[6],col=2,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[7],col=2,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[8],col=2,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[9],col=2,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[10],col=2,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[11],col=2,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[12],col=2,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[13],col=2,pch=20,cex=2.3)
points(sample3$pct_target_04x_mean[14],col=2,pch=20,cex=2.3)
mean(sample3$pct_target_04x_mean); sd(sample3$pct_target_04x_mean)
sd(sample3$pct_target_04x_mean)/mean(sample3$pct_target_04x_mean)

summary(dati_new$pct_target_100x)
boxplot(dati_new$pct_target_100x,main="Target 100x", ylab="%")
legend("bottomleft",legend=c("MD13","500x"),pch=c(20,20),col=c(3,2),cex=0.75,pt.cex=c(1.75,1.75))
points(sample3$pct_target_100x[1],col=3,pch=20,cex=2.3)
points(sample3$pct_target_100x[2],col=2,pch=20,cex=2.3)
points(sample3$pct_target_100x[3],col=2,pch=20,cex=2.3)
points(sample3$pct_target_100x[4],col=2,pch=20,cex=2.3)
points(sample3$pct_target_100x[5],col=2,pch=20,cex=2.3)
points(sample3$pct_target_100x[6],col=2,pch=20,cex=2.3)
points(sample3$pct_target_100x[7],col=2,pch=20,cex=2.3)
points(sample3$pct_target_100x[8],col=2,pch=20,cex=2.3)
points(sample3$pct_target_100x[9],col=2,pch=20,cex=2.3)
points(sample3$pct_target_100x[10],col=2,pch=20,cex=2.3)
points(sample3$pct_target_100x[11],col=2,pch=20,cex=2.3)
points(sample3$pct_target_100x[12],col=2,pch=20,cex=2.3)
points(sample3$pct_target_100x[13],col=2,pch=20,cex=2.3)
points(sample3$pct_target_100x[14],col=2,pch=20,cex=2.3)
mean(sample3$pct_target_100x); sd(sample3$pct_target_100x)
sd(sample3$pct_target_100x)/mean(sample3$pct_target_100x)

summary(dati_new$pct_target_250x)
boxplot(dati_new$pct_target_250x,main="Target 250x", ylab="%")
legend("bottomleft",legend=c("MD13","500x"),pch=c(20,20),col=c(3,2),cex=0.75,pt.cex=c(1.75,1.75))
points(sample3$pct_target_250x[1],col=3,pch=20,cex=2.3)
points(sample3$pct_target_250x[2],col=2,pch=20,cex=2.3)
points(sample3$pct_target_250x[3],col=2,pch=20,cex=2.3)
points(sample3$pct_target_250x[4],col=2,pch=20,cex=2.3)
points(sample3$pct_target_250x[5],col=2,pch=20,cex=2.3)
points(sample3$pct_target_250x[6],col=2,pch=20,cex=2.3)
points(sample3$pct_target_250x[7],col=2,pch=20,cex=2.3)
points(sample3$pct_target_250x[8],col=2,pch=20,cex=2.3)
points(sample3$pct_target_250x[9],col=2,pch=20,cex=2.3)
points(sample3$pct_target_250x[10],col=2,pch=20,cex=2.3)
points(sample3$pct_target_250x[11],col=2,pch=20,cex=2.3)
points(sample3$pct_target_250x[12],col=2,pch=20,cex=2.3)
points(sample3$pct_target_250x[13],col=2,pch=20,cex=2.3)
points(sample3$pct_target_250x[14],col=2,pch=20,cex=2.3)
mean(sample3$pct_target_250x); sd(sample3$pct_target_250x)
sd(sample3$pct_target_250x)/mean(sample3$pct_target_250x)



#######################
####ANALISI ESTESA#####
#######################

#install.packages("vioplot")
#install.packages("coin")

library(vioplot)
library(ggplot2)
library(coin)


dati<-read.csv2(file.choose(),sep="\t",h=T,dec=".",na.strings="NA") #importa tab2.tab
head(dati)
dim(dati) #112 47
colnames(dati)
summary(dati)
table(dati$run) #numero di run: 5

summary(dati$rna_median_cv__gene_500x_perc)
dati$rna_median_cv__gene_500x_perc[dati$rna_median_cv__gene_500x_perc==117.08]<-100
summary(dati$rna_median_cv__gene_500x_perc)

###verificare i potential CNVs valutando i contamination score, p-value e reads chimeriche



########################################################################
###Valutare i campioni solo DNA, solo RNA, e con entrambe le features###
########################################################################

colnames(dati)

#########
###DNA###
#########

datiDNAsample<-data.frame(dati$contamination_score,
dati$contamination_p_value,
dati$total_pf_reads,
dati$median_target_coverage,
dati$pct_chimeric_reads,
dati$pct_exon_100x,
dati$pct_read_enrichment,
dati$pct_usable_umi_reads,
dati$mean_target_coverage,
dati$pct_aligned_reads,
dati$pct_contamination_est,
dati$pct_pf_uq_reads,
dati$pct_target_04x_mean,
dati$pct_target_100x,
dati$pct_target_250x,
dati$median_insert_size,
dati$median_exon_coverage_count,
dati$pct_exon_50x_perc)

dim(datiDNAsample) #112 18 
summary(datiDNAsample)


####################################################
###DNA Expanded Metrics with MSI sites_guidelines###
####################################################


###DA INSERIRE IN WORD (#DA INCLUDERE IN PRESENTAZIONE PPT)
summary(dati$pct_chimeric_reads)
tapply(dati$pct_chimeric_reads,dati$unstable_msi_sites_guidelines,summary)
tapply(dati$pct_chimeric_reads,dati$unstable_msi_sites_guidelines,length)
wilcox.test(dati$pct_chimeric_reads~dati$unstable_msi_sites_guidelines, conf.int = TRUE)
wilcox_test(pct_chimeric_reads~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)
boxplot(dati$pct_chimeric_reads~dati$unstable_msi_sites_guidelines, main="pct_chimeric_reads",ylab="pct_chimeric_reads", xlab="Unstable MSI sites guidelines")

par(mfrow=c(1,2))
vioplot(dati$pct_chimeric_reads[dati$unstable_msi_sites_guidelines==0], ylab="pct_chimeric_reads",names = "unstable msi_sites guidelines==0",main="Violin plot of pct_chimeric_reads")
vioplot(dati$pct_chimeric_reads[dati$unstable_msi_sites_guidelines==1], ylab="pct_chimeric_reads",names = "unstable msi_sites guidelines==1",main="Violin plot of pct_chimeric_reads")



###DA INSERIRE IN WORD (#DA INCLUDERE IN PRESENTAZIONE PPT)
summary(dati$pct_aligned_reads)
tapply(dati$pct_aligned_reads,dati$unstable_msi_sites_guidelines,summary)
tapply(dati$pct_aligned_reads,dati$unstable_msi_sites_guidelines,length)
wilcox.test(dati$pct_aligned_reads~dati$unstable_msi_sites_guidelines)
wilcox_test(pct_aligned_reads~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)


###DA INSERIRE IN WORD (#DA INCLUDERE IN PRESENTAZIONE PPT)
summary(dati$pct_pf_uq_reads)
tapply(dati$pct_pf_uq_reads,dati$unstable_msi_sites_guidelines,summary)
tapply(dati$pct_pf_uq_reads,dati$unstable_msi_sites_guidelines,length)
wilcox.test(dati$pct_pf_uq_reads~dati$unstable_msi_sites_guidelines)
wilcox_test(pct_pf_uq_reads~factor(unstable_msi_sites_guidelines), data = dati,distribution = "exact", conf.int = TRUE)



###MATRICE MISSING DNA###

colnames(dati)

#https://cran.r-project.org/web/packages/naniar/vignettes/naniar-visualisation.html
install.packages("naniar")
library(naniar)

#at-a-glance ggplot of the missingness
dati[,2:46]
vis_miss(dati[,2:46])

NA_datiDNAsample<-apply(datiDNAsample,1,function(x) 1*is.element(NA,x)); NA_datiDNAsample


#########
###RNA###
#########

datiRNAsample<-data.frame(dati$rna_total_pf_reads,
dati$rna_scaled_median_gene_coverage,
dati$rna_pct_chimeric_reads,
dati$rna_pct_on_target_reads,
dati$rna_median_cv__gene_500x_perc)

dim(datiRNAsample) #112 5
summary(datiRNAsample)





###MATRICE MISSING RNA###
NA_datiRNAsample<-apply(datiRNAsample,1,function(x) 1*is.element(NA,x)); NA_datiRNAsample


###PLOT MISSING###
vis_miss(data.frame(datiDNAsample,datiRNAsample))


NA_both<-ifelse(NA_datiDNAsample==0 & NA_datiRNAsample==0,"DNA & RNA sample",NA)
NA_both<-ifelse(NA_datiDNAsample==0 & NA_datiRNAsample==1,"DNA sample",NA_both)
NA_both<-ifelse(NA_datiDNAsample==1 & NA_datiRNAsample==0,"RNA sample",NA_both)
NA_both<-ifelse(NA_datiDNAsample==1 & NA_datiRNAsample==1,"Failed sample",NA_both)

dati[NA_both=="Failed sample",]


NA_both<-factor(NA_both,levels=c("DNA sample","RNA sample","DNA & RNA sample","Failed sample"))

counts_ass<-table(NA_both); counts_ass
sum(table(NA_both))

counts_perc<-counts_ass/sum(counts_ass)*100; counts_perc
sum(counts_perc)


par(mfrow=c(1,2))
bp_ass<-barplot(counts_ass, main="Absolute frequencies of the samples tipology",
xlab="Sample tipology", col=c("blue","red","forestgreen","black"),ylab="Absolute frequency", beside=TRUE, ylim=c(0,60),cex.main=1.3,cex.names=1.5,cex.lab=1.5)
#legend("topright",legend=c("DNA sample","RNA sample","DNA & RNA sample","Failed sample"),cex=0.75,fill=c("blue","red","forestgreen","black"))
text(x = bp_ass, y = as.numeric(counts_ass), label = round(as.numeric(counts_ass),1), pos = 3, cex = 1.8, col = c("blue","red","forestgreen","black"))

bp_perc<-barplot(counts_perc, main="Percentage frequencies of the sample tipology", beside=TRUE, ylim=c(0,100),cex.main=1.3,
xlab="Sample tipology", col=c("blue","red","forestgreen","black"),ylab="Percent frequency",cex.names=1.5,cex.lab=1.5)
#legend("topright",legend = c("DNA sample","RNA sample","DNA & RNA sample","Failed sample"),cex=0.75,fill=c("blue","red","forestgreen","black"))
text(x = bp_perc, y = as.numeric(counts_perc), label = round(as.numeric(counts_perc),1), pos = 3, cex = 1.8, col = c("blue","red","forestgreen","black"))



 

#########################################################
#####Radar plot per le features in giallo di DNA#########
#########################################################

## The default radar chart 
#par(mar=c(2, 0, 2, 0))
#radarchart(data, vlabels=c("Dimension 1\n(from ...)", 
#"Dimension 2\n(from ...)", 
#"Dimension 3\n(from ...)", 
#"Dimension 4\n(from ...)"),vlcex=0.65,axistype=4,axislabcol="black",calcex=0.7)
#title("Radar chart of the performance \n(subject level)",cex.main=0.85)


summary(dati$pct_chimeric_reads)
summary(dati$pct_exon_100x)
summary(dati$pct_read_enrichment)
summary(dati$pct_aligned_reads)
summary(dati$pct_target_100x)
summary(dati$pct_target_250x)


length(na.omit(dati$pct_chimeric_reads))
length(na.omit(dati$pct_exon_100x))
length(na.omit(dati$pct_read_enrichment))
length(na.omit(dati$pct_aligned_reads))
length(na.omit(dati$pct_target_100x))
length(na.omit(dati$pct_target_250x))



inf_DNA<-rep(0,6); inf_DNA
sup_DNA<-rep(100,6); sup_DNA


med_DNA<-c(median(dati$pct_chimeric_reads,na.rm=T),median(dati$pct_exon_100x,na.rm=T),median(dati$pct_read_enrichment,na.rm=T),
           median(dati$pct_aligned_reads,na.rm=T),median(dati$pct_target_100x,na.rm=T),median(dati$pct_target_250x,na.rm=T)); med_DNA

min_DNA<-c(min(dati$pct_chimeric_reads,na.rm=T),min(dati$pct_exon_100x,na.rm=T),min(dati$pct_read_enrichment,na.rm=T),
           min(dati$pct_aligned_reads,na.rm=T),min(dati$pct_target_100x,na.rm=T),min(dati$pct_target_250x,na.rm=T)); min_DNA

max_DNA<-c(max(dati$pct_chimeric_reads,na.rm=T),max(dati$pct_exon_100x,na.rm=T),max(dati$pct_read_enrichment,na.rm=T),
           max(dati$pct_aligned_reads,na.rm=T),max(dati$pct_target_100x,na.rm=T),max(dati$pct_target_250x,na.rm=T)); max_DNA


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
title("Radar chart of the DNA performance \n(sample level, n=71)",cex.main=1.25)
legend("bottomleft",lty=c(1,2,3),col=c(1,2,3),lwd=c(2,2,2),c("min","median","max"),cex=1)




#############################################################################################################
#####Radar plot per le features in giallo di RNA (con anche la feature RNA median cv gene 500x perc)#########
#############################################################################################################

summary(dati$rna_pct_chimeric_reads)
summary(dati$rna_median_cv__gene_500x_perc)
summary(dati$rna_pct_on_target_reads)

length(na.omit(dati$rna_pct_chimeric_reads))
length(na.omit(dati$rna_median_cv__gene_500x_perc))
length(na.omit(dati$rna_pct_on_target_reads))


inf_RNA<-rep(0,3); inf_RNA
sup_RNA<-rep(100,3); sup_RNA

med_RNA<-c(median(dati$rna_pct_chimeric_reads,na.rm=T),median(dati$rna_median_cv__gene_500x_perc,na.rm=T),median(dati$rna_pct_on_target_reads,na.rm=T))
med_RNA


min_RNA<-c(min(dati$rna_pct_chimeric_reads,na.rm=T),min(dati$rna_median_cv__gene_500x_perc,na.rm=T),min(dati$rna_pct_on_target_reads,na.rm=T))
min_RNA


max_RNA<-c(max(dati$rna_pct_chimeric_reads,na.rm=T),max(dati$rna_median_cv__gene_500x_perc,na.rm=T),max(dati$rna_pct_on_target_reads,na.rm=T))
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
title("Radar chart of the RNA performance \n(sample level, n=59)",cex.main=1.25)
legend("bottomleft",lty=c(1,2,3),col=c(1,2,3),lwd=c(2,2,2),c("min","median","max"),cex=1)





###########################
###COPY NUMBER VARIATION###
###########################

#N.b. APPROFONDIRE SOPRATTUTTO LE STATISTICHE SUI COPY NUMBER VARIATION#

summary(dati$contamination_score)
summary(dati$contamination_p_value)
summary(dati$pct_aligned_reads)
summary(dati$pct_read_enrichment)
summary(dati$pct_chimeric_reads)

length(na.omit(dati$contamination_score[dati$Contamination==0]))

summary(dati$contamination_score[dati$Contamination==0])
summary(dati$contamination_p_value[dati$Contamination==0])
summary(dati$pct_aligned_reads[dati$Contamination==0])
summary(dati$pct_read_enrichment[dati$Contamination==0])
summary(dati$pct_chimeric_reads[dati$Contamination==0])

length(na.omit(dati$contamination_score[dati$Contamination==1]))

summary(dati$contamination_score[dati$Contamination==1])
summary(dati$contamination_p_value[dati$Contamination==1])
summary(dati$pct_aligned_reads[dati$Contamination==1])
summary(dati$pct_read_enrichment[dati$Contamination==1])
summary(dati$pct_chimeric_reads[dati$Contamination==1])


potent_CNV<-subset(dati, Contamination==1 &  contamination_score >3106); potent_CNV
dim(potent_CNV) #15 47
summary(potent_CNV)
summary(potent_CNV$contamination_score)
summary(potent_CNV$contamination_p_value)
summary(potent_CNV$pct_aligned_reads)
summary(potent_CNV$pct_read_enrichment)
summary(potent_CNV$pct_chimeric_reads)

#write.table(potent_CNV,"potent_CNV.csv",sep=";",col.names=T,row.names=F,dec=",",na="")


###CREO LA VARIABILE BINARIA COPY NUMBER VARIATION (1=CNV, 0=NO CNV)
potCNV<-ifelse((dati$Contamination==1 &  dati$contamination_score >3106 & dati$Contamination!=0),1,0)
potCNV[dati$Contamination==0]<-NA
table(potCNV)
summary(potCNV)

tapply(dati$contamination_score,potCNV,summary)
tapply(dati$contamination_score,potCNV,length)
boxplot(dati$contamination_score~potCNV,ylab="Contamination score", xlab="Potential CNV", names=c("No","Yes"),main="Contamination score") #DA INCLUDERE IN PRESENTAZIONE PPT
text(1.5,65,"P-value<0.001")
wilcox.test(dati$contamination_score~potCNV)

tapply(dati$contamination_p_value,potCNV,summary)
tapply(dati$contamination_p_value,potCNV,length)
boxplot(dati$contamination_p_value~potCNV,ylab="contamination_p_value", xlab="Potential CNV", names=c("No","Yes"),main="contamination_p_value") #DA INCLUDERE IN PRESENTAZIONE PPT
text(1.5,65,"P-value=0.881")
wilcox.test(dati$contamination_p_value~potCNV)



par(mfrow=c(1,2))
tapply(dati$pct_aligned_reads,potCNV,summary)
tapply(dati$pct_aligned_reads,potCNV,length)
boxplot(dati$pct_aligned_reads~potCNV,ylab="Aligned reads (%)", xlab="Potential CNV", names=c("No","Yes"),main="Aligned read by potential CNVs") #DA INCLUDERE IN PRESENTAZIONE PPT
text(1.5,65,"P-value=0.813")
wilcox.test(dati$pct_aligned_reads~potCNV)

tapply(dati$pct_read_enrichment,potCNV,summary)
tapply(dati$pct_read_enrichment,potCNV,length)
boxplot(dati$pct_read_enrichment~potCNV)
wilcox.test(dati$pct_read_enrichment~potCNV)

tapply(dati$pct_chimeric_reads,potCNV,summary)
tapply(dati$pct_chimeric_reads,potCNV,length)
boxplot(dati$pct_chimeric_reads~potCNV,ylab="Chimeric reads (%)", xlab="Potential CNV", names=c("No","Yes"),main="Chimeric read by potential CNVs") #DA INCLUDERE IN PRESENTAZIONE PPT
text(1.5,4.5,"P-value=0.078")
wilcox.test(dati$pct_chimeric_reads~potCNV)




#I POTENZIALI COPY NUMBER VARIATION VANNO INCROCIATI (VALUTARE ASSOCIAZIONE) CON I CAMPIONI IN CUI E' STATA RILEVATA UN COPY NUMBER

summary(potent_CNV$coverage_mad_count)
summary(potent_CNV$median_bin_count_cnv_target_count)
summary(potent_CNV$pct_pf_uq_reads) #pct_pf_uq_reads =qualità delle reads
summary(potent_CNV$pct_aligned_reads)
summary(potent_CNV$pct_aligned_reads)


table(dati$coverage_mad_guidelines)
prop.test(table(dati$coverage_mad_guidelines))

summary(dati$coverage_mad_count)
summary(dati$median_bin_count_cnv_target_count)


table(dati$coverage_mad_guidelines)
table(is.na(dati$coverage_mad_guidelines))







#####################################################################
###FARE DESCRITTIVE DNA PER GRUPPI DI SOGGETTI CREATI DA FERNANDO####
#####################################################################


dati_cov_pat<-read.csv2(file.choose(),h=T,sep=";",dec=",",na.strings="") #importo TSO500_SNV_coverage_byPatient.csv
dim(dati_cov_pat) #72 26
head(dati_cov_pat)
colnames(dati_cov_pat)
summary(dati_cov_pat)


###########################
#######Agreement###########
###########################

########################
###Variabili continue###
########################

#DA CONSIDERARE IL Q1 PER OGNI LIVELLO DI PROFONDITA' (DEPT)#

data_long<-data.frame(ID=rep(dati_cov_pat$PatientID,4),prof=c(rep(50,72),rep(100,72),rep(250,72),rep(500,72)),
Q1_value=c(dati_cov_pat$Q1x50,dati_cov_pat$Q1x100,dati_cov_pat$Q1x250,dati_cov_pat$Q1x500),
Q2_value=c(dati_cov_pat$Q2x50,dati_cov_pat$Q2x100,dati_cov_pat$Q2x250,dati_cov_pat$Q2x500))
dim(data_long) #288   4
head(data_long)
colnames(data_long)
summary(data_long)
table(data_long$ID)

by(data_long$Q1_value,factor(data_long$prof),summary) #mediane per valore di Q1
by(data_long$Q2_value,factor(data_long$prof),summary) #mediane per valore di Q2

par(mfrow=c(1,2))
plot(factor(data_long$prof), data_long$Q1_value,xlab="Depth (read count)",ylab="Q1 Exon Coverage (%)", main="Boxplots of Q1 Exon Coverage by dept",names=c("50X","100X","250X","500X"))
abline(h=75,col=2,lty=2)
plot(factor(data_long$prof), data_long$Q2_value,xlab="Depth (read count)",ylab="Median Exon Coverage (%)", main="Boxplots of median Exon Coverage by dept",names=c("50X","100X","250X","500X"))
abline(h=75,col=2,lty=2)

#75% è il limite minimo di coverage consentito



#########################
####BLAND-ALTMAN PLOT####
#########################

colnames(dati_cov_pat)

x1=dati_cov_pat$Q1x50
x2=dati_cov_pat$Q1x100
x3=dati_cov_pat$Q1x250
x4=dati_cov_pat$Q1x500

d12=x1-x2 #differenza (asse y del plot di bland-altman)
mu12=(x1+x2)/2; mu12 #medie per coppia di punti, i.e. media tra strumenti (asse x del plot di bland altman)
#(x1[1]+x2[1])/2

plot((x1+x2)/2, x1-x2, main="Mean-Difference-Plot")


#install.packages("BlandAltmanLeh")
library(BlandAltmanLeh)

#bland.altman.plot(x1, x2, x3, x4, main="This is a Bland Altman Plot", xlab="Means", ylab="Differences")

bland.altman.plot(x1, x2, main="This is a Bland Altman Plot", xlab="Means", ylab="Differences",silent=FALSE)
bland.altman.plot(x2, x3, main="This is a Bland Altman Plot", xlab="Means", ylab="Differences",silent=FALSE)
bland.altman.plot(x3, x4, main="This is a Bland Altman Plot", xlab="Means", ylab="Differences",silent=FALSE)

par(mfrow=c(1,3))
bland.altman.plot(x1, x2, main="This is a Bland Altman Plot", xlab="Means", ylab="Differences",silent=FALSE)
bland.altman.plot(x1, x3, main="This is a Bland Altman Plot", xlab="Means", ylab="Differences",silent=FALSE)
bland.altman.plot(x1, x4, main="This is a Bland Altman Plot", xlab="Means", ylab="Differences",silent=FALSE)

?bland.altman.stats
ba.stats <- bland.altman.stats(x1, x2); ba.stats
ba.stats <- bland.altman.stats(x2, x3); ba.stats
ba.stats <- bland.altman.stats(x3, x4); ba.stats


#Bland-Altman limits by hand
MEANd12=mean(d12); MEANd12
SDd12=sd(d12); SDd12
lower=MEANd12-1.96*sd(d12); lower
upper=MEANd12+1.96*sd(d12); upper

#usando il quantile della T-student (anziché della normale)
qt(c(.025, .975), 98)
lowerT=MEANd12-1.984467*sd(d12); lowerT
upperT=MEANd12+1.984467*sd(d12); upperT


library(ggplot2)
pl <- bland.altman.plot(x1, x2, graph.sys = "ggplot2")
pl <- bland.altman.plot(x2, x3, graph.sys = "ggplot2")
pl <- bland.altman.plot(x3, x4, graph.sys = "ggplot2")

print(pl)
str(pl)
#calcolando i limiti di confidenza sui limiti di Bland-Altman
bland.altman.plot(x1, x2, conf.int=.95, pch=19)



#################################################
####Coefficiente ICC (Intra Class Correlation)###
#################################################

#install.packages("irr")
library(irr)

#Inter-rater reliability
http://www.cookbook-r.com/Statistical_analysis/Inter-rater_reliability/ (Cookbook for R)

#esempio
#r1 <- round(rnorm(20, 10, 4))
#r2 <- round(r1 + 10 + rnorm(20, 0, 2))
#r3 <- round(r1 + 20 + rnorm(20, 0, 2))
#icc(cbind(r1, r2, r3), "twoway")              # High consistency
#icc(cbind(r1, r2, r3), "twoway", "agreement") # Low agreement


x1=dati_cov_pat$Q1x50 
x2=dati_cov_pat$Q1x100
x3=dati_cov_pat$Q1x250
x4=dati_cov_pat$Q1x500

X<-data.frame(x1,x2)
icc(X, model="twoway", type="agreement")

X<-data.frame(x1,x3)
icc(X, model="twoway", type="agreement")

X<-data.frame(x1,x4)
icc(X, model="twoway", type="agreement")

X<-data.frame(x1,x2,x3,x4)
icc(X, model="twoway", type="agreement")



##########################################
###Variabili categoriche (tricotomiche)###
##########################################


#############################
#######Kappa di Cohen########
#############################

#Info su Cohen's Kappa
#http://www.statisticshowto.com/cohens-kappa-statistic/

library(psych) #per la kappa (e 95%CI) e per la phi
library(boot) #per la kappa di cohen bootstrap

?cohen.kappa # si trova anche il riferimento per l'uso della phi


####################
###Weighted Kappa###
####################

#https://www.datanovia.com/en/lessons/weighted-kappa-in-r-for-two-ordinal-variables/

#install.packages("vcd")
library(vcd)

colnames(dati_cov_pat)
head(dati_cov_pat)

x1<-dati_cov_pat$FLAGx50
x2<-dati_cov_pat$FLAGx100
x3<-dati_cov_pat$FLAGx250
x4<-dati_cov_pat$FLAGx500


###x1 vs x2###

X<-table(x1,x2); X
X1 <- as.table(
  rbind(
    c(67, 4, 0), c(0, 0, 0),
    c(0, 0, 1))
)
dimnames(X1) <- list(
  x1 = c("0", "1", "2"),
  x2 = c("0", "1", "2")
)
X1

res.k <- Kappa(X1); res.k

# Confidence intervals
confint(res.k)

# Summary showing the weights assigned to each cell
summary(res.k)


###x1 vs x3###

X<-table(x1,x3); X
X2 <- as.table(
  rbind(
    c(55, 13, 3), c(0, 0, 0),
    c(0, 0, 1))
)
dimnames(X2) <- list(
  x1 = c("0", "1", "2"),
  x3 = c("0", "1", "2")
)
X2

res.k <- Kappa(X2); res.k

# Confidence intervals
confint(res.k)

# Summary showing the weights assigned to each cell
summary(res.k)



###x1 vs x4###

X<-table(x1,x4); X
X3 <- as.table(
  rbind(
    c(35, 23, 13), c(0, 0, 0),
    c(0, 0, 1))
)
dimnames(X3) <- list(
  x1 = c("0", "1", "2"),
  x4 = c("0", "1", "2")
)
X3

res.k <- Kappa(X3); res.k

# Confidence intervals
confint(res.k)

# Summary showing the weights assigned to each cell
summary(res.k)


https://www.datanovia.com/en/lessons/inter-rater-agreement-chart-in-r/
par(mar = c(4, 2, 2, 2))
library(vcd)
# Create the plot
p <- agreementplot(X1)
p <- agreementplot(X2)
p <- agreementplot(X3)




##########################
###95%CI AGREEMENT PLOT###
##########################

library(ggplot2)

pippo<-data.frame(agree=c(0.6,0.5,0.4),min=c(0.4,0.3,0.2),max=c(0.8,0.7,0.6),depth=c("100","250","500"))
pippo

pd = position_dodge(.2)    ### How much to jitter the points on the plot

ggplot(pippo,                ### The data frame to use.
       aes(x     = depth,
           y     = agree)) +

    geom_point(shape = 15,
               size  = 4,
             position = pd) +

    geom_errorbar(aes(ymin  = agree - min,
                      ymax  = agree + max),
                      width = 0.2,
                      size  = 0.7,
                      position = pd) +
    theme_bw() +
    theme(axis.title = element_text(face = "bold")) +

    ylab("EMPC") +

   labs(title = "Agreement coefficient", subtitle = "zhao")



##################
###Fleiss Kappa###
##################

#https://www.datanovia.com/en/lessons/fleiss-kappa-in-r-for-multiple-categorical-variables/

x1<-dati_cov_pat$FLAGx50
x2<-dati_cov_pat$FLAGx100
x3<-dati_cov_pat$FLAGx250
x4<-dati_cov_pat$FLAGx500

X<-data.frame(x1,x2,x3); X
kappam.fleiss(X)
kappam.fleiss(X, detail = TRUE)

X<-data.frame(x2,x3,x4); X
kappam.fleiss(X)
kappam.fleiss(X, detail = TRUE)

X<-data.frame(x1,x2,x3,x4); X
kappam.fleiss(X)
kappam.fleiss(X, detail = TRUE)


##########################
###95%CI AGREEMENT PLOT###
##########################

library(ggplot2)

pippo<-data.frame(agree=c(0.6,0.5,0.4),min=c(0.4,0.3,0.2),max=c(0.8,0.7,0.6),depth=c("100","250","500"))
pippo

pd = position_dodge(.2)    ### How much to jitter the points on the plot

ggplot(pippo,                ### The data frame to use.
       aes(x     = depth,
           y     = agree)) +

    geom_point(shape = 15,
               size  = 4,
             position = pd) +

    geom_errorbar(aes(ymin  = agree - min,
                      ymax  = agree + max),
                      width = 0.2,
                      size  = 0.7,
                      position = pd) +
    theme_bw() +
    theme(axis.title = element_text(face = "bold")) +

    ylab("EMPC") +

   labs(title = "Agreement coefficient", subtitle = "zhao")


###########################################################
#####METRICHE DNA-EXTENDED SUI CAMPIONI CHE FALLISCONO#####
###########################################################

rbind(
dati[dati$Sample=="NR_RNA",],
dati[dati$Sample=="CTRL_DIL_DNA",])





#############################################################
######FARE ANALISI STRATIFICATE PER RUN E PER HAMILTON#######
#############################################################

#MERGIARE PER ID SAMPLE E RUN CON IL DATASET summary_samples.csv#

dati_sum_sample<-read.csv2(file.choose(),h=T,sep=";",dec=",",na.strings="") #importo summary_samples.csv
dim(dati_sum_sample) #132  20
head(dati_sum_sample)
colnames(dati_sum_sample)
summary(dati_sum_sample)

colnames(dati)

dati_all<-merge(dati,dati_sum_sample,by.x=c("Sample","run"),by.y=c("Sample_ID","Run"),all=F)
dim(dati_all) # 90 65
head(dati_all)
colnames(dati_all)
summary(dati_all)




#############################
#####RECOVERY ANALYSIS#######
#############################


#ANALISI CON INCOGNITE DEPTH: IN UN CAMPIONE CHE E' FALLITO, A CHE PROFONDITA' POSSO DIRE CHE UNA MUTAZIONE PUO' ESSERE CONSIDERATA ?



########################
###DETERMINISTIC PLOT###
########################

f1 <- function(t) 20 * (1-exp(-0.001 * -2*t))+100
abline(h=75,col=2,lty=2) #threshold di validità per la dept (da leggere sulle x)
curve(f1, from = 0, to = 1000, xlab="Depth (X)",ylab="(Exon or Gene) Coverage [EMPC-Exon Median Percent Coverage]", ylim=c(0,100))
title("Coverage by Depth")



###############################
####DOSE RESPONSE MODELLING####
###############################




###########################################################################################################
######DOMANDA DI RICERCA: DISTINZIONE DI UN CAMPIONE (TYPE) BUONO (1) DA UNO CATTIVO (0): COME FARE?####### 
###########################################################################################################

#Dose-Response Analysis Using R (2020), by Christian Ritz et al.

#install.packages("drc")
library(drc)
library(sandwich)
library(drcData)
library(boot)
library(lmtest)
library(metafor)


EMPC50_0<-runif(100,85,95)
EMPC50_1<-runif(100,90,100)

EMPC100_0<-runif(100,75,95)
EMPC100_1<-runif(100,80,100)

EMPC250_0<-runif(100,60,80)
EMPC250_1<-runif(100,70,90)

EMPC500_0<-runif(100,10,50)
EMPC500_1<-runif(100,20,60)


EMPC_0<-c(EMPC50_0,EMPC100_0,EMPC250_0,EMPC500_0)
EMPC_1<-c(EMPC50_1,EMPC100_1,EMPC250_1,EMPC500_1)

Dose_0<-c(rep(50,100),rep(100,100),rep(250,100),rep(500,100))
Dose_1<-c(rep(50,100),rep(100,100),rep(250,100),rep(500,100))

type_0<-rep(0,400)
type_1<-rep(1,400)


#Dose<-c(rep(50,100),rep(100,100),rep(250,100),rep(500,100))
#type<-rbinom(100,1,0.2)

dati_EMPC<-data.frame(EMPC=c(EMPC_0,EMPC_1),Dose=c(Dose_0,Dose_1),type=c(type_0,type_1))
head(dati_EMPC)
summary(dati_EMPC)

###1° MODEL

#page 3

secalonic.LL.4 <- drm(EMPC ~ Dose,
data = dati_EMPC,curveid = type,
fct = LL.4())

summary(secalonic.LL.4)
coeftest(secalonic.LL.4,vcov = sandwich)

confint(secalonic.LL.4)

#plot(secalonic.LL.4,
#bp = 1e-3, broken = TRUE,
#ylim = c(0, 100),
#xlab = "Depth (X)",
#ylab = "(Exon or Gene) Coverage [EMPC-Exon Median Percent Coverage]")
#title("Coverage by Depth")

plot(secalonic.LL.4,
ylim = c(0, 100),
type = "all",
xlab = "Depth (X)",
ylab = "(Exon or Gene) Coverage [EMPC-Exon Median Percent Coverage]")
abline(h=75,lty=2,col=2)
title("Coverage by Depth")




#Estimation of arbitrary ED values
ED(secalonic.LL.4, c(10, 20))


#One way of understanding the calculation is by reversing
#it using the function predict to predict the fitted value corresponding
#to the obtained dose:

predict(secalonic.LL.4,
data.frame(dose = ED(secalonic.LL.4, c(10), display = FALSE)))





###2° MODEL

?drm

barley.LL.4<- drm(EMPC ~ Dose,
data = dati_EMPC,
fct = LL.4(),
na.action = na.omit
)

summary(barley.LL.4.log)
title("Coverage by Depth")

plot(barley.LL.4,
bp = 125,
broken = TRUE,
type = "all",
xlab = "Depth (X)",
ylab = "(Exon or Gene) Coverage [EMPC-Exon Median Percent Coverage]")
title("Coverage by Depth")



#######################################################################################################
#############FARE UNO STUDIO CON INCOGNITE DEPTH E DURATA IN PARAFFINA SU EXON COVERAGE################
#######################################################################################################





