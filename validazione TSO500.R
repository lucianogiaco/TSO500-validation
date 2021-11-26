
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
"Dimension 4\n(from ...)"),vlcex=0.65,axistype=4,axislabcol="black",calcex=0.7)
title("Radar chart of the performance \n(subject level)",cex.main=0.85)



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



colnames(dati)











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


RA_DNA<-dati_new[dati_new$Sample=="RA_DNA",]
#write.table(RA_DNA,"RA_DNA.csv",sep=";",col.names=T,row.names=F,dec=",",na="")




