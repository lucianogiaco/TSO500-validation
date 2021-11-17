
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
#1= SODDISFA LE LINEE GUIDA; 0=NON SODDISFA LE LINEE GUIDA


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

install.packages("vioplot")
library(vioplot)

library(ggplot2)


dati<-read.csv2(file.choose(),sep="\t",h=T,dec=".",na.strings="NA") #importa tab.tab
head(dati)
dim(dati) #63 47
colnames(dati)
summary(dati)




###DNA CONTAMINATION###
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
  labs(title="Plot of length  by dose",x="Dose (mg)", y = "Length")+
  geom_boxplot(width=0.1)+
  theme_classic()

subset(dati, Contamination==0)

table(dati$Contamination==0)
tapply(dati$contamination_score,dati$Contamination,summary)
tapply(dati$contamination_score,dati$Contamination,length)

